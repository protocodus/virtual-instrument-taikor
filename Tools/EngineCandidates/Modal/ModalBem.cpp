// Assembly and solution of the axisymmetric exterior problem (ModalBem.h).

#include "ModalBem.h"
#include "ModalInternal.h"

#include <algorithm>
#include <cmath>
#include <map>

namespace taikor::candidates::modal
{
namespace
{
using cplx = std::complex<double>;

struct Rule
{
    std::vector<double> x, w;
};

const Rule& rule (int n)
{
    thread_local std::map<int, Rule> rules;
    auto& r = rules[n];
    if (r.x.empty())
        special::gaussLegendre (n, r.x, r.w);
    return r;
}

// Nodes t and weights along a panel for a field point (rx, zx); self = the
// point is the panel's own midpoint.
void panelNodes (const BemPanel& p, double rx, double zx, bool self, int refine, std::vector<double>& t,
                 std::vector<double>& w)
{
    t.clear();
    w.clear();
    const double len = p.length;
    if (self)
    {
        // Log singularity at the midpoint: t = len/2 +- (len/2) u^3 on each half.
        const Rule& g = rule (12 * refine);
        for (int side = -1; side <= 1; side += 2)
            for (std::size_t i = 0; i < g.x.size(); ++i)
            {
                const double u = 0.5 * (g.x[i] + 1.0);
                t.push_back (0.5 * len + side * 0.5 * len * u * u * u);
                w.push_back (0.5 * g.w[i] * 0.5 * len * 3.0 * u * u);
            }
        return;
    }
    const double along = std::clamp ((rx - p.r0) * p.tr + (zx - p.z0) * p.tz, 0.0, len);
    const double dist = std::hypot (rx - p.r (along), zx - p.z (along));
    if (dist < 1.5 * len)
    {
        // Near: sinh grading towards the nearest point on each side of it.
        const Rule& g = rule (8 * refine);
        const double d = std::max (dist, 1.0e-6 * len);
        for (int side = -1; side <= 1; side += 2)
        {
            const double piece = side < 0 ? along : len - along;
            if (piece <= 0.0)
                continue;
            const double uMax = std::asinh (piece / d);
            for (std::size_t i = 0; i < g.x.size(); ++i)
            {
                const double u = 0.5 * uMax * (g.x[i] + 1.0);
                t.push_back (along + side * d * std::sinh (u));
                w.push_back (0.5 * uMax * g.w[i] * d * std::cosh (u));
            }
        }
        return;
    }
    const Rule& g = rule (4 * refine);
    for (std::size_t i = 0; i < g.x.size(); ++i)
    {
        t.push_back (0.5 * len * (g.x[i] + 1.0));
        w.push_back (0.5 * len * g.w[i]);
    }
}

// Solves min |A x - B| for A (rows x n, column-major, rows >= n) and B (rows x
// nrhs, column-major) by Householder QR; X (n x nrhs) is left in the first n
// rows of each column of B. Returns false if A is rank deficient.
bool leastSquares (int rows, int n, std::vector<cplx>& a, std::vector<cplx>& b, int nrhs)
{
    const auto at = [&] (int i, int j) -> cplx& { return a[static_cast<std::size_t> (j) * rows + i]; };
    const auto bt = [&] (int i, int j) -> cplx& { return b[static_cast<std::size_t> (j) * rows + i]; };
    std::vector<cplx> v (static_cast<std::size_t> (rows));
    for (int k = 0; k < n; ++k)
    {
        double norm = 0.0;
        for (int i = k; i < rows; ++i)
            norm += std::norm (at (i, k));
        norm = std::sqrt (norm);
        if (norm == 0.0)
            return false;
        const cplx x0 = at (k, k);
        const cplx alpha = -(std::abs (x0) > 0.0 ? x0 / std::abs (x0) : cplx (1.0)) * norm;
        for (int i = k; i < rows; ++i)
            v[static_cast<std::size_t> (i)] = at (i, k);
        v[static_cast<std::size_t> (k)] -= alpha;
        double vn = 0.0;
        for (int i = k; i < rows; ++i)
            vn += std::norm (v[static_cast<std::size_t> (i)]);
        vn = std::sqrt (vn);
        if (vn == 0.0)
            continue;
        for (int i = k; i < rows; ++i)
            v[static_cast<std::size_t> (i)] /= vn;
        // H = I - 2 v v^H on the remaining columns of A and on B.
        const auto reflect = [&] (auto&& column)
        {
            cplx dot = 0.0;
            for (int i = k; i < rows; ++i)
                dot += std::conj (v[static_cast<std::size_t> (i)]) * column (i);
            dot *= 2.0;
            for (int i = k; i < rows; ++i)
                column (i) -= v[static_cast<std::size_t> (i)] * dot;
        };
        for (int j = k; j < n; ++j)
            reflect ([&] (int i) -> cplx& { return at (i, j); });
        for (int j = 0; j < nrhs; ++j)
            reflect ([&] (int i) -> cplx& { return bt (i, j); });
    }
    for (int j = 0; j < nrhs; ++j)
        for (int i = n - 1; i >= 0; --i)
        {
            cplx sum = bt (i, j);
            for (int l = i + 1; l < n; ++l)
                sum -= at (i, l) * bt (l, j);
            if (std::abs (at (i, i)) == 0.0)
                return false;
            bt (i, j) = sum / at (i, i);
        }
    return true;
}
// Solves the square system A X = B (A n x n, B n x nrhs, both column-major) by
// LU with partial pivoting; X replaces B.
bool luSolve (int n, std::vector<cplx>& a, std::vector<cplx>& b, int nrhs)
{
    const auto at = [&] (int i, int j) -> cplx& { return a[static_cast<std::size_t> (j) * n + i]; };
    const auto bt = [&] (int i, int j) -> cplx& { return b[static_cast<std::size_t> (j) * n + i]; };
    for (int k = 0; k < n; ++k)
    {
        int pivot = k;
        for (int i = k + 1; i < n; ++i)
            if (std::abs (at (i, k)) > std::abs (at (pivot, k)))
                pivot = i;
        if (std::abs (at (pivot, k)) == 0.0)
            return false;
        if (pivot != k)
        {
            for (int j = 0; j < n; ++j)
                std::swap (at (k, j), at (pivot, j));
            for (int j = 0; j < nrhs; ++j)
                std::swap (bt (k, j), bt (pivot, j));
        }
        const cplx inv = 1.0 / at (k, k);
        for (int i = k + 1; i < n; ++i)
            at (i, k) *= inv;
        for (int j = k + 1; j < n; ++j)
        {
            const cplx f = at (k, j);
            if (f == 0.0)
                continue;
            cplx* col = &at (0, j);
            const cplx* piv = &at (0, k);
            for (int i = k + 1; i < n; ++i)
                col[i] -= piv[i] * f;
        }
        for (int j = 0; j < nrhs; ++j)
        {
            const cplx f = bt (k, j);
            if (f == 0.0)
                continue;
            cplx* col = &bt (0, j);
            const cplx* piv = &at (0, k);
            for (int i = k + 1; i < n; ++i)
                col[i] -= piv[i] * f;
        }
    }
    for (int j = 0; j < nrhs; ++j)
        for (int i = n - 1; i >= 0; --i)
        {
            cplx sum = bt (i, j);
            for (int l = i + 1; l < n; ++l)
                sum -= at (i, l) * bt (l, j);
            bt (i, j) = sum / at (i, i);
        }
    return true;
}
} // namespace

BemSolver::BemSolver (BemGeometry geometry, int refineIn) : geo (std::move (geometry)), refine (std::max (1, refineIn)) {}

void BemSolver::solve (double k, double density, double soundSpeed, std::vector<BemProblem>& problems, bool useChief,
                       int threads) const
{
    const int n = static_cast<int> (geo.panels.size());
    int chiefCount = 0;
    for (const auto& pr : problems)
        if (useChief && pr.chief)
            chiefCount = static_cast<int> (geo.chief.size());
    // Field points are those of the first problem; every problem gets them.
    int maxOrder = 0;
    for (const auto& pr : problems)
        maxOrder = std::max (maxOrder, pr.m);
    const std::vector<std::array<double, 2>> noPoints;
    const auto& fieldPoints = problems.empty() ? noPoints : problems.front().fieldPoints;
    const int fieldCount = static_cast<int> (fieldPoints.size());
    const int rows = n + chiefCount;
    const int allRows = rows + fieldCount;
    const auto nProblems = problems.size();
    // Active index of each panel per problem (-1 if its part carries no velocity).
    std::vector<std::vector<int>> active (nProblems, std::vector<int> (static_cast<std::size_t> (n), -1));
    for (std::size_t q = 0; q < nProblems; ++q)
    {
        int index = 0;
        for (int j = 0; j < n; ++j)
            if (std::find (problems[q].activeParts.begin(), problems[q].activeParts.end(), geo.panels[static_cast<std::size_t> (j)].part)
                != problems[q].activeParts.end())
                active[q][static_cast<std::size_t> (j)] = index++;
    }
    // Per problem: A (allRows x n, column-major by row block) and single-layer RHS (allRows x columns).
    std::vector<std::vector<cplx>> amat (nProblems), rhs (nProblems);
    for (std::size_t q = 0; q < nProblems; ++q)
    {
        amat[q].assign (static_cast<std::size_t> (allRows) * static_cast<std::size_t> (n), 0.0);
        rhs[q].assign (static_cast<std::size_t> (allRows) * static_cast<std::size_t> (problems[q].columns), 0.0);
    }
    parallelFor (allRows, threads, [&] (int row)
    {
        double rx = 0.0, zx = 0.0;
        int selfPanel = -1;
        if (row < n)
        {
            const auto& p = geo.panels[static_cast<std::size_t> (row)];
            rx = p.r (0.5 * p.length);
            zx = p.z (0.5 * p.length);
            selfPanel = row;
        }
        else if (row < rows)
        {
            rx = geo.chief[static_cast<std::size_t> (row - n)][0];
            zx = geo.chief[static_cast<std::size_t> (row - n)][1];
        }
        else
        {
            rx = fieldPoints[static_cast<std::size_t> (row - rows)][0];
            zx = fieldPoints[static_cast<std::size_t> (row - rows)][1];
        }
        std::vector<double> t, w;
        std::vector<cplx> g (static_cast<std::size_t> (maxOrder + 1)), h (g.size());
        std::vector<cplx> aSum (g.size()), sSum (g.size() * bemVelocityNodes);
        for (int j = 0; j < n; ++j)
        {
            const BemPanel& p = geo.panels[static_cast<std::size_t> (j)];
            panelNodes (p, rx, zx, j == selfPanel, refine, t, w);
            const auto nodes = bemNodeOffsets (p);
            std::fill (aSum.begin(), aSum.end(), 0.0);
            std::fill (sSum.begin(), sSum.end(), 0.0);
            bool isActive = false;
            for (std::size_t q = 0; q < nProblems; ++q)
                isActive = isActive || active[q][static_cast<std::size_t> (j)] >= 0;
            for (std::size_t i = 0; i < t.size(); ++i)
            {
                const double ry = p.r (t[i]), zy = p.z (t[i]);
                ringKernels (rx, zx, ry, zy, p.nr, p.nz, k, maxOrder, refine, g.data(), h.data());
                const double wr = w[i] * ry;
                for (int m = 0; m <= maxOrder; ++m)
                    aSum[static_cast<std::size_t> (m)] += wr * h[static_cast<std::size_t> (m)];
                if (! isActive)
                    continue;
                for (int qn = 0; qn < bemVelocityNodes; ++qn)
                {
                    double lagrange = 1.0;
                    for (int pn = 0; pn < bemVelocityNodes; ++pn)
                        if (pn != qn)
                            lagrange *= (t[i] - nodes[static_cast<std::size_t> (pn)])
                                      / (nodes[static_cast<std::size_t> (qn)] - nodes[static_cast<std::size_t> (pn)]);
                    for (int m = 0; m <= maxOrder; ++m)
                        sSum[static_cast<std::size_t> (m * bemVelocityNodes + qn)] += (wr * lagrange) * g[static_cast<std::size_t> (m)];
                }
            }
            for (std::size_t q = 0; q < nProblems; ++q)
            {
                const auto& pr = problems[q];
                const auto um = static_cast<std::size_t> (pr.m);
                amat[q][static_cast<std::size_t> (j) * static_cast<std::size_t> (allRows) + static_cast<std::size_t> (row)] = aSum[um];
                const int ai = active[q][static_cast<std::size_t> (j)];
                if (ai < 0)
                    continue;
                for (int qn = 0; qn < bemVelocityNodes; ++qn)
                {
                    const cplx s = sSum[um * bemVelocityNodes + static_cast<std::size_t> (qn)];
                    const double* vel = &pr.velocity[static_cast<std::size_t> (ai * bemVelocityNodes + qn) * static_cast<std::size_t> (pr.columns)];
                    for (int c = 0; c < pr.columns; ++c)
                        rhs[q][static_cast<std::size_t> (c) * static_cast<std::size_t> (allRows) + static_cast<std::size_t> (row)] += s * vel[c];
                }
            }
        }
    });
    const cplx jwr (0.0, k * soundSpeed * density);
    parallelFor (static_cast<int> (nProblems), threads, [&] (int qi)
    {
        const auto q = static_cast<std::size_t> (qi);
        auto& pr = problems[q];
        const int cols = pr.columns;
        const int used = pr.chief && useChief ? rows : n; // rows of this problem's system
        std::vector<cplx> lhs (static_cast<std::size_t> (used) * static_cast<std::size_t> (n));
        std::vector<cplx> b (static_cast<std::size_t> (used) * static_cast<std::size_t> (cols));
        // CHIEF rows enter with their natural weight: for high orders the interior
        // field, and so the row, is small, and scaling it up would amplify the
        // quadrature error into spurious constraints.
        for (int i = 0; i < used; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                cplx value = -amat[q][static_cast<std::size_t> (j) * static_cast<std::size_t> (allRows) + static_cast<std::size_t> (i)];
                if (i == j)
                    value += 0.5;
                lhs[static_cast<std::size_t> (j) * static_cast<std::size_t> (used) + static_cast<std::size_t> (i)] = value;
            }
            for (int c = 0; c < cols; ++c)
                b[static_cast<std::size_t> (c) * static_cast<std::size_t> (used) + static_cast<std::size_t> (i)]
                    = jwr * rhs[q][static_cast<std::size_t> (c) * static_cast<std::size_t> (allRows) + static_cast<std::size_t> (i)];
        }
        if (used > n && useChief)
            leastSquares (used, n, lhs, b, cols);
        else
            luSolve (n, lhs, b, cols);
        pr.pressure.assign (static_cast<std::size_t> (n) * static_cast<std::size_t> (cols), 0.0);
        for (int j = 0; j < n; ++j)
            for (int c = 0; c < cols; ++c)
                pr.pressure[static_cast<std::size_t> (j) * static_cast<std::size_t> (cols) + static_cast<std::size_t> (c)]
                    = b[static_cast<std::size_t> (c) * static_cast<std::size_t> (used) + static_cast<std::size_t> (j)];
        // Field points: p(x) = sum_j A_j(x) p_j + j w rho S(x) v.
        pr.field.assign (fieldPoints.size() * static_cast<std::size_t> (cols), 0.0);
        for (std::size_t f = 0; f < fieldPoints.size(); ++f)
        {
            const auto row = static_cast<std::size_t> (rows) + f;
            for (int c = 0; c < cols; ++c)
            {
                cplx sum = jwr * rhs[q][static_cast<std::size_t> (c) * static_cast<std::size_t> (allRows) + row];
                for (int j = 0; j < n; ++j)
                    sum += amat[q][static_cast<std::size_t> (j) * static_cast<std::size_t> (allRows) + row]
                         * pr.pressure[static_cast<std::size_t> (j) * static_cast<std::size_t> (cols) + static_cast<std::size_t> (c)];
                pr.field[f * static_cast<std::size_t> (cols) + static_cast<std::size_t> (c)] = sum;
            }
        }
    });
}

} // namespace taikor::candidates::modal
