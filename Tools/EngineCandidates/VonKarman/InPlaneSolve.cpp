// Smoother, V-cycle and preconditioned conjugate gradients of the in-plane
// solver (see InPlaneMultigrid.cpp for the hierarchy).

#include "InPlaneSolver.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace taikor::candidates::vk
{


double InPlaneSolver::rayleighEstimate (int l, int iterations) const
{
    // Power iteration on D^-1 K from a deterministic random start, finished
    // with its Rayleigh quotient in the D inner product (a lower bound on the
    // largest eigenvalue).
    const InPlaneLevel& L = levels_[static_cast<std::size_t> (l)];
    const auto s = static_cast<std::size_t> (L.size());
    std::vector<double> x (s, 0.0), y (s, 0.0);
    std::uint64_t state = 0x2545f4914f6cdd1dULL + static_cast<std::uint64_t> (l);
    for (int c = 0; c < 2; ++c)
        for (int k = 0; k < L.N; ++k)
        {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            const auto idx = static_cast<std::size_t> (L.base (c) + k);
            x[idx] = L.inverseDiagonal[idx] != 0.0 ? static_cast<double> (state >> 11) * 1.1102230246251565e-16 - 0.5 : 0.0;
        }
    double rayleigh = 0.0;
    for (int it = 0; it < iterations; ++it)
    {
        apply (l, x.data(), y.data());
        double xkx = 0.0, xdx = 0.0;
        for (std::size_t k = 0; k < s; ++k)
            if (L.inverseDiagonal[k] != 0.0)
            {
                xkx += x[k] * y[k];
                xdx += x[k] * x[k] / L.inverseDiagonal[k];
            }
        rayleigh = xdx > 0.0 ? xkx / xdx : 0.0;
        double norm = 0.0;
        for (std::size_t k = 0; k < s; ++k)
        {
            x[k] = y[k] * L.inverseDiagonal[k];
            norm += x[k] * x[k];
        }
        norm = std::sqrt (norm);
        if (norm == 0.0)
            break;
        for (auto& v : x)
            v /= norm;
    }
    return rayleigh;
}

void InPlaneSolver::estimateLambdaMax (int l)
{
    // The Rayleigh estimate after 200 power iterations raised by 10 % (the
    // safety margin of Adams et al. 2003), capped by the Gershgorin bound
    // where that is smaller. tests/vonkarman "solver" checks that the margin
    // covers the estimate after 800 iterations on every level (60 iterations
    // left only 0.3 % to spare).
    InPlaneLevel& L = levels_[static_cast<std::size_t> (l)];
    const double rayleigh = rayleighEstimate (l, 200);
    double gershgorin = 0.0;
    if (! L.matrixFree)
    {
        // Regular rows share the stencil; the rim rows carry their own taps.
        for (int c = 0; c < 2; ++c)
        {
            double sum = 0.0, diag = 0.0;
            for (const auto& t : L.stencil[c])
            {
                sum += std::abs (t.coefficient);
                if (t.offset == 0)
                    diag = t.coefficient;
            }
            if (diag > 0.0)
                gershgorin = std::max (gershgorin, sum / diag);
        }
        // Rim rows: stencil plus their corrections, bounded by the triangle
        // inequality (a larger Gershgorin radius is still a bound).
        for (std::size_t r = 0; r < L.irregularRow.size(); ++r)
        {
            const int c = L.irregularRow[r] >= L.base (1) ? 1 : 0;
            double sum = 0.0;
            for (const auto& t : L.stencil[c])
                sum += std::abs (t.coefficient);
            for (int slot = 0; slot < L.irregularWidth; ++slot)
                sum += std::abs (L.irregularCoefficient[static_cast<std::size_t> (slot) * L.irregularRow.size() + r]);
            gershgorin = std::max (gershgorin, sum * L.inverseDiagonal[static_cast<std::size_t> (L.irregularRow[r])]);
        }
    }
    else
    {
        // Finest operator: |row| <= 2 (2 C11 + 2 C66) + 4 (C12 + C66) terms
        // bounded by the stencil itself; use the exact stencil sum.
        const double diag = 2.0 * (c_.c11 + c_.c66);
        const double off = 2.0 * c_.c11 + 2.0 * c_.c66 + 4.0 * (std::abs (c_.c12) + c_.c66);
        gershgorin = (diag + off) / diag;
    }
    L.lambdaMax = std::min (1.1 * rayleigh, gershgorin);
}

void InPlaneSolver::chebyshev (int l, const float* b, float* x, bool zeroInitial)
{
    // Chebyshev iteration for the Jacobi-scaled system (Saad, Iterative
    // Methods for Sparse Linear Systems, 2nd ed. 2003, Algorithm 12.1), in
    // single precision.
    const InPlaneLevel& L = levels_[static_cast<std::size_t> (l)];
    Work& W = work_[static_cast<std::size_t> (l)];
    const double lmax = L.lambdaMax;
    const double lmin = lmax / 30.0;
    const double theta = 0.5 * (lmax + lmin);
    const double delta = 0.5 * (lmax - lmin);
    const double sigma = theta / delta;
    double rho = 1.0 / sigma;
    const float* __restrict dinv = L.inverseDiagonalF.data();
    float* __restrict r = W.c.data();
    float* __restrict d = W.d.data();
    float* __restrict t = W.t.data();
    float* __restrict xs = x;
    const float* __restrict bs = b;
    const auto invTheta = static_cast<float> (1.0 / theta);
    // The component blocks and the pad between them are contiguous, so one
    // loop covers both (the pad entries stay zero: b, t and dinv are zero there).
    const long o = L.base (0);
    const long count = L.base (1) + L.N - o;
    if (zeroInitial)
    {
        for (long k = o; k < o + count; ++k)
        {
            const float rk = bs[k];
            r[k] = rk;
            const float dk = dinv[k] * rk * invTheta;
            d[k] = dk;
            xs[k] = dk;
        }
    }
    else
    {
        apply (l, x, t);
        for (long k = o; k < o + count; ++k)
        {
            const float rk = bs[k] - t[k];
            r[k] = rk;
            const float dk = dinv[k] * rk * invTheta;
            d[k] = dk;
            xs[k] += dk;
        }
    }
    const int degree = l == 0 ? settings_.chebyshevDegree : settings_.coarseChebyshevDegree;
    for (int step = 1; step < degree; ++step)
    {
        apply (l, d, t);
        const double rhoNew = 1.0 / (2.0 * sigma - rho);
        const auto c1 = static_cast<float> (rhoNew * rho);
        const auto c2 = static_cast<float> (2.0 * rhoNew / delta);
        for (long k = o; k < o + count; ++k)
        {
            const float rk = r[k] - t[k];
            r[k] = rk;
            const float dk = c1 * d[k] + c2 * dinv[k] * rk;
            d[k] = dk;
            xs[k] += dk;
        }
        rho = rhoNew;
    }
}

void InPlaneSolver::vcycle (int l, const float* r, float* z)
{
    const int last = static_cast<int> (levels_.size()) - 1;
    const InPlaneLevel& L = levels_[static_cast<std::size_t> (l)];
    if (l == last)
    {
        // Exact coarsest solve by the stored inverse, column by column.
        const auto m = coarsestIndex_.size();
        float* y = coarsestWork_.data();
        std::fill (y, y + m, 0.0f);
        for (std::size_t j = 0; j < m; ++j)
        {
            const float rj = r[coarsestIndex_[j]];
            const float* __restrict col = coarsestInverse_.data() + j * m;
            for (std::size_t i = 0; i < m; ++i)
                y[i] += col[i] * rj;
        }
        std::fill (z, z + L.size(), 0.0f);
        for (std::size_t i = 0; i < m; ++i)
            z[coarsestIndex_[i]] = y[i];
        return;
    }
    Work& W = work_[static_cast<std::size_t> (l)];
    Work& Wc = work_[static_cast<std::size_t> (l + 1)];
    chebyshev (l, r, z, true);
    float* t = W.t.data();
    apply (l, z, t);
    {
        const long o = L.base (0), e = L.base (1) + L.N;
        for (long k = o; k < e; ++k)
            t[k] = r[k] - t[k];
    }
    restrict_ (l + 1, t, Wc.b.data());
    vcycle (l + 1, Wc.b.data(), Wc.x.data());
    prolong (l + 1, Wc.x.data(), t);
    {
        const long o = L.base (0), e = L.base (1) + L.N;
        for (long k = o; k < e; ++k)
            z[k] += t[k];
    }
    chebyshev (l, r, z, false);
}

void InPlaneSolver::precondition (const double* r, double* z)
{
    const InPlaneLevel& L = levels_.front();
    Work& W = work_.front();
    float* rf = W.b.data();
    float* zf = W.x.data();
    for (int c = 0; c < 2; ++c)
    {
        const long o = L.base (c);
        for (int k = 0; k < L.N; ++k)
            rf[o + k] = static_cast<float> (r[o + k]);
    }
    vcycle (0, rf, zf);
    for (int c = 0; c < 2; ++c)
    {
        const long o = L.base (c);
        for (int k = 0; k < L.N; ++k)
            z[o + k] = static_cast<double> (zf[o + k]);
    }
}

int InPlaneSolver::solve (const double* b, double* x, double& relativeResidual, double tolerance)
{
    // Flexible preconditioned CG, beta = r_new . (z_new - z_old) / (r_old .
    // z_old) (Polak-Ribiere; Y. Notay, "Flexible conjugate gradients", SIAM J.
    // Sci. Comput. 22, 2000): exact CG for a symmetric preconditioner and
    // robust to the single-precision V-cycle's rounding-level asymmetry. The
    // vector updates are fused so an iteration makes four passes besides the
    // operator and the V-cycle; the V-cycle's float input is written by the
    // residual update and its float output is read in place.
    const InPlaneLevel& L = levels_.front();
    Work& W = work_.front();
    const long o = L.base (0);
    const long e = L.base (1) + L.N; // both blocks and the zero pad between them
    double* __restrict r = pcgR_.data();
    double* __restrict p = pcgP_.data();
    double* __restrict q = pcgQ_.data();
    float* __restrict rf = W.b.data();
    float* zf = W.x.data();
    float* zOld = zOldF_.data();
    double bb = 0.0;
    for (long k = o; k < e; ++k)
        bb += b[k] * b[k];
    const double bnorm = std::sqrt (bb);
    if (bnorm == 0.0)
    {
        std::fill (x, x + L.size(), 0.0);
        relativeResidual = 0.0;
        lastInitialResidual_ = 0.0;
        return 0;
    }
    apply (0, x, q);
    double rr = 0.0;
    for (long k = o; k < e; ++k)
    {
        const double rk = b[k] - q[k];
        r[k] = rk;
        rf[k] = static_cast<float> (rk);
        rr += rk * rk;
    }
    const double target = (tolerance > 0.0 ? tolerance : settings_.relativeTolerance) * bnorm;
    double rnorm = std::sqrt (rr);
    lastInitialResidual_ = rnorm / bnorm;
    int it = 0;
    if (rnorm > target)
    {
        vcycle (0, rf, zf);
        double rz = 0.0;
        for (long k = o; k < e; ++k)
        {
            p[k] = zf[k];
            rz += r[k] * zf[k];
        }
        while (it < settings_.maxIterations)
        {
            ++it;
            apply (0, p, q);
            double pq = 0.0;
            for (long k = o; k < e; ++k)
                pq += p[k] * q[k];
            const double alpha = rz / pq;
            rr = 0.0;
            for (long k = o; k < e; ++k)
            {
                x[k] += alpha * p[k];
                const double rk = r[k] - alpha * q[k];
                r[k] = rk;
                rf[k] = static_cast<float> (rk);
                rr += rk * rk;
            }
            rnorm = std::sqrt (rr);
            if (rnorm <= target)
                break;
            std::swap (zf, zOld);
            vcycle (0, rf, zf);
            double rzNew = 0.0, rzOld = 0.0;
            for (long k = o; k < e; ++k)
            {
                rzNew += r[k] * zf[k];
                rzOld += r[k] * zOld[k];
            }
            const double beta = (rzNew - rzOld) / rz;
            for (long k = o; k < e; ++k)
                p[k] = zf[k] + beta * p[k];
            rz = rzNew;
        }
    }
    relativeResidual = rnorm / bnorm;
    return it;
}

// ---- Team versions (Team.h): the same arithmetic per entry as the serial
// functions above, split by rows; sums over the members in member order.

void InPlaneSolver::chebyshev (int l, const float* b, float* x, bool zeroInitial, Team& team)
{
    const InPlaneLevel& L = levels_[static_cast<std::size_t> (l)];
    Work& W = work_[static_cast<std::size_t> (l)];
    const double lmax = L.lambdaMax;
    const double lmin = lmax / 30.0;
    const double theta = 0.5 * (lmax + lmin);
    const double delta = 0.5 * (lmax - lmin);
    const double sigma = theta / delta;
    double rho = 1.0 / sigma;
    const float* __restrict dinv = L.inverseDiagonalF.data();
    float* __restrict r = W.c.data();
    float* __restrict d = W.d.data();
    float* __restrict t = W.t.data();
    float* __restrict xs = x;
    const float* __restrict bs = b;
    const auto invTheta = static_cast<float> (1.0 / theta);
    int j0, j1;
    team.rows (L.n, j0, j1);
    const long span[2][2] = { { L.base (0) + static_cast<long> (j0) * L.n, L.base (0) + static_cast<long> (j1) * L.n },
                              { L.base (1) + static_cast<long> (j0) * L.n, L.base (1) + static_cast<long> (j1) * L.n } };
    if (zeroInitial)
    {
        for (const auto& sp : span)
            for (long k = sp[0]; k < sp[1]; ++k)
            {
                const float rk = bs[k];
                r[k] = rk;
                const float dk = dinv[k] * rk * invTheta;
                d[k] = dk;
                xs[k] = dk;
            }
    }
    else
    {
        applyRows (l, x, t, j0, j1, team.member);
        team.sync(); // every member done reading x before any updates it
        for (const auto& sp : span)
            for (long k = sp[0]; k < sp[1]; ++k)
            {
                const float rk = bs[k] - t[k];
                r[k] = rk;
                const float dk = dinv[k] * rk * invTheta;
                d[k] = dk;
                xs[k] += dk;
            }
    }
    const int degree = l == 0 ? settings_.chebyshevDegree : settings_.coarseChebyshevDegree;
    for (int step = 1; step < degree; ++step)
    {
        team.sync(); // d complete
        applyRows (l, d, t, j0, j1, team.member);
        team.sync(); // every member done reading d before any updates it
        const double rhoNew = 1.0 / (2.0 * sigma - rho);
        const auto c1 = static_cast<float> (rhoNew * rho);
        const auto c2 = static_cast<float> (2.0 * rhoNew / delta);
        for (const auto& sp : span)
            for (long k = sp[0]; k < sp[1]; ++k)
            {
                const float rk = r[k] - t[k];
                r[k] = rk;
                const float dk = c1 * d[k] + c2 * dinv[k] * rk;
                d[k] = dk;
                xs[k] += dk;
            }
        rho = rhoNew;
    }
    team.sync(); // x complete
}

void InPlaneSolver::vcycle (int l, const float* r, float* z, Team& team)
{
    const int last = static_cast<int> (levels_.size()) - 1;
    if (team.size == 1 || l >= teamLevels || l == last)
    {
        if (team.member == 0)
            vcycle (l, r, z);
        team.sync();
        return;
    }
    const InPlaneLevel& L = levels_[static_cast<std::size_t> (l)];
    Work& W = work_[static_cast<std::size_t> (l)];
    Work& Wc = work_[static_cast<std::size_t> (l + 1)];
    int j0, j1;
    team.rows (L.n, j0, j1);
    const long span[2][2] = { { L.base (0) + static_cast<long> (j0) * L.n, L.base (0) + static_cast<long> (j1) * L.n },
                              { L.base (1) + static_cast<long> (j0) * L.n, L.base (1) + static_cast<long> (j1) * L.n } };
    chebyshev (l, r, z, true, team);
    float* t = W.t.data();
    applyRows (l, z, t, j0, j1, team.member);
    for (const auto& sp : span)
        for (long k = sp[0]; k < sp[1]; ++k)
            t[k] = r[k] - t[k];
    team.sync();
    const int nc = levels_[static_cast<std::size_t> (l + 1)].n;
    if (l + 1 < teamLevels && l + 1 != last)
    {
        int c0, c1;
        team.rows (nc, c0, c1);
        restrictRows (l + 1, t, Wc.b.data(), c0, c1);
    }
    else if (team.member == 0)
        restrictRows (l + 1, t, Wc.b.data(), 0, nc);
    team.sync();
    vcycle (l + 1, Wc.b.data(), Wc.x.data(), team); // ends synced
    prolongRows (l + 1, Wc.x.data(), t, j0, j1);
    for (const auto& sp : span)
        for (long k = sp[0]; k < sp[1]; ++k)
            z[k] += t[k];
    team.sync();
    chebyshev (l, r, z, false, team);
}

int InPlaneSolver::solve (const double* b, double* x, double& relativeResidual, double tolerance, Team& team)
{
    // The flexible PCG of the serial solve() above, pass for pass.
    const InPlaneLevel& L = levels_.front();
    Work& W = work_.front();
    int j0, j1;
    team.rows (L.n, j0, j1);
    const long span[2][2] = { { L.base (0) + static_cast<long> (j0) * L.n, L.base (0) + static_cast<long> (j1) * L.n },
                              { L.base (1) + static_cast<long> (j0) * L.n, L.base (1) + static_cast<long> (j1) * L.n } };
    double* __restrict r = pcgR_.data();
    double* __restrict p = pcgP_.data();
    double* __restrict q = pcgQ_.data();
    float* __restrict rf = W.b.data();
    float* zf = W.x.data();
    float* zOld = zOldF_.data();
    applyRows (0, x, q, j0, j1, team.member);
    double bbLocal = 0.0, rrLocal = 0.0;
    for (const auto& sp : span)
        for (long k = sp[0]; k < sp[1]; ++k)
        {
            bbLocal += b[k] * b[k];
            const double rk = b[k] - q[k];
            r[k] = rk;
            rf[k] = static_cast<float> (rk);
            rrLocal += rk * rk;
        }
    const auto sums = team.sum2 (bbLocal, rrLocal);
    const double bnorm = std::sqrt (sums[0]);
    if (bnorm == 0.0)
    {
        for (const auto& sp : span)
            std::fill (x + sp[0], x + sp[1], 0.0);
        if (team.member == 0)
            lastInitialResidual_ = 0.0;
        team.sync();
        relativeResidual = 0.0;
        return 0;
    }
    const double target = (tolerance > 0.0 ? tolerance : settings_.relativeTolerance) * bnorm;
    double rnorm = std::sqrt (sums[1]);
    if (team.member == 0)
        lastInitialResidual_ = rnorm / bnorm;
    int it = 0;
    if (rnorm > target)
    {
        vcycle (0, rf, zf, team);
        double rzLocal = 0.0;
        for (const auto& sp : span)
            for (long k = sp[0]; k < sp[1]; ++k)
            {
                p[k] = zf[k];
                rzLocal += r[k] * zf[k];
            }
        double rz = team.sum (rzLocal); // p complete
        while (it < settings_.maxIterations)
        {
            ++it;
            applyRows (0, p, q, j0, j1, team.member);
            double pqLocal = 0.0;
            for (const auto& sp : span)
                for (long k = sp[0]; k < sp[1]; ++k)
                    pqLocal += p[k] * q[k];
            const double alpha = rz / team.sum (pqLocal);
            rrLocal = 0.0;
            for (const auto& sp : span)
                for (long k = sp[0]; k < sp[1]; ++k)
                {
                    x[k] += alpha * p[k];
                    const double rk = r[k] - alpha * q[k];
                    r[k] = rk;
                    rf[k] = static_cast<float> (rk);
                    rrLocal += rk * rk;
                }
            rnorm = std::sqrt (team.sum (rrLocal)); // x and rf complete
            if (rnorm <= target)
                break;
            std::swap (zf, zOld);
            vcycle (0, rf, zf, team);
            double rzNewLocal = 0.0, rzOldLocal = 0.0;
            for (const auto& sp : span)
                for (long k = sp[0]; k < sp[1]; ++k)
                {
                    rzNewLocal += r[k] * zf[k];
                    rzOldLocal += r[k] * zOld[k];
                }
            const auto rzs = team.sum2 (rzNewLocal, rzOldLocal);
            const double beta = (rzs[0] - rzs[1]) / rz;
            for (const auto& sp : span)
                for (long k = sp[0]; k < sp[1]; ++k)
                    p[k] = zf[k] + beta * p[k];
            rz = rzs[0];
            team.sync(); // p complete
        }
    }
    relativeResidual = rnorm / bnorm;
    return it;
}

} // namespace taikor::candidates::vk
