// Multigrid-preconditioned conjugate gradients for the in-plane problem
// K U = b of InPlaneSolver.h.
//
// Hierarchy: MAC-type 2:1 coarsening with the bilinear staggered prolongation
// P of InPlaneOperator.cpp, masked to the free dofs, and Galerkin coarse
// operators K_c = P^T K P. The Galerkin product follows the staircase rim
// exactly on every level, which rediscretisation on a coarser staircase would
// not: with rediscretised coarse levels the preconditioned CG converges at
// 0.59 per iteration instead of 0.24 (tests/vonkarman "solver"). It is
// assembled by probing: the coarse stencil couples dofs at most two coarse
// cells apart, so applying P^T K P to the indicator of every fifth dof in
// each direction (25 colours per component) recovers every coefficient
// without overlap; build() verifies the assembled operator against the
// product on random vectors. It is stored as the 31-tap stencil of the
// interior rows plus sparse corrections on the rows near the rim.
//
// Smoother: Chebyshev iteration on the Jacobi-scaled operator over
// [lambda_max / 30, lambda_max] (Adams, Brezina, Hu & Tuminaro, "Parallel
// multigrid smoothing: polynomial versus Gauss-Seidel", J. Comput. Phys. 188,
// 2003, whose interval this is). A polynomial smoother is its own adjoint,
// so the V-cycle with equal pre- and post-smoothing and an exact coarsest
// solve is a symmetric positive definite preconditioner for CG. Degree 3 on
// the finest level and 1 below cost least per decade of residual among the
// degrees 1-4 x 1-3 measured on the 48 kHz grid (about 0.17 per iteration).

#include "InPlaneSolver.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <tuple>

namespace taikor::candidates::vk
{

namespace
{
// Deterministic uniform numbers in [-1, 1) for the power iteration and the
// probe check (a 64-bit linear congruential generator, Knuth MMIX constants).
struct Lcg
{
    std::uint64_t state = 0x9e3779b97f4a7c15ULL;
    double next() noexcept
    {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double> (state >> 11) * (2.0 / 9007199254740992.0) - 1.0;
    }
};

void randomFree (const InPlaneLevel& L, double* x, Lcg& rng)
{
    std::fill (x, x + L.size(), 0.0);
    for (int k = 0; k < L.N; ++k)
    {
        if (L.freeU[static_cast<std::size_t> (k)] != 0)
            x[L.base (0) + k] = rng.next();
        if (L.freeV[static_cast<std::size_t> (k)] != 0)
            x[L.base (1) + k] = rng.next();
    }
}

long countFree (const InPlaneLevel& L)
{
    long c = 0;
    for (int k = 0; k < L.N; ++k)
        c += (L.freeU[static_cast<std::size_t> (k)] != 0) + (L.freeV[static_cast<std::size_t> (k)] != 0);
    return c;
}
} // namespace

bool InPlaneSolver::build (int n, const std::vector<std::uint8_t>& interior,
                           const PlaneStress& stiffness, const Settings& settings,
                           std::string& error)
{
    c_ = stiffness;
    settings_ = settings;
    levels_.clear();
    work_.clear();
    if (n < 8 || static_cast<long> (interior.size()) != static_cast<long> (n) * n)
    {
        error = "von Karman in-plane solver: bad grid";
        return false;
    }

    InPlaneLevel fineLevel;
    fineLevel.n = n;
    fineLevel.N = n * n;
    fineLevel.pad = 3 * n + 8; // restriction reads rows -3..n of a block
    fineLevel.freeU.assign (static_cast<std::size_t> (fineLevel.N), 0);
    fineLevel.freeV.assign (static_cast<std::size_t> (fineLevel.N), 0);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
        {
            const int k = j * n + i;
            if (interior[static_cast<std::size_t> (k)] == 0)
                continue;
            if (i == 0 || j == 0 || i == n - 1 || j == n - 1)
            {
                error = "von Karman in-plane solver: interior node on the grid border";
                return false;
            }
            fineLevel.freeU[static_cast<std::size_t> (k)] = interior[static_cast<std::size_t> (k + 1)];
            fineLevel.freeV[static_cast<std::size_t> (k)] = interior[static_cast<std::size_t> (k + n)];
        }
    // Every free dof of the finest operator has the same diagonal: two normal
    // strains (C11 each) and two shear corners (C66 each).
    fineLevel.inverseDiagonal.assign (static_cast<std::size_t> (fineLevel.size()), 0.0);
    const double diag = 2.0 * (c_.c11 + c_.c66);
    for (int k = 0; k < fineLevel.N; ++k)
    {
        if (fineLevel.freeU[static_cast<std::size_t> (k)] != 0)
            fineLevel.inverseDiagonal[static_cast<std::size_t> (fineLevel.base (0) + k)] = 1.0 / diag;
        if (fineLevel.freeV[static_cast<std::size_t> (k)] != 0)
            fineLevel.inverseDiagonal[static_cast<std::size_t> (fineLevel.base (1) + k)] = 1.0 / diag;
    }
    fineLevel.matrixFree = true;
    levels_.push_back (std::move (fineLevel));
    addWork (levels_[0]);
    estimateLambdaMax (0);

    while (countFree (levels_.back()) > settings_.coarsestMaxFree)
    {
        if (! buildCoarse (static_cast<int> (levels_.size()) - 1, error))
            return false;
        estimateLambdaMax (static_cast<int> (levels_.size()) - 1);
    }
    if (! buildCoarsest (error))
        return false;
    for (auto& L : levels_)
        L.inverseDiagonalF.assign (L.inverseDiagonal.begin(), L.inverseDiagonal.end());

    const auto s = static_cast<std::size_t> (levels_[0].size());
    for (auto* v : { &pcgR_, &pcgZ_, &pcgP_, &pcgQ_ })
        v->assign (s, 0.0);
    zOldF_.assign (s, 0.0f);
    return true;
}

void InPlaneSolver::addWork (const InPlaneLevel& L)
{
    Work w;
    const auto s = static_cast<std::size_t> (L.size());
    for (auto* v : { &w.b, &w.x, &w.t, &w.d, &w.c, &w.syyF, &w.sxyF })
        v->assign (s, 0.0f);
    // The matrix-free apply's six rolling rows, one set per team member (Team.h).
    const auto rolling = std::max (s, static_cast<std::size_t> (maxTeam) * 6 * static_cast<std::size_t> (L.n));
    w.sxxF.assign (rolling, 0.0f);
    if (L.matrixFree)
    {
        w.sxx.assign (rolling, 0.0);
        for (auto* v : { &w.syy, &w.sxy })
            v->assign (s, 0.0);
    }
    work_.push_back (std::move (w));
}

bool InPlaneSolver::buildCoarse (int fl, std::string& error)
{
    InPlaneLevel C;
    {
        const InPlaneLevel& F = levels_[static_cast<std::size_t> (fl)];
        C.n = F.n / 2 + 3;
        C.N = C.n * C.n;
        C.pad = 3 * C.n + 8;
        markCoarseFree (F, C);
    }
    C.matrixFree = ! settings_.galerkinCoarse;
    levels_.push_back (std::move (C));
    addWork (levels_.back());
    const int cl = fl + 1;
    const InPlaneLevel& F = levels_[static_cast<std::size_t> (fl)];
    InPlaneLevel& L = levels_[static_cast<std::size_t> (cl)];
    const int nc = L.n;

    if (! settings_.galerkinCoarse)
    {
        // Rediscretisation: the same staggered operator at spacing 2h on the
        // coarse free set (an SPD operator on any free set).
        L.inverseDiagonal.assign (static_cast<std::size_t> (L.size()), 0.0);
        for (int k = 0; k < L.N; ++k)
        {
            if (L.freeU[static_cast<std::size_t> (k)] != 0)
                L.inverseDiagonal[static_cast<std::size_t> (L.base (0) + k)] = 1.0 / (2.0 * (c_.c11 + c_.c66));
            if (L.freeV[static_cast<std::size_t> (k)] != 0)
                L.inverseDiagonal[static_cast<std::size_t> (L.base (1) + k)] = 1.0 / (2.0 * (c_.c11 + c_.c66));
        }
        return true;
    }

    constexpr int period = 5; // > twice the coarse stencil radius (2)
    std::vector<double> e (static_cast<std::size_t> (L.size()));
    std::vector<double> tf (static_cast<std::size_t> (F.size()));
    std::vector<double> yf (static_cast<std::size_t> (F.size()));
    std::vector<double> yc (static_cast<std::size_t> (L.size()));
    std::map<std::tuple<int, int, int, int>, std::vector<double>> bandMap;
    auto isFree = [&L] (int comp, int k)
    {
        return (comp == 0 ? L.freeU : L.freeV)[static_cast<std::size_t> (k)] != 0;
    };

    for (int colComp = 0; colComp < 2; ++colComp)
        for (int cb = 0; cb < period; ++cb)
            for (int ca = 0; ca < period; ++ca)
            {
                std::fill (e.begin(), e.end(), 0.0);
                bool any = false;
                for (int J = cb; J < nc; J += period)
                    for (int I = ca; I < nc; I += period)
                        if (isFree (colComp, J * nc + I))
                        {
                            e[static_cast<std::size_t> (L.base (colComp) + J * nc + I)] = 1.0;
                            any = true;
                        }
                if (! any)
                    continue;
                prolong (cl, e.data(), tf.data());
                apply (fl, tf.data(), yf.data());
                restrict_ (cl, yf.data(), yc.data());
                for (int rowComp = 0; rowComp < 2; ++rowComp)
                    for (int J = 0; J < nc; ++J)
                        for (int I = 0; I < nc; ++I)
                        {
                            const int k = J * nc + I;
                            const double value = yc[static_cast<std::size_t> (L.base (rowComp) + k)];
                            if (value == 0.0 || ! isFree (rowComp, k))
                                continue;
                            // The unique column of this colour within +-2 of the row.
                            const int dI = ((ca - I) % period + period + 2) % period - 2;
                            const int dJ = ((cb - J) % period + period + 2) % period - 2;
                            const int Ic = I + dI, Jc = J + dJ;
                            if (Ic < 0 || Jc < 0 || Ic >= nc || Jc >= nc || ! isFree (colComp, Jc * nc + Ic))
                            {
                                error = "von Karman in-plane solver: coarse stencil wider than the probe period";
                                return false;
                            }
                            auto& coef = bandMap[{ rowComp, colComp, dI, dJ }];
                            if (coef.empty())
                                coef.assign (static_cast<std::size_t> (L.N), 0.0);
                            coef[static_cast<std::size_t> (k)] = value;
                        }
            }

    // Split the operator into the constant stencil of the interior rows and
    // sparse corrections on the rows near the rim. The reference row is the
    // grid centre; a coefficient counts as equal to it within 1e-12 of the
    // diagonal (they agree to rounding in the interior, where P^T K P is
    // translation invariant).
    L.inverseDiagonal.assign (static_cast<std::size_t> (L.size()), 0.0);
    for (auto& [key, coef] : bandMap)
    {
        const auto [rowComp, colComp, dI, dJ] = key;
        if (rowComp == colComp && dI == 0 && dJ == 0)
            for (int k = 0; k < L.N; ++k)
                if (coef[static_cast<std::size_t> (k)] != 0.0)
                    L.inverseDiagonal[static_cast<std::size_t> (L.base (rowComp) + k)] = 1.0 / coef[static_cast<std::size_t> (k)];
    }
    const int centre = (nc / 2) * nc + nc / 2;
    if (! isFree (0, centre) || ! isFree (1, centre))
    {
        error = "von Karman in-plane solver: coarse grid centre is not free";
        return false;
    }
    L.irregularRow.clear();
    std::vector<std::vector<InPlaneLevel::Tap>> rimTaps;
    for (int rowComp = 0; rowComp < 2; ++rowComp)
    {
        L.stencil[rowComp].clear();
        double diagonal = 0.0;
        for (auto& [key, coef] : bandMap)
            if (std::get<0> (key) == rowComp && std::get<1> (key) == rowComp && std::get<2> (key) == 0 && std::get<3> (key) == 0)
                diagonal = std::abs (coef[static_cast<std::size_t> (centre)]);
        struct Band { long offset; int colComp, dI, dJ; const std::vector<double>* coef; };
        std::vector<Band> bands;
        for (auto& [key, coef] : bandMap)
        {
            const auto [rc, colComp, dI, dJ] = key;
            if (rc != rowComp)
                continue;
            const long offset = L.base (colComp) - L.base (rowComp) + static_cast<long> (dJ) * nc + dI;
            bands.push_back ({ offset, colComp, dI, dJ, &coef });
            if (coef[static_cast<std::size_t> (centre)] != 0.0)
                L.stencil[rowComp].push_back ({ offset, coef[static_cast<std::size_t> (centre)] });
        }
        // Rim rows keep only their differences from the stencil, on free
        // columns (x is zero on the others).
        for (int k = 0; k < L.N; ++k)
        {
            if (! isFree (rowComp, k))
                continue;
            std::vector<InPlaneLevel::Tap> taps;
            for (const auto& b : bands)
            {
                const double delta = (*b.coef)[static_cast<std::size_t> (k)] - (*b.coef)[static_cast<std::size_t> (centre)];
                const int column = k + b.dJ * nc + b.dI;
                if (std::abs (delta) > 1e-12 * diagonal && column >= 0 && column < L.N && isFree (b.colComp, column))
                    taps.push_back ({ b.offset, delta });
            }
            if (taps.empty())
                continue;
            L.irregularRow.push_back (L.base (rowComp) + k);
            rimTaps.push_back (std::move (taps));
        }
        L.stencilF[rowComp].clear();
        for (const auto& t : L.stencil[rowComp])
            L.stencilF[rowComp].push_back (static_cast<float> (t.coefficient));
    }
    const std::size_t rows = L.irregularRow.size();
    L.irregularWidth = 0;
    for (const auto& taps : rimTaps)
        L.irregularWidth = std::max (L.irregularWidth, static_cast<int> (taps.size()));
    L.irregularOffset.assign (rows * static_cast<std::size_t> (L.irregularWidth), 0);
    L.irregularCoefficient.assign (L.irregularOffset.size(), 0.0);
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t t = 0; t < rimTaps[r].size(); ++t)
        {
            L.irregularOffset[t * rows + r] = rimTaps[r][t].offset;
            L.irregularCoefficient[t * rows + r] = rimTaps[r][t].coefficient;
        }
    L.irregularCoefficientF.assign (L.irregularCoefficient.begin(), L.irregularCoefficient.end());
    work_[static_cast<std::size_t> (cl)].rimF.assign (rows, 0.0f);
    work_[static_cast<std::size_t> (cl)].rim.assign (rows, 0.0);

    // Verify the probed operator against P^T K P on random vectors: a stencil
    // wider than the probe period would alias silently otherwise.
    Lcg rng;
    std::vector<double> x (static_cast<std::size_t> (L.size()));
    std::vector<double> yb (static_cast<std::size_t> (L.size()));
    for (int trial = 0; trial < 2; ++trial)
    {
        randomFree (L, x.data(), rng);
        prolong (cl, x.data(), tf.data());
        apply (fl, tf.data(), yf.data());
        restrict_ (cl, yf.data(), yc.data());
        apply (cl, x.data(), yb.data());
        double num = 0.0, den = 0.0;
        for (std::size_t k = 0; k < yc.size(); ++k)
        {
            num += (yc[k] - yb[k]) * (yc[k] - yb[k]);
            den += yc[k] * yc[k];
        }
        if (! (num <= 1.0e-24 * den))
        {
            error = "von Karman in-plane solver: probed coarse operator does not match P^T K P";
            return false;
        }
    }
    return true;
}

bool InPlaneSolver::buildCoarsest (std::string& error)
{
    const int l = static_cast<int> (levels_.size()) - 1;
    const InPlaneLevel& L = levels_[static_cast<std::size_t> (l)];
    coarsestIndex_.clear();
    for (int comp = 0; comp < 2; ++comp)
        for (int k = 0; k < L.N; ++k)
            if ((comp == 0 ? L.freeU : L.freeV)[static_cast<std::size_t> (k)] != 0)
                coarsestIndex_.push_back (L.base (comp) + k);
    const auto m = coarsestIndex_.size();
    coarsestFactor_.assign (m * m, 0.0);
    coarsestRhs_.assign (m, 0.0);
    std::vector<double> e (static_cast<std::size_t> (L.size()), 0.0);
    std::vector<double> y (static_cast<std::size_t> (L.size()), 0.0);
    for (std::size_t c = 0; c < m; ++c)
    {
        e[static_cast<std::size_t> (coarsestIndex_[c])] = 1.0;
        apply (l, e.data(), y.data());
        e[static_cast<std::size_t> (coarsestIndex_[c])] = 0.0;
        for (std::size_t r = 0; r < m; ++r)
            coarsestFactor_[r * m + c] = y[static_cast<std::size_t> (coarsestIndex_[r])];
    }
    // In-place Cholesky, lower triangle.
    double* A = coarsestFactor_.data();
    for (std::size_t j = 0; j < m; ++j)
    {
        double d = A[j * m + j];
        for (std::size_t k = 0; k < j; ++k)
            d -= A[j * m + k] * A[j * m + k];
        if (! (d > 0.0))
        {
            error = "von Karman in-plane solver: coarsest operator not positive definite";
            return false;
        }
        d = std::sqrt (d);
        A[j * m + j] = d;
        for (std::size_t i = j + 1; i < m; ++i)
        {
            double s = A[i * m + j];
            for (std::size_t k = 0; k < j; ++k)
                s -= A[i * m + k] * A[j * m + k];
            A[i * m + j] = s / d;
        }
    }
    // Explicit inverse (m <= coarsestMaxFree), column j = A^-1 e_j by two
    // triangular solves; applied in single precision in the V-cycle.
    coarsestInverse_.assign (m * m, 0.0f);
    coarsestWork_.assign (m, 0.0f);
    std::vector<double> col (m);
    for (std::size_t j = 0; j < m; ++j)
    {
        for (std::size_t i = 0; i < m; ++i)
        {
            double sum = i == j ? 1.0 : 0.0;
            for (std::size_t k = 0; k < i; ++k)
                sum -= A[i * m + k] * col[k];
            col[i] = sum / A[i * m + i];
        }
        for (std::size_t ii = m; ii-- > 0;)
        {
            double sum = col[ii];
            for (std::size_t k = ii + 1; k < m; ++k)
                sum -= A[k * m + ii] * col[k];
            col[ii] = sum / A[ii * m + ii];
        }
        for (std::size_t i = 0; i < m; ++i)
            coarsestInverse_[j * m + i] = static_cast<float> (col[i]);
    }
    return true;
}

} // namespace taikor::candidates::vk
