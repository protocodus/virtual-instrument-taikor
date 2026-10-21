// The in-plane solver's algebra: grid transfers are adjoint, every level
// operator is symmetric, the V-cycle is a symmetric positive definite
// preconditioner, and PCG converges at a grid-independent rate.

#include "InVacuoHead.h"

#include <algorithm>

namespace vktest
{

namespace
{
void randomLevel (const vk::InPlaneLevel& L, std::vector<double>& x, std::uint64_t& state)
{
    x.assign (static_cast<std::size_t> (L.size()), 0.0);
    for (int c = 0; c < 2; ++c)
        for (int k = 0; k < L.N; ++k)
        {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            const bool free = (c == 0 ? L.freeU : L.freeV)[static_cast<std::size_t> (k)] != 0;
            x[static_cast<std::size_t> (L.base (c) + k)] = free ? static_cast<double> (state >> 11) * 1.1102230246251565e-16 - 0.5 : 0.0;
        }
}

double dotv (const std::vector<double>& a, const std::vector<double>& b)
{
    double s = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k)
        s += a[k] * b[k];
    return s;
}
} // namespace

void testSolver()
{
    const double h = 1.001 * InVacuoHead::minimumSpacing (batter(), 48000.0);
    for (bool galerkin : { true, false })
    {
        const MembraneGrid grid = makeDiskGrid (headRadius(), h);
        VonKarmanPotential vk;
        vk.solverSettings().galerkinCoarse = galerkin;
        std::string error;
        if (! vk.prepare (grid, batter(), headRadius(), error))
            return check (false, "prepare: " + error);
        auto& solver = const_cast<vk::InPlaneSolver&> (vk.solver());
        std::uint64_t state = 99;
        double worstAdjoint = 0.0, worstSymmetry = 0.0;
        for (int l = 0; l < solver.levelCount(); ++l)
        {
            const auto& L = solver.level (l);
            std::vector<double> x, y, kx (static_cast<std::size_t> (L.size())), ky (kx.size());
            randomLevel (L, x, state);
            randomLevel (L, y, state);
            solver.apply (l, x.data(), kx.data());
            solver.apply (l, y.data(), ky.data());
            worstSymmetry = std::max (worstSymmetry, std::abs (dotv (y, kx) - dotv (x, ky)) / std::abs (dotv (x, kx)));
            if (l > 0)
            {
                const auto& F = solver.level (l - 1);
                std::vector<double> xf, px (static_cast<std::size_t> (F.size())), ry (static_cast<std::size_t> (L.size()));
                randomLevel (F, xf, state);
                solver.prolong (l, x.data(), px.data());
                solver.restrict_ (l, xf.data(), ry.data());
                worstAdjoint = std::max (worstAdjoint, std::abs (dotv (xf, px) - dotv (ry, x)) / std::abs (dotv (xf, px)));
            }
        }
        const char* kind = galerkin ? "Galerkin" : "rediscretised";
        double worstMargin = 1e9;
        std::string spectrum;
        for (int l = 0; l < solver.levelCount(); ++l)
        {
            const double r200 = solver.rayleighEstimate (l, 200), r800 = solver.rayleighEstimate (l, 800);
            worstMargin = std::min (worstMargin, solver.level (l).lambdaMax / r800);
            spectrum += format ("L%d %.4f/%.4f/%.4f ", l, r200, r800, solver.level (l).lambdaMax);
        }
        note (std::string (kind) + " Rayleigh 200/800 its and lambda_max used: " + spectrum);
        check (worstMargin >= 1.0, format ("%s: smoother bound lambda_max >= the 800-iteration Rayleigh estimate on every level (worst ratio %.3f)",
                                          kind, worstMargin));
        check (worstAdjoint < 1e-13, format ("%s: <P x, y> = <x, P^T y> on every level (worst %.1e)", kind, worstAdjoint));
        check (worstSymmetry < 1e-12, format ("%s: every level operator symmetric (worst %.1e)", kind, worstSymmetry));
        const auto& L0 = solver.level (0);
        std::vector<double> x, y, bx (static_cast<std::size_t> (L0.size())), by (bx.size());
        randomLevel (L0, x, state);
        randomLevel (L0, y, state);
        solver.precondition (x.data(), bx.data());
        solver.precondition (y.data(), by.data());
        const double asym = std::abs (dotv (y, bx) - dotv (x, by)) / std::sqrt (dotv (x, bx) * dotv (y, by));
        // Single-precision V-cycle: symmetric to float rounding.
        check (asym < 1e-5 && dotv (x, bx) > 0 && dotv (y, by) > 0,
               format ("%s: V-cycle symmetric to single-precision rounding (%.1e) and positive", kind, asym));
        const int shippingDegree = solver.settings().chebyshevDegree;
        for (int degree : { 1, 2, 3, 4 })
        {
            solver.settings().chebyshevDegree = degree;
            solver.settings().relativeTolerance = 1e-8;
            std::vector<double> b, sol (static_cast<std::size_t> (L0.size()), 0.0);
            randomLevel (L0, b, state);
            double residual = 0.0;
            const double t0 = cpuSeconds();
            const int its = solver.solve (b.data(), sol.data(), residual);
            const double dt = cpuSeconds() - t0;
            note (format ("%s, finest-level Chebyshev degree %d (coarse levels %d): cold PCG to 1e-8 in %d its (rate %.3f/it), %.2f ms CPU",
                          kind, degree, solver.settings().coarseChebyshevDegree, its, std::pow (residual, 1.0 / its), dt * 1e3));
            if (galerkin && degree == shippingDegree)
                check (std::pow (residual, 1.0 / its) < 0.25,
                       format ("Galerkin, shipping smoother: PCG contracts the residual by %.3f per iteration (< 0.25)", std::pow (residual, 1.0 / its)));
        }
        solver.settings().chebyshevDegree = shippingDegree;
    }
}

} // namespace vktest
