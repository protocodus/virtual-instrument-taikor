// Test 4: units. The drum's nodal equation is in newtons (sigma h^2 d_tt w_l
// = ... - dV/dw_l) and its energy balance in joules, so both potentials must
// return V in J and dV/dw in N on the drum's own batter grid. Checked against
// the continuum for the clamped (0,1) membrane shape w = A J0(j01 r / a):
//  - Berger: V_B = C / (8 pi a^2) (int |grad w|^2 dA)^2, C = E h / (1 - nu^2)
//    (FdtdContracts.h; xi = C / (2 pi a^2) and V = xi / 4 G^2);
//  - von Karman with an immovable edge: V = min_u pi C int [eps_r^2 +
//    eps_t^2 + 2 nu eps_r eps_t] r dr, eps_r = u' + w'^2 / 2, eps_t = u / r,
//    u(0) = u(a) = 0, minimised here by linear finite elements in r (an
//    independent solve: 1-D, nothing shared with the 2-D in-plane solver);
//  - the gradient in newtons: the directional derivative of V along w, which
//    for a quartic potential is exactly 4 V (Euler), against grad V . w.
// On the staircase grid the energies carry its O(h/a) rim error (the grid's
// disk is the staircase-corrected inclusion radius, HeadCalibration.h), so the
// limits are 5 % on each energy and 1 % on their ratio, whose rim errors
// largely cancel (tests/vonkarman "ratios": 0.7825 against 0.7855 at 4.8 mm).
// A units error is a factor of h^2 ~ 2e-5, a/h ~ 80 or 2 pi.

#include "IntegrationSupport.h"

#include "SpecialFunctions.h"
#include "VonKarmanPotential.h"

#include <cmath>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::integration;

constexpr double pi = 3.14159265358979323846;

// Continuum stretching energies (J) of w = A J0(j01 r / a): {von Karman, Berger}.
std::pair<double, double> continuumEnergies (double A, double a, const HeadSpec& head)
{
    namespace sf = taikor::candidates::special;
    const double nu = head.poissonRatio;
    const double C = head.youngsModulus * head.thickness / (1.0 - nu * nu);
    const double j01 = sf::besselJZeros (0, 3.0)[0][0];
    const auto slope = [&] (double r) { return -A * j01 / a * sf::besselJ (1, j01 * r / a); };
    // Numerical: 20000 elements (midpoint rule, O(dr^2): ~1e-8 relative).
    const int M = 20000;
    const double dr = a / M;
    // Energy per element e: pi C dr r_m [(u' + q)^2 + (u_m / r_m)^2 + 2 nu (u' + q) u_m / r_m],
    // u' = (u_{e+1} - u_e) / dr, u_m = (u_e + u_{e+1}) / 2, q = w'(r_m)^2 / 2: a
    // quadratic in (u_e, u_{e+1}); its minimiser solves a tridiagonal system.
    std::vector<double> diag (static_cast<std::size_t> (M + 1), 0.0), off (static_cast<std::size_t> (M), 0.0),
        rhs (static_cast<std::size_t> (M + 1), 0.0);
    double constant = 0.0, gSum = 0.0;
    for (int e = 0; e < M; ++e)
    {
        const double rm = (e + 0.5) * dr, q = 0.5 * slope (rm) * slope (rm);
        const double wgt = pi * C * dr * rm;
        // strain = (a0 u_e + a1 u_{e+1}) + q for eps_r and (b0 u_e + b1 u_{e+1}) for eps_t
        const double a0 = -1.0 / dr, a1 = 1.0 / dr, b0 = 0.5 / rm, b1 = 0.5 / rm;
        const double c[2][2] = { { a0 * a0 + b0 * b0 + 2 * nu * a0 * b0, a0 * a1 + b0 * b1 + nu * (a0 * b1 + a1 * b0) },
                                 { 0.0, a1 * a1 + b1 * b1 + 2 * nu * a1 * b1 } };
        diag[static_cast<std::size_t> (e)] += wgt * c[0][0];
        diag[static_cast<std::size_t> (e + 1)] += wgt * c[1][1];
        off[static_cast<std::size_t> (e)] += wgt * c[0][1];
        rhs[static_cast<std::size_t> (e)] -= wgt * q * (a0 + nu * b0);
        rhs[static_cast<std::size_t> (e + 1)] -= wgt * q * (a1 + nu * b1);
        constant += wgt * q * q;
        gSum += 2.0 * pi * rm * dr * slope (rm) * slope (rm);
    }
    // u_0 = u_M = 0: solve the interior tridiagonal system (Thomas).
    std::vector<double> u (static_cast<std::size_t> (M + 1), 0.0), cp (static_cast<std::size_t> (M + 1), 0.0),
        dp (static_cast<std::size_t> (M + 1), 0.0);
    for (int i = 1; i < M; ++i)
    {
        const double lower = i > 1 ? off[static_cast<std::size_t> (i - 1)] : 0.0;
        const double m = diag[static_cast<std::size_t> (i)] - lower * cp[static_cast<std::size_t> (i - 1)];
        cp[static_cast<std::size_t> (i)] = off[static_cast<std::size_t> (i)] / m;
        dp[static_cast<std::size_t> (i)] = (rhs[static_cast<std::size_t> (i)] - lower * dp[static_cast<std::size_t> (i - 1)]) / m;
    }
    for (int i = M - 1; i >= 1; --i)
        u[static_cast<std::size_t> (i)] = dp[static_cast<std::size_t> (i)] - cp[static_cast<std::size_t> (i)] * u[static_cast<std::size_t> (i + 1)];
    // At the minimiser the energy is constant - (1/2) u^T rhs... evaluate it directly instead.
    double vk = 0.0;
    for (int e = 0; e < M; ++e)
    {
        const double rm = (e + 0.5) * dr, q = 0.5 * slope (rm) * slope (rm);
        const double er = (u[static_cast<std::size_t> (e + 1)] - u[static_cast<std::size_t> (e)]) / dr + q;
        const double et = 0.5 * (u[static_cast<std::size_t> (e)] + u[static_cast<std::size_t> (e + 1)]) / rm;
        vk += pi * C * dr * rm * (er * er + et * et + 2.0 * nu * er * et);
    }
    (void) constant;
    const double berger = C / (8.0 * pi * a * a) * gSum * gSum;
    return { vk, berger };
}
} // namespace

INTEGRATION_TEST (potentialUnits)
{
    namespace sf = taikor::candidates::special;
    const auto& spec = nagado();
    // The letter's batter grid, as the take builds it (sized for the Don).
    FdtdOptions options;
    options.selfChecks = false;
    FdtdEngine<float> engine (std::make_unique<VonKarmanPotential>(), options);
    std::string error;
    if (! engine.prepare (spec.drum, spec.mics, 48000.0, error))
    {
        report.check ("prepare", false, error);
        return;
    }
    engine.strike (spec.stroke);
    const auto& grid = engine.batterState().grid.grid();
    const double a = spec.drum.radius, A = 1.0e-3;
    const double j01 = sf::besselJZeros (0, 3.0)[0][0];
    std::vector<double> w (static_cast<std::size_t> (grid.size()), 0.0), gradient (w.size(), 0.0);
    for (int j = 0; j < grid.ny; ++j)
        for (int i = 0; i < grid.nx; ++i)
        {
            const auto l = static_cast<std::size_t> (grid.index (i, j));
            const double r = std::hypot (grid.x (i), grid.y (j));
            w[l] = grid.interior[l] != 0 && r < a ? A * sf::besselJ (0, j01 * r / a) : 0.0;
        }
    VonKarmanPotential vk;
    vk.prepare (grid, spec.drum.batter, a, error);
    vk.solverSettings().relativeTolerance = 1e-12;
    vk.solverSettings().maxIterations = 2000;
    auto berger = makeBergerPotential();
    berger->prepare (grid, spec.drum.batter, a, error);
    const double vVk = vk.evaluate (w.data(), gradient.data());
    double euler = 0.0;
    for (std::size_t l = 0; l < w.size(); ++l)
        euler += gradient[l] * w[l];
    const double vB = berger->evaluate (w.data(), gradient.data());
    double eulerB = 0.0;
    for (std::size_t l = 0; l < w.size(); ++l)
        eulerB += gradient[l] * w[l];
    const auto [cVk, cB] = continuumEnergies (A, a, spec.drum.batter);
    report.check ("Berger energy in joules on the batter grid", std::abs (vB / cB - 1.0) < 0.05,
                  format ("A = 1 mm (0,1): grid %.6e J, continuum %.6e J (%+.2f %%), limit 5 %%", vB, cB, 100.0 * (vB / cB - 1.0)));
    report.check ("von Karman energy in joules on the batter grid", std::abs (vVk / cVk - 1.0) < 0.05,
                  format ("grid %.6e J, continuum (1-D finite elements, immovable edge) %.6e J (%+.2f %%), limit 5 %%", vVk,
                          cVk, 100.0 * (vVk / cVk - 1.0)));
    const double ratio = vVk / vB, cRatio = cVk / cB;
    report.check ("von Karman over Berger on the batter grid", std::abs (ratio / cRatio - 1.0) < 0.01,
                  format ("%.5f against the continuum %.5f (%+.2f %%), limit 1 %%", ratio, cRatio, 100.0 * (ratio / cRatio - 1.0)));
    report.check ("gradients in newtons (Euler: grad V . w = 4 V for a quartic V)",
                  std::abs (euler / (4.0 * vVk) - 1.0) < 1e-6 && std::abs (eulerB / (4.0 * vB) - 1.0) < 1e-12,
                  format ("von Karman %.3e, Berger %.3e relative, limits 1e-6 (the in-plane solve to 1e-12) and 1e-12",
                          euler / (4.0 * vVk) - 1.0, eulerB / (4.0 * vB) - 1.0));
    report.note (format ("batter grid h %.4f mm, %d x %d nodes", 1e3 * grid.h, grid.nx, grid.ny));
}
