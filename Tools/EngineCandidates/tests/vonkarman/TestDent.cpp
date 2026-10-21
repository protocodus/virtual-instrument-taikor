// (5b) The local dent against its continuum value. The dent of TestRatios is
// a discrete Poisson solve on each grid, so its change from grid to grid
// mixes the load's sampling with the potential's discretisation, and it has
// no continuum reference. Here the dent is the continuum one, sampled at the
// nodes: the static deflection of the tensioned membrane (w(a) = 0) under a
// uniform disk load of the test stroke's footprint radius b, centred on the
// head (-T Lap w = p inside r < b; amplitude-free, both energies are quartic):
//   w(r) = (1 - r^2 / b^2) / 2 + ln(a / b)  (r < b),   ln(a / r)  (r >= b),
//   w'(r) = -r / b^2                         (r < b),   -1 / r     (r >= b).
// Its continuum von Karman energy with an immovable edge comes from the
// radial quadrature of TestAnalytic.cpp (checked there against the closed
// forms); an independent P2 finite-element minimisation over u(r), a
// different method, gave V / V_B = 16.168 in verification. Continuum and
// grid use the same exact Berger reference V_B = C I^2 / (8 pi a^2).
//
// The discrete V_h / V_B must converge to it: monotonically, at the order
// the smooth interior sets (about 2; the rim staircase is first order, but
// the dent's energy sits under the stick, far from the rim), and the
// Richardson extrapolation of the two finest grids must land on it.

#include "TestUtil.h"

#include <algorithm>

namespace vktest
{

void testDent()
{
    const double a = headRadius();
    const double nu = batter().poissonRatio;
    const double b = shipping().stroke.footprintRadius;
    const auto slope = [a, b] (double r) { return r < b ? -r / (b * b) : -1.0 / r; };
    const auto continuum = axisymmetricEnergy (slope, a, nu, c11());
    const double exact = continuum.first / continuum.second;
    note (format ("continuum centred dent (b = %.2f mm): V/V_B = %.5f (radial quadrature)", 1e3 * b, exact));
    check (std::abs (exact / 16.168 - 1.0) < 1e-3,
           format ("radial quadrature %.5f against the independent finite-element value 16.168 (< 0.1 %%)", exact));

    // Numerical: spacings from the production batter's (4.9 mm) down by
    // factors of sqrt 2, to a/h ~ 227; solve to 1e-11 (exact to plotting).
    const std::vector<double> spacings { 0.0069, 0.0049, 0.00345, 0.00244, 0.00172 };
    std::vector<double> ratios;
    for (const double h : spacings)
    {
        Fixture f = makeFixture (a, h, 1.0e-11);
        if (! f.vk)
            return check (false, "prepare: " + f.error);
        f.vk->solverSettings().maxIterations = 2000;
        auto w = sampleField (f.grid, [&] (double x, double y)
        {
            const double r = std::hypot (x, y);
            return r < b ? 0.5 * (1.0 - r * r / (b * b)) + std::log (a / b) : std::log (a / r);
        });
        double peak = 0.0;
        for (const double v : w)
            peak = std::max (peak, std::abs (v));
        for (double& v : w)
            v *= 0.001 / peak; // 1 mm peak
        const double ratio = f.vk->evaluate (w.data(), nullptr) / bergerEnergy (f.grid, w.data(), a);
        ratios.push_back (ratio);
        note (format ("h = %.2f mm (a/h = %5.1f): V_h/V_B = %.5f, error %+.3f %%", 1e3 * h, a / h, ratio,
                      100.0 * (ratio / exact - 1.0)));
    }
    bool monotone = true;
    for (std::size_t i = 1; i < ratios.size(); ++i)
        monotone = monotone && std::abs (ratios[i] - exact) < std::abs (ratios[i - 1] - exact);
    const std::size_t n = ratios.size();
    // Observed order from the three finest (spacing ratio sqrt 2), then
    // Richardson on the two finest.
    const double d1 = ratios[n - 2] - ratios[n - 3], d2 = ratios[n - 1] - ratios[n - 2];
    const double order = std::log (d1 / d2) / std::log (spacings[n - 2] / spacings[n - 1]);
    const double factor = std::pow (spacings[n - 2] / spacings[n - 1], order);
    const double extrapolated = ratios[n - 1] + d2 / (factor - 1.0);
    check (monotone && order > 1.5 && order < 2.5,
           format ("monotone convergence at observed order %.2f (1.5-2.5) toward the continuum", order));
    check (std::abs (ratios[n - 1] / exact - 1.0) < 5e-3 && std::abs (extrapolated / exact - 1.0) < 2e-3,
           format ("finest grid %.5f (%+.3f %%, limit 0.5 %%), Richardson %.5f (%+.3f %%, limit 0.2 %%) against %.5f",
                   ratios[n - 1], 100.0 * (ratios[n - 1] / exact - 1.0), extrapolated,
                   100.0 * (extrapolated / exact - 1.0), exact));
}

} // namespace vktest
