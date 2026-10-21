// Test 1: the whole drum with the von Karman batter (both heads, air inside
// and out in a reflecting box, von Karman batter and Berger rear through SAV,
// the bachi collision) conserves its discrete energy when lossless, whatever
// the in-plane solve's accuracy, and with the hide losses and the barrel
// wall's boundary layer (WallLoss.h) the energy never increases and energy +
// dissipated work stays constant. The batter's grid is
// sized for the stroke as in the take (FdtdEngineSizing.cpp).

#include "IntegrationSupport.h"

#include "VonKarmanPotential.h"

#include <cmath>
#include <cstdlib>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::integration;

double lossPower (const FdtdEngine<double>& engine)
{
    const long long s = engine.stepsTaken();
    double total = 0.0;
    for (const auto* head : { &engine.batterState(), &engine.rearState() })
        total += head->grid.lossPower (head->w[s % 3].data(), head->w[(s + 1) % 3].data());
    return total;
}

struct Run
{
    double worstDrift = 0.0, worstRise = 0.0, e0 = 0.0, finalStick = 0.0, design = 0.0, spacing = 0.0;
    double peakSav = 0.0;
    long long steps = 0;
    int episodes = 0;
};

Run run (bool lossless, int steps, double tolerance)
{
    // Numerical settings of this test only (as tests/fdtd energyConservation):
    // a 2-cell air margin inside the reflecting box keeps the run short;
    // order 6 makes the coupling solve exact to round-off (nu^7 ~ 1e-22).
    // They stay set through strike(), whose sizing runs read them too.
    setenv ("TAIKOR_FDTD_AIR_MARGIN", "0.025", 1);
    setenv ("TAIKOR_FDTD_AIR_MARGIN_ABOVE", "0.025", 1);
    setenv ("TAIKOR_FDTD_NEUMANN_ORDER", "6", 1);
    FdtdOptions options;
    options.lossless = lossless;
    options.rigidBox = true;
    options.selfChecks = false;
    options.trackWallLoss = true;
    auto potential = std::make_unique<VonKarmanPotential>();
    potential->solverSettings().relativeTolerance = tolerance;
    FdtdEngine<double> engine (std::move (potential), options);
    std::string error;
    engine.prepare (nagado().drum, nagado().mics, 48000.0, error);
    engine.strike (nagado().stroke);
    unsetenv ("TAIKOR_FDTD_AIR_MARGIN");
    unsetenv ("TAIKOR_FDTD_AIR_MARGIN_ABOVE");
    unsetenv ("TAIKOR_FDTD_NEUMANN_ORDER");
    const auto total = [&engine]
    { return engine.energyHeads() + engine.energyAir() + engine.energyStick() + engine.energySav(); };
    engine.advance (1);
    Run result;
    result.e0 = total();
    double dissipated = 0.0, previous = result.e0;
    for (int n = 1; n < steps; ++n)
    {
        engine.advance (1);
        dissipated += lossPower (engine) / engine.internalRate();
        const double e = total();
        result.worstDrift = std::max (result.worstDrift, std::abs ((e + dissipated + engine.dissipatedWall()) / result.e0 - 1.0));
        result.worstRise = std::max (result.worstRise, (e - previous) / result.e0);
        previous = e;
    }
    result.steps = engine.stepsTaken();
    result.finalStick = engine.energyStick() / result.e0;
    result.episodes = engine.contact().episodes;
    result.design = engine.designTension (0);
    result.spacing = engine.batterState().grid.spacing();
    result.peakSav = engine.batterState().peakSavTension;
    return result;
}
} // namespace

INTEGRATION_TEST (energyConservation)
{
    // 8 ms at the internal rate: the whole first contact and its release.
    const int steps = 768;
    const auto lossless = run (true, steps, 3.0e-3);
    report.check ("lossless energy (heads + air + stick + contact + SAV) conserved, von Karman batter",
                  lossless.worstDrift < 1e-10,
                  format ("max |E/E0 - 1| = %.2e over %lld steps (E0 = %.4f J; stick keeps %.3f of it at the end; %d "
                          "contact episode(s); batter h %.4f mm sized for %.0f N/m, peak SAV tension %.0f N/m), "
                          "limit 1e-10",
                          lossless.worstDrift, lossless.steps, lossless.e0, lossless.finalStick, lossless.episodes,
                          1e3 * lossless.spacing, lossless.design, lossless.peakSav));
    // The SAV step conserves its energy for any gradient (Bilbao, Ducceschi &
    // Zama, J. Comput. Phys. 472, 2023): a solve 100 times looser than the
    // letter's leaves the balance exact.
    const auto loose = run (true, steps, 0.3);
    report.check ("lossless energy conserved with a loose in-plane solve", loose.worstDrift < 1e-10,
                  format ("in-plane relative tolerance 0.3: max |E/E0 - 1| = %.2e, limit 1e-10", loose.worstDrift));
    const auto lossy = run (false, steps, 3.0e-3);
    report.check ("with hide and wall losses the energy never increases", lossy.worstRise <= 1e-13,
                  format ("largest single-step rise %.2e of E0, limit 1e-13", lossy.worstRise));
    report.check ("with hide and wall losses E + dissipated work is constant", lossy.worstDrift < 1e-10,
                  format ("max |(E + D)/E0 - 1| = %.2e, limit 1e-10", lossy.worstDrift));
}
