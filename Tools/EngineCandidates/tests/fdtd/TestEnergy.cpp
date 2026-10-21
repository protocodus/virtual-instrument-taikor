// Test 3: the whole coupled drum (both heads, air inside and out in a
// reflecting box, Berger on both heads through SAV, the bachi collision)
// conserves its discrete energy when lossless, and with the hide losses and
// the barrel wall's boundary layer on (WallLoss.h; its ring of head-plane
// cells enters the heads' coupling solve) the energy never increases and
// energy + dissipated work stays constant.

#include "TestSupport.h"

#include "FdtdEngine.h"

#include <cmath>
#include <cstdlib>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::fdtd::tests;

// Heads' loss power at the last step n (w^{n+1}, w^{n-1} available).
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
    double worstDrift = 0.0, worstRise = 0.0, e0 = 0.0, finalStick = 0.0, wallShare = 0.0;
    long long steps = 0;
    int episodes = 0;
    std::vector<double> batter; // w at the last step
};

Run run (bool lossless, int steps, const char* epsilon = nullptr)
{
    if (epsilon)
        setenv ("TAIKOR_FDTD_SAV_EPSILON", epsilon, 1);
    // Numerical settings of this test only: a 2-cell air margin inside the
    // reflecting box keeps the run short; order 6 makes the coupling solve
    // exact to round-off (nu^7 ~ 1e-22).
    setenv ("TAIKOR_FDTD_AIR_MARGIN", "0.025", 1);
    setenv ("TAIKOR_FDTD_AIR_MARGIN_ABOVE", "0.025", 1);
    setenv ("TAIKOR_FDTD_NEUMANN_ORDER", "6", 1);
    FdtdOptions options;
    options.lossless = lossless;
    options.rigidBox = true;
    options.selfChecks = false;
    options.trackWallLoss = true;
    FdtdEngine<double> engine (makeBergerPotential(), options);
    std::string error;
    engine.prepare (nagado().drum, nagado().mics, 48000.0, error);
    unsetenv ("TAIKOR_FDTD_AIR_MARGIN");
    unsetenv ("TAIKOR_FDTD_AIR_MARGIN_ABOVE");
    unsetenv ("TAIKOR_FDTD_NEUMANN_ORDER");
    unsetenv ("TAIKOR_FDTD_SAV_EPSILON");
    engine.strike (nagado().stroke);
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
    result.wallShare = engine.dissipatedWall() / std::max (1e-300, dissipated + engine.dissipatedWall());
    result.finalStick = engine.energyStick() / result.e0;
    result.episodes = engine.contact().episodes;
    result.batter = engine.batterState().w[engine.stepsTaken() % 3];
    return result;
}
} // namespace

FDTD_TEST (energyConservation)
{
    // 8 ms at the internal rate: the whole first contact and its release.
    const auto lossless = run (true, 768);
    report.check ("lossless energy (heads + air + stick + contact + SAV) conserved",
                  lossless.worstDrift < 1e-10,
                  format ("max |E/E0 - 1| = %.2e over %lld steps (E0 = %.4f J, the stroke's kinetic energy; stick "
                          "keeps %.3f of it at the end; %d contact episode(s)), limit 1e-10",
                          lossless.worstDrift, lossless.steps, lossless.e0, lossless.finalStick, lossless.episodes));
    const auto lossy = run (false, 768);
    report.check ("with hide losses the energy never increases", lossy.worstRise <= 1e-13,
                  format ("largest single-step rise %.2e of E0, limit 1e-13", lossy.worstRise));
    report.check ("with hide and wall losses E + dissipated work is constant", lossy.worstDrift < 1e-10,
                  format ("max |(E + D)/E0 - 1| = %.2e, limit 1e-10 (the wall's share of D %.2e)", lossy.worstDrift,
                          lossy.wallShare));
    // The SAV regulariser (FdtdNumerics::savEpsilon, 1e-18 J against a
    // 0.43 J stroke) is numerical: halving it must not change the result.
    const auto halved = run (false, 768, "5e-19");
    double difference = 0.0, norm = 0.0;
    for (std::size_t l = 0; l < lossy.batter.size(); ++l)
    {
        difference += (halved.batter[l] - lossy.batter[l]) * (halved.batter[l] - lossy.batter[l]);
        norm += lossy.batter[l] * lossy.batter[l];
    }
    const double relative = std::sqrt (difference / norm);
    report.check ("halving the SAV regulariser leaves the batter unchanged", relative < 1e-8,
                  format ("relative L2 change of w after 8 ms: %.2e, limit 1e-8", relative));
}
