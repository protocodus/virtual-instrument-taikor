// Test 3: the forte take does not depend on the in-plane solve's tolerance.
// The von Karman potential stops its in-plane solve at a relative residual
// of 3e-3 (InPlaneSolver::Settings, chosen in tests/vonkarman "tolerance" on
// a head in vacuo); here the whole drum plays the Don at that tolerance and
// at 1e-8, on the same grids (the second run takes the first one's sized
// batter), and the capsule pressures and the contact are compared. Double
// air, so the comparison measures the solve and not float air's rounding
// (which by itself reaches -58 dB over 40 ms, tests/fdtd floatNullTest).
// TAIKOR_INTEGRATION_TOLERANCE_MS lengthens the take (default 30 ms: the
// contact, its re-contact and the attack's largest stretching).

#include "IntegrationSupport.h"

#include "VonKarmanPotential.h"

#include <cmath>
#include <cstdlib>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::integration;
} // namespace

INTEGRATION_TEST (toleranceInsensitivity)
{
    const auto& spec = nagado();
    const double ms = std::getenv ("TAIKOR_INTEGRATION_TOLERANCE_MS") ? std::atof (std::getenv ("TAIKOR_INTEGRATION_TOLERANCE_MS")) : 30.0;
    const int steps = static_cast<int> (std::ceil (ms * 1e-3 * 96000.0));
    FdtdOptions options;
    options.selfChecks = false;
    std::string error;
    auto loose = std::make_unique<VonKarmanPotential>();
    auto* looseStats = loose.get();
    auto a = runTake<double> (std::move (loose), options, spec.stroke, steps, error);
    if (! error.empty())
    {
        report.check ("prepare", false, error);
        return;
    }
    const std::string looseReport = looseStats->diagnostics();
    options.batterDesignTension = a.engine->designTension (0);
    options.batterInclusion = a.engine->batterState().calibration.inclusionRadius;
    options.rearInclusion = a.engine->rearState().calibration.inclusionRadius;
    auto tight = std::make_unique<VonKarmanPotential>();
    tight->solverSettings().relativeTolerance = 1.0e-8;
    auto* tightStats = tight.get();
    auto b = runTake<double> (std::move (tight), options, spec.stroke, steps, error);
    const std::string tightReport = tightStats->diagnostics();
    const bool sameGrid = a.engine->batterState().grid.spacing() == b.engine->batterState().grid.spacing();
    // Windows: the contact (0-5 ms), the rest of the attack, the tail.
    const std::size_t e1 = static_cast<std::size_t> (0.005 * 96000.0), e2 = static_cast<std::size_t> (0.015 * 96000.0);
    const std::size_t end = static_cast<std::size_t> (steps);
    double worst = -300.0;
    std::string windows;
    for (const auto& [begin, stop] : { std::pair { std::size_t (0), e1 }, std::pair { e1, e2 }, std::pair { e2, end } })
    {
        const double left = differenceDb (b.left, a.left, begin, stop);
        const double right = differenceDb (b.right, a.right, begin, stop);
        worst = std::max ({ worst, left, right });
        windows += format ("%.1f-%.1f ms: %.1f / %.1f dB; ", 1e3 * begin / 96000.0, 1e3 * stop / 96000.0, left, right);
    }
    const auto& ca = a.engine->contact();
    const auto& cb = b.engine->contact();
    const double dPeak = ca.peakForce / cb.peakForce - 1.0, dImpulse = ca.impulse / cb.impulse - 1.0;
    report.check ("tightening the in-plane tolerance from 3e-3 to 1e-8 leaves the capsules unchanged",
                  sameGrid && worst < -60.0,
                  format ("difference energy re the 1e-8 take, left / right: %s worst %.1f dB, limit -60 dB (0.1 %% in "
                          "amplitude); same batter grid %s (h %.4f mm)",
                          windows.c_str(), worst, sameGrid ? "yes" : "NO", 1e3 * a.engine->batterState().grid.spacing()));
    report.check ("... and the contact", std::abs (dPeak) < 1e-3 && std::abs (dImpulse) < 1e-3,
                  format ("peak force %.3f vs %.3f N (%+.2e), impulse %.6f vs %.6f N s (%+.2e), limit 1e-3", ca.peakForce,
                          cb.peakForce, dPeak, ca.impulse, cb.impulse, dImpulse));
    report.note ("3e-3: " + looseReport);
    report.note ("1e-8: " + tightReport);
}
