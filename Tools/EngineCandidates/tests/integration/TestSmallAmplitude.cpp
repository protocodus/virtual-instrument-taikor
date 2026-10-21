// Test 2: at small amplitude the von Karman and Berger letters converge to
// the same linear drum. Both nonlinearities are quartic potentials of the
// batter's displacement, so their forces are cubic and, relative to the
// linear response, their effect on the take falls as the square of the
// impact speed; the two letters' difference must fall the same way and
// vanish with it.
//
// Setup: the production drum (lossy heads, PML, production margins) in
// double-precision air, so the differences are not buried under float air's
// rounding (-58 dB over 40 ms, tests/fdtd floatNullTest); both letters on the
// same grids (the Berger-bound spacing for the full stroke's energy, which
// covers the von Karman stress of these quiet strokes: its peak SAV tension
// is 867 N/m = 0.14 T at a quarter of the Don's speed, 57 N/m at 1/16, against
// the 2558 N/m = 0.42 T the grid admits); the rear head Berger in every run. Each speed also runs with the
// batter linear, the reference both letters converge to. The observable is
// the left capsule's pressure over the first 20 ms (sizingWindow: the contact,
// which lengthens at low speed as v^-1/5 under Hertz, and the attack).

#include "IntegrationSupport.h"

#include "VonKarmanPotential.h"

#include <cmath>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::integration;

double relative (const std::vector<double>& reference, const std::vector<double>& other)
{
    return std::sqrt (differenceEnergy (reference, other, 0, reference.size()) / energy (reference, 0, reference.size()));
}
} // namespace

INTEGRATION_TEST (smallAmplitudeConvergence)
{
    const auto& spec = nagado();
    FdtdOptions options;
    options.selfChecks = false;
    options.sizeBatterByStress = false;
    options.strokeEnergy = fdtd::strokeEnergy (spec.stroke);
    const int steps = static_cast<int> (std::ceil (sizingWindow * 96000.0));
    const double fractions[3] = { 0.25, 1.0 / 16.0, 1.0 / 64.0 };
    double eLetters[3], eVk[3], eBerger[3], peak[3];
    std::string error;
    for (int s = 0; s < 3; ++s)
    {
        StrokeSpec stroke = spec.stroke;
        stroke.impactSpeed *= fractions[s];
        auto berger = runTake<double> (makeBergerPotential(), options, stroke, steps, error);
        if (! error.empty())
        {
            report.check ("prepare", false, error);
            return;
        }
        // Later runs reuse the calibrated grids (identical by construction).
        options.batterInclusion = berger.engine->batterState().calibration.inclusionRadius;
        options.rearInclusion = berger.engine->rearState().calibration.inclusionRadius;
        auto vk = runTake<double> (std::make_unique<VonKarmanPotential>(), options, stroke, steps, error);
        FdtdOptions linearOptions = options;
        linearOptions.withBatterNonlinearity = false;
        auto linear = runTake<double> (makeBergerPotential(), linearOptions, stroke, steps, error);
        eLetters[s] = relative (berger.left, vk.left);
        eBerger[s] = relative (linear.left, berger.left);
        eVk[s] = relative (linear.left, vk.left);
        peak[s] = vk.engine->batterState().peakSavTension;
        report.note (format ("speed %.4f m/s (1/%.0f of the Don): |vK - Berger| / |Berger| = %.3e; |Berger - linear| / "
                             "|linear| = %.3e; |vK - linear| / |linear| = %.3e; vK peak SAV tension %.1f N/m (Berger %.2f); "
                             "batter h %.4f / %.4f mm",
                             stroke.impactSpeed, 1.0 / fractions[s], eLetters[s], eBerger[s], eVk[s], peak[s],
                             berger.engine->batterState().peakSavTension, 1e3 * berger.engine->batterState().grid.spacing(),
                             1e3 * vk.engine->batterState().grid.spacing()));
    }
    // Log-log slopes against the speed between successive speeds (ratio 4).
    const auto slope = [&] (const double* e, int i) { return std::log (e[i] / e[i + 1]) / std::log (4.0); };
    double worst = 0.0;
    for (int i = 0; i < 2; ++i)
        for (const double* e : { eLetters, eBerger, eVk })
            worst = std::max (worst, std::abs (slope (e, i) - 2.0));
    report.check ("the letters' difference and each letter's nonlinear effect fall as the speed squared", worst < 0.3,
                  format ("log-log slopes vs speed (1/4 -> 1/16, 1/16 -> 1/64): vK-Berger %.3f, %.3f; Berger-linear "
                          "%.3f, %.3f; vK-linear %.3f, %.3f (cubic forces: 2); worst deviation %.3f, limit 0.3",
                          slope (eLetters, 0), slope (eLetters, 1), slope (eBerger, 0), slope (eBerger, 1),
                          slope (eVk, 0), slope (eVk, 1), worst));
    report.check ("at 1/64 of the Don's speed the letters agree with each other and with the linear drum",
                  eLetters[2] < 1e-3 && eVk[2] < 1e-3 && eBerger[2] < 1e-3,
                  format ("relative L2 differences of the left capsule over 20 ms: vK-Berger %.2e, vK-linear %.2e, "
                          "Berger-linear %.2e, limit 1e-3 each",
                          eLetters[2], eVk[2], eBerger[2]));
    // Physics, not a check: von Karman's effect on the capsule is several
    // times Berger's at every speed, because its stress stays under the
    // stick, where the slopes are (peak SAV tension ~40 times Berger's
    // uniform one), while Berger spreads the same stretch over the head
    // (tests/vonkarman "ratios": a static dent under the footprint stores
    // 16-18 times Berger's energy). The units of both potentials are checked
    // against the continuum in potentialUnits.
    report.note (format ("(vK - linear) / (Berger - linear): %.2f, %.2f, %.2f at 1/4, 1/16, 1/64 of the Don's speed",
                         eVk[0] / eBerger[0], eVk[1] / eBerger[1], eVk[2] / eBerger[2]));
}
