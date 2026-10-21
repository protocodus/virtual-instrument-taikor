// Test 8 (heavy): what the batter's spacing alone does to the take. The two
// finite-difference letters differ in the batter's potential and, because
// the von Karman batter is sized by the stretching its stroke measures
// (FdtdEngineSizing.cpp) rather than by Berger's energy bound, in the
// batter's spacing (on the Don 4.9187 against 4.8330 mm, 1.8 %; rear head,
// air grid and output filter are identical). Three takes of the Don, 0.1 s,
// float air, default threads:
//  - V: the von Karman letter (its sized grid);
//  - B: Berger on V's batter grid (FdtdOptions::batterDesignTension = V's
//    design, the same staircase calibration);
//  - L: the Berger letter (its own grid).
// L - B is the spacing's effect, V - B the potential's on one grid. Both are
// measured on the capsules low-passed at the output FIR's pass edge
// (4.9 kHz, 4th-order Butterworth forward and backward, at the internal
// rate), over 0-5, 5-20, 20-50 and 50-100 ms. Run by name:
// TaikorCandidateTests_integration gridParity.

#include "IntegrationSupport.h"

#include "VonKarmanPotential.h"

#include <cmath>
#include <numbers>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::integration;

// y = x through a 4th-order Butterworth low-pass (two RBJ-cookbook biquads,
// Q = 1 / (2 cos(pi/8)) and 1 / (2 cos(3 pi/8))), then again time-reversed.
std::vector<double> lowPass (const std::vector<double>& x, double cutoff, double rate)
{
    std::vector<double> y = x;
    const double w0 = 2.0 * std::numbers::pi * cutoff / rate;
    for (const double q : { 1.0 / (2.0 * std::cos (std::numbers::pi / 8.0)), 1.0 / (2.0 * std::cos (3.0 * std::numbers::pi / 8.0)) })
    {
        const double alpha = std::sin (w0) / (2.0 * q), c = std::cos (w0), a0 = 1.0 + alpha;
        const double b0 = (1.0 - c) / 2.0 / a0, b1 = (1.0 - c) / a0, b2 = b0, a1 = -2.0 * c / a0, a2 = (1.0 - alpha) / a0;
        for (int pass = 0; pass < 2; ++pass)
        {
            double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0;
            const auto n = y.size();
            for (std::size_t s = 0; s < n; ++s)
            {
                const std::size_t i = pass == 0 ? s : n - 1 - s;
                const double in = y[i], out = b0 * in + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
                x2 = x1;
                x1 = in;
                y2 = y1;
                y1 = out;
                y[i] = out;
            }
        }
    }
    return y;
}
} // namespace

INTEGRATION_HEAVY_TEST (gridParity)
{
    const auto& spec = nagado();
    const int steps = static_cast<int> (std::lround (0.1 * 96000.0));
    FdtdOptions options;
    options.selfChecks = false;
    std::string error;
    auto v = runTake<float> (std::make_unique<VonKarmanPotential>(), options, spec.stroke, steps, error);
    if (! error.empty())
        return report.check ("prepare", false, error);
    FdtdOptions onVonKarmanGrid = options;
    onVonKarmanGrid.batterDesignTension = v.engine->designTension (0);
    onVonKarmanGrid.batterInclusion = v.engine->batterState().calibration.inclusionRadius;
    onVonKarmanGrid.rearInclusion = v.engine->rearState().calibration.inclusionRadius;
    auto b = runTake<float> (makeBergerPotential(), onVonKarmanGrid, spec.stroke, steps, error);
    auto l = runTake<float> (makeBergerPotential(), options, spec.stroke, steps, error);
    if (! error.empty())
        return report.check ("prepare", false, error);
    const double rate = v.engine->internalRate(), cutoff = 4900.0;
    const bool sameOther = b.engine->batterState().grid.spacing() == v.engine->batterState().grid.spacing()
                           && l.engine->rearState().grid.spacing() == b.engine->rearState().grid.spacing();
    const auto filtered = [&] (const std::vector<double>& x) { return lowPass (x, cutoff, rate); };
    const std::vector<double> sides[3][2] = { { filtered (v.left), filtered (v.right) },
                                              { filtered (b.left), filtered (b.right) },
                                              { filtered (l.left), filtered (l.right) } };
    std::string text;
    double worstMargin = 300.0, worstGrid = -300.0;
    const double edges[] = { 0.0, 0.005, 0.020, 0.050, 0.100 };
    for (int w = 0; w < 4; ++w)
    {
        const auto begin = static_cast<std::size_t> (edges[w] * rate), end = static_cast<std::size_t> (edges[w + 1] * rate);
        double grid = -300.0, potential = -300.0;
        for (int side = 0; side < 2; ++side)
        {
            grid = std::max (grid, differenceDb (sides[1][side], sides[2][side], begin, end));
            potential = std::max (potential, differenceDb (sides[1][side], sides[0][side], begin, end));
        }
        worstGrid = std::max (worstGrid, grid);
        worstMargin = std::min (worstMargin, potential - grid);
        text += format ("%.0f-%.0f ms: spacing %.1f dB, potential %.1f dB; ", 1e3 * edges[w], 1e3 * edges[w + 1], grid, potential);
    }
    const auto& cb = b.engine->contact();
    const auto& cl = l.engine->contact();
    report.note (format ("batter h: von Karman letter %.4f mm, Berger on it %.4f mm, Berger letter %.4f mm; contact of "
                         "Berger on the two grids: peak %.2f / %.2f N, impulse %.5f / %.5f N s",
                         1e3 * v.engine->batterState().grid.spacing(), 1e3 * b.engine->batterState().grid.spacing(),
                         1e3 * l.engine->batterState().grid.spacing(), cb.peakForce, cl.peakForce, cb.impulse, cl.impulse));
    // Numerical limits: the spacing's effect must be negligible beside the
    // mechanism the letters compare, i.e. 20 dB (1 % of the difference
    // energy) below the potential's effect in every window, and below -30 dB
    // re the take.
    report.check ("the batter's spacing changes the take far less than the potential does",
                  sameOther && worstGrid < -30.0 && worstMargin > 20.0,
                  format ("difference energy re Berger on the von Karman grid, below %.1f kHz, worst of the two "
                          "capsules: %s worst spacing effect %.1f dB (limit -30 dB), smallest margin to the "
                          "potential's effect %.1f dB (limit 20 dB)",
                          1e-3 * cutoff, text.c_str(), worstGrid, worstMargin));
}
