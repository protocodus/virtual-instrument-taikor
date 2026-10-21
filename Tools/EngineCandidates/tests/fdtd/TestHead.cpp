// Test 1: the head alone, lossless, in vacuo: discrete eigenfrequencies of
// the production grid (after the staircase correction) against the analytic
// clamped stiff membrane, and the spurious doublet split of the staircase.
// Test 6: the Berger potential's gradient against finite differences of V.

#include "TestSupport.h"

#include "Analytic.h"
#include "FdtdDrumConfig.h"
#include "HeadCalibration.h"

#include <cmath>
#include <random>
#include <thread>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::fdtd::tests;

struct ModeResult
{
    StiffMode mode;
    bool sine = false;
    double discreteHz = 0.0;
};

std::vector<ModeResult> measureModes (const HeadSpec& spec, double radius, double k, double spacing,
                                      double inclusion, double maxHz, int maxOrder, int maxN, double seconds)
{
    HeadGrid head;
    head.build (spec, radius, k, spacing, inclusion, true);
    std::vector<ModeResult> results;
    // Every mode with m <= maxOrder and n <= maxN, plus a sparse set reaching
    // 2 kHz (each probe costs a projection per step).
    const auto extra = [] (int m, int n)
    {
        return (m == 0 && n % 3 == 0) || (n == 1 && m % 4 == 0) || (m == 10 && n == 5) || (m == 5 && n == 10)
               || (m == 16 && n == 4) || (m == 3 && n == 13);
    };
    for (const auto& mode : clampedStiffModes (spec, radius, 40, maxHz))
    {
        if (! ((mode.m <= maxOrder && mode.n <= maxN) || extra (mode.m, mode.n)))
            continue;
        results.push_back ({ mode, false, 0.0 });
        if (mode.m > 0 && mode.m % 2 == 0)
            results.push_back ({ mode, true, 0.0 });
    }
    std::vector<std::vector<double>> shapes;
    std::vector<double> initial (static_cast<std::size_t> (head.size()), 0.0);
    for (const auto& r : results)
    {
        auto shape = modeShapeOnGrid (head, r.mode, r.sine);
        double norm = std::sqrt (head.dot (shape.data(), shape.data()));
        for (auto& v : shape)
            v /= norm;
        for (std::size_t l = 0; l < shape.size(); ++l)
            initial[l] += shape[l];
        shapes.push_back (std::move (shape));
    }
    const int steps = static_cast<int> (seconds / k);
    const auto traces = simulateInVacuo (head, initial, shapes, steps);
    for (std::size_t i = 0; i < results.size(); ++i)
    {
        const double f = results[i].mode.frequency;
        results[i].discreteHz = spectralPeakFrequency (traces[i], k, 0.94 * f, 1.03 * f);
    }
    return results;
}

void headCase (Report& report, const char* label, const HeadSpec& spec, double radius)
{
    const FdtdNumerics numerics = FdtdNumerics::fromEnvironment();
    const double k = 1.0 / (48000.0 * numerics.oversample);
    // The production spacing: the stability bound at the tension the Don's
    // energy allows the SAV step (FdtdDrumConfig.h headSpacing).
    const double spacing = headSpacing (spec, k, bergerTensionBound (spec, radius, strokeEnergy (nagado().stroke)));
    const auto calibration = calibrateInclusionRadius (spec, radius, k, spacing);
    report.note (format ("%s: fs %.0f Hz, h %.4f mm (a/h %.2f), inclusion radius a - %.3f h; "
                         "(0,1) analytic %.4f Hz, staircase %.4f Hz (%+.3f %%), corrected %.4f Hz",
                         label, 1.0 / k, 1e3 * spacing, radius / spacing, (radius - calibration.inclusionRadius) / spacing,
                         calibration.analyticHz, calibration.uncorrectedHz,
                         100.0 * (calibration.uncorrectedHz / calibration.analyticHz - 1.0), calibration.correctedHz));
    report.check (format ("%s (0,1) after staircase correction", label),
                  std::abs (calibration.correctedHz / calibration.analyticHz - 1.0) < 5.0e-4,
                  format ("%.4f vs %.4f Hz (%+.4f %%), limit 0.05 %%", calibration.correctedHz, calibration.analyticHz,
                          100.0 * (calibration.correctedHz / calibration.analyticHz - 1.0)));
    const auto modes = measureModes (spec, radius, k, spacing, calibration.inclusionRadius, 2000.0, 8, 4, 0.3);
    double worst1k = 0.0, worst2k = 0.0, worstSplit = 0.0, worstSplitAll = 0.0, highest = 0.0;
    std::string worstSplitMode;
    std::string table;
    for (std::size_t i = 0; i < modes.size(); ++i)
    {
        const auto& r = modes[i];
        const double error = r.discreteHz / r.mode.frequency - 1.0;
        if (r.mode.frequency <= 1000.0)
            worst1k = std::max (worst1k, std::abs (error));
        worst2k = std::max (worst2k, std::abs (error));
        highest = std::max (highest, r.mode.frequency);
        table += format ("(%d,%d%s) %.2f/%.2f Hz %+.3f%%; ", r.mode.m, r.mode.n, r.sine ? "s" : "",
                         r.discreteHz, r.mode.frequency, 100.0 * error);
        if (r.sine && i > 0 && modes[i - 1].mode.m == r.mode.m && modes[i - 1].mode.n == r.mode.n)
        {
            const double split = std::abs (r.discreteHz - modes[i - 1].discreteHz) / r.mode.frequency;
            if (r.mode.m <= 8)
                worstSplit = std::max (worstSplit, split);
            if (split > worstSplitAll)
            {
                worstSplitAll = split;
                worstSplitMode = format ("(%d,%d) at %.0f Hz", r.mode.m, r.mode.n, r.mode.frequency);
            }
        }
    }
    report.note (table);
    report.check (format ("%s modes below 1 kHz", label), worst1k < 0.01,
                  format ("worst |error| %.3f %% (%zu modes measured, m <= 8 and n <= 4 complete), limit 1 %%", 100.0 * worst1k, modes.size()));
    report.check (format ("%s modes below 2 kHz", label), worst2k < 0.025,
                  format ("worst |error| %.3f %% up to %.0f Hz (head dispersion at this rate), limit 2.5 %%", 100.0 * worst2k, highest));
    // The staircase disk has only the square's symmetry, so cos(m theta) and
    // sin(m theta) of even m are different discrete modes (Bilbao, JASA 131(1),
    // 2012, Fig. 6). Bound: below the shipping engine's own doublet-split
    // prior of 0.28 %, so the grid cannot fake more beating than A already has.
    // The gate covers the complete set m <= 8: the Don's footprint (radius b
    // at r_s = 0.15 a) samples a mode of order m >= 10 through J_m(alpha r)
    // with alpha r <= 2.6, i.e. at least 100 dB below the low modes, so the
    // rim-hugging (m, 1) modes of the extended set, whose split follows the
    // pseudo-random rim staircase, are reported only.
    report.check (format ("%s spurious doublet split (even m <= 8)", label), worstSplit < 0.0028,
                  format ("worst cos/sin split %.4f %% of f, limit 0.28 %%; over every measured mode %.4f %% at %s",
                          100.0 * worstSplit, 100.0 * worstSplitAll, worstSplitMode.c_str()));
}
} // namespace

FDTD_TEST (headEigenfrequencies)
{
    const auto& drum = nagado().drum;
    Report batter, rear;
    std::thread other ([&] { headCase (rear, "rear", drum.rear, drum.radius); });
    headCase (batter, "batter", drum.batter, drum.radius);
    other.join();
    for (const auto* part : { &batter, &rear })
        for (const auto& c : part->checks)
            report.checks.push_back (c);
}

FDTD_TEST (bergerGradient)
{
    const auto& drum = nagado().drum;
    const double k = 1.0 / 48000.0;
    HeadGrid head;
    head.build (drum.batter, drum.radius, k, headSpacing (drum.batter, k), drum.radius, true);
    auto potential = makeBergerPotential();
    std::string error;
    potential->prepare (head.grid(), drum.batter, drum.radius, error);
    std::mt19937 random (7);
    std::normal_distribution<double> normal;
    std::vector<double> w (static_cast<std::size_t> (head.size()), 0.0), d (w), g (w), scratch (w);
    // A smooth field of ~5 mm plus a rough perturbation direction.
    for (int l = 0; l < head.size(); ++l)
        if (head.isInterior (l))
        {
            const double x = head.grid().x (l % head.grid().nx), y = head.grid().y (l / head.grid().nx);
            const double r = std::hypot (x, y) / drum.radius;
            w[static_cast<std::size_t> (l)] = 5e-3 * (1.0 - r * r) * (1.0 + 0.3 * x / drum.radius) + 1e-4 * normal (random);
            d[static_cast<std::size_t> (l)] = normal (random);
        }
    const double v0 = potential->evaluate (w.data(), g.data());
    double directional = 0.0;
    for (std::size_t l = 0; l < w.size(); ++l)
        directional += g[l] * d[l];
    // V is quartic, so the central difference has an eps^2 error term;
    // Richardson extrapolation of eps and eps/10 removes it, leaving round-off.
    std::vector<double> unused (w.size());
    std::vector<double> estimates;
    for (const double eps : { 1e-6, 1e-7 })
    {
        for (std::size_t l = 0; l < w.size(); ++l)
            scratch[l] = w[l] + eps * d[l];
        const double vp = potential->evaluate (scratch.data(), unused.data());
        for (std::size_t l = 0; l < w.size(); ++l)
            scratch[l] = w[l] - eps * d[l];
        const double vm = potential->evaluate (scratch.data(), unused.data());
        estimates.push_back ((vp - vm) / (2.0 * eps));
        report.note (format ("eps %.0e: finite difference %.12e vs gradient %.12e", eps, estimates.back(), directional));
    }
    const double extrapolated = (100.0 * estimates[1] - estimates[0]) / 99.0;
    const double worst = std::abs (extrapolated / directional - 1.0);
    report.check ("Berger dV/dw is the exact derivative of V", worst < 1e-9,
                  format ("V = %.6e J; Richardson-extrapolated finite difference vs gradient: relative mismatch %.2e (limit 1e-9)", v0, worst));
}
