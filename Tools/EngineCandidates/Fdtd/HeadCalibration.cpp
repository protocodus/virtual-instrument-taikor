#include "HeadCalibration.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <thread>

namespace taikor::candidates::fdtd
{

std::vector<double> modeShapeOnGrid (const HeadGrid& head, const StiffMode& mode, bool sine)
{
    const auto& g = head.grid();
    std::vector<double> shape (static_cast<std::size_t> (g.size()), 0.0);
    for (int j = 0; j < g.ny; ++j)
        for (int i = head.rowBegin (j); i < head.rowEnd (j); ++i)
        {
            const double x = g.x (i), y = g.y (j);
            const double r = std::min (std::hypot (x, y), head.radius());
            const double theta = std::atan2 (y, x);
            const double angular = sine ? std::sin (mode.m * theta) : std::cos (mode.m * theta);
            shape[static_cast<std::size_t> (g.index (i, j))] = clampedStiffShape (mode, r, head.radius()) * angular;
        }
    return shape;
}

double spectralPeakFrequency (const std::vector<double>& trace, double k, double lowHz, double highHz)
{
    return spectralPeakFrequency (std::vector<std::vector<double>> { trace }, k, lowHz, highHz);
}

double spectralPeakFrequency (const std::vector<std::vector<double>>& traces, double k, double lowHz, double highHz)
{
    const std::size_t n = traces.front().size();
    std::vector<std::vector<double>> windowed (traces.size(), std::vector<double> (n));
    for (std::size_t t = 0; t < traces.size(); ++t)
        for (std::size_t i = 0; i < n; ++i)
            windowed[t][i] = traces[t][i] * (0.5 - 0.5 * std::cos (2.0 * special::pi * (i + 0.5) / static_cast<double> (n)));
    const auto power = [&] (double f)
    {
        // Recursive phasor: e^{-i 2 pi f k} applied n times.
        const std::complex<double> step = std::polar (1.0, -2.0 * special::pi * f * k);
        double total = 0.0;
        for (const auto& w : windowed)
        {
            std::complex<double> phasor (1.0, 0.0), sum (0.0, 0.0);
            for (std::size_t i = 0; i < n; ++i)
            {
                sum += w[i] * phasor;
                phasor *= step;
                if ((i & 1023u) == 1023u)
                    phasor /= std::abs (phasor);
            }
            total += std::norm (sum);
        }
        return total;
    };
    const double duration = static_cast<double> (n) * k;
    const double coarse = 0.25 / duration;
    double best = lowHz, bestPower = -1.0;
    for (double f = lowHz; f <= highHz; f += coarse)
    {
        const double p = power (f);
        if (p > bestPower)
        {
            bestPower = p;
            best = f;
        }
    }
    // Golden-section search on the bracketing cell.
    double a = best - coarse, b = best + coarse;
    const double ratio = 0.5 * (std::sqrt (5.0) - 1.0);
    double c = b - ratio * (b - a), d = a + ratio * (b - a);
    double pc = power (c), pd = power (d);
    for (int it = 0; it < 60 && (b - a) > 1.0e-9 * best; ++it)
    {
        if (pc > pd)
        {
            b = d;
            d = c;
            pd = pc;
            c = b - ratio * (b - a);
            pc = power (c);
        }
        else
        {
            a = c;
            c = d;
            pc = pd;
            d = a + ratio * (b - a);
            pd = power (d);
        }
    }
    return 0.5 * (a + b);
}

double discreteFundamental (const HeadSpec& spec, double radius, double k, double spacing, double inclusionRadius)
{
    HeadGrid head;
    head.build (spec, radius, k, spacing, inclusionRadius, true);
    const auto modes = clampedStiffModes (spec, radius, 0, 400.0);
    const StiffMode& fundamental = modes.front();
    const auto shape = modeShapeOnGrid (head, fundamental, false);
    // Numerical: 24 periods give a Hann main lobe of +/- f/12, far narrower
    // than the gap to the next axisymmetric mode (f02 / f01 ~ 2.3).
    const int steps = static_cast<int> (std::ceil (24.0 / (fundamental.frequency * k)));
    const auto traces = simulateInVacuo (head, shape, { shape }, steps);
    return spectralPeakFrequency (traces.front(), k, 0.9 * fundamental.frequency, 1.02 * fundamental.frequency);
}

InclusionCalibration calibrateInclusionRadius (const HeadSpec& spec, double radius, double k, double spacing)
{
    InclusionCalibration result;
    result.analyticHz = clampedStiffModes (spec, radius, 0, 400.0).front().frequency;
    const double r1 = radius, r2 = radius - spacing;
    double f1 = 0.0, f2 = 0.0;
    {
        std::thread other ([&] { f2 = discreteFundamental (spec, radius, k, spacing, r2); });
        f1 = discreteFundamental (spec, radius, k, spacing, r1);
        other.join();
    }
    result.uncorrectedHz = f1;
    // Regula falsi on the bracket [a - h, a]: the discrete frequency is not
    // linear in the radius (the node set changes in small steps), so one
    // linear interpolation leaves up to ~0.06 %. Numerical: at most four
    // refinements, stopping within 1e-4 of the analytic value (test 1 allows
    // 5e-4); the best radius seen is kept.
    const double target = result.analyticHz;
    double ra = r1, fa = f1, rb = r2, fb = f2;
    result.inclusionRadius = r1;
    result.correctedHz = f1;
    for (int iteration = 0; iteration < 4 && std::abs (fb - fa) > 0.0; ++iteration)
    {
        double chosen = ra + (target - fa) * (rb - ra) / (fb - fa);
        chosen = std::clamp (chosen, radius - 1.5 * spacing, radius + 0.5 * spacing);
        const double fc = discreteFundamental (spec, radius, k, spacing, chosen);
        if (std::abs (fc - target) < std::abs (result.correctedHz - target))
        {
            result.inclusionRadius = chosen;
            result.correctedHz = fc;
        }
        if (std::abs (fc / target - 1.0) < 1.0e-4)
            break;
        if ((fc - target) * (fa - target) > 0.0)
        {
            ra = chosen;
            fa = fc;
        }
        else
        {
            rb = chosen;
            fb = fc;
        }
    }
    return result;
}

} // namespace taikor::candidates::fdtd
