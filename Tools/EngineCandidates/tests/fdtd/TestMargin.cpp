// Test 5 (second part): the absorbing layer does not colour the drum. The
// production take is rendered with the production air margins and with both
// margins doubled; the lowest partials' frequencies must agree within 1 % and
// their decay times within 2 %. Numerical estimator: each partial's decay is
// the least-squares slope of log |X(f)| over 0.4 s Hann windows sliding by
// 0.05 s through 0.05-1.0 s, so a neighbour 10 Hz away sits 4 bins out, past
// the window's main lobe. With 60 ms windows over a 0.4 s take (the earlier
// form of this test) the Don's 174.1 Hz partial (T60 0.85 s) beat against
// the 163.7 Hz one (T60 3.4 s, of similar amplitude over the first second)
// inside the window's +-33 Hz main lobe, and the "decay" measured was the beat
// pattern's: 2.66 against 2.90 s (-8.4 %) where a joint fit of two damped
// sinusoids gives 0.849 against 0.852 s (+0.35 %). The remaining changes are
// of both signs and about 1 % (estimation-limited: T60 up to ~3.7 s, so the
// 0.95 s span holds ~8-25 dB of decay). Heavy: runs only with
// TAIKOR_FDTD_HEAVY.

#include "TestSupport.h"

#include "FdtdEngine.h"
#include "HeadCalibration.h"

#include <cmath>
#include <complex>
#include <cstdlib>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::fdtd::tests;

std::vector<double> take (double margin, double marginAbove, double seconds, double& rate)
{
    setenv ("TAIKOR_FDTD_AIR_MARGIN", format ("%.6f", margin).c_str(), 1);
    setenv ("TAIKOR_FDTD_AIR_MARGIN_ABOVE", format ("%.6f", marginAbove).c_str(), 1);
    FdtdOptions options;
    options.selfChecks = false;
    FdtdEngine<float> engine (makeBergerPotential(), options);
    std::string error;
    engine.prepare (nagado().drum, nagado().mics, 48000.0, error);
    unsetenv ("TAIKOR_FDTD_AIR_MARGIN");
    unsetenv ("TAIKOR_FDTD_AIR_MARGIN_ABOVE");
    engine.strike (nagado().stroke);
    rate = engine.internalRate();
    std::vector<double> trace;
    const int steps = static_cast<int> (seconds * rate);
    for (int n = 0; n < steps; ++n)
    {
        engine.advance (1);
        trace.push_back (engine.micPressure (0));
    }
    return trace;
}

// |DTFT| of trace[begin, end) at f, Hann-windowed.
double magnitude (const std::vector<double>& trace, double rate, double f, std::size_t begin, std::size_t end)
{
    const std::complex<double> step = std::polar (1.0, -2.0 * 3.14159265358979323846 * f / rate);
    std::complex<double> phasor (1.0, 0.0), sum (0.0, 0.0);
    const double length = static_cast<double> (end - begin);
    for (std::size_t i = begin; i < end; ++i)
    {
        const double w = 0.5 - 0.5 * std::cos (2.0 * 3.14159265358979323846 * (i - begin + 0.5) / length);
        sum += w * trace[i] * phasor;
        phasor *= step;
    }
    return std::abs (sum);
}

// Decay rate from a least-squares line through log |X| over sliding windows
// (numerical: 0.4 s Hann windows, 0.05 s apart; see the head comment).
double decayRate (const std::vector<double>& trace, double rate, double f, double start, double stop)
{
    constexpr double span = 0.4, hop = 0.05;
    const std::size_t window = static_cast<std::size_t> (span * rate);
    std::vector<double> t, y;
    for (double t0 = start; t0 + span <= stop + 1e-9; t0 += hop)
    {
        const auto begin = static_cast<std::size_t> (t0 * rate);
        if (begin + window > trace.size())
            break;
        t.push_back (t0 + 0.5 * span);
        y.push_back (std::log (magnitude (trace, rate, f, begin, begin + window)));
    }
    double mt = 0.0, my = 0.0;
    for (std::size_t i = 0; i < t.size(); ++i)
    {
        mt += t[i];
        my += y[i];
    }
    mt /= t.size();
    my /= y.size();
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < t.size(); ++i)
    {
        num += (t[i] - mt) * (y[i] - my);
        den += (t[i] - mt) * (t[i] - mt);
    }
    return -num / den;
}
} // namespace

FDTD_TEST (marginDoubling)
{
    if (! std::getenv ("TAIKOR_FDTD_HEAVY"))
        return;
    const auto numerics = FdtdNumerics::fromEnvironment();
    const double seconds = 1.0;
    double rate = 0.0;
    const auto base = take (numerics.airMargin, numerics.airMarginAbove, seconds, rate);
    const auto wide = take (2.0 * numerics.airMargin, 2.0 * numerics.airMarginAbove, seconds, rate);
    // Partials of the wide take between 40 and 400 Hz within 30 dB of the strongest.
    const std::size_t from = static_cast<std::size_t> (0.02 * rate);
    std::vector<double> grid, level;
    double strongest = 0.0;
    for (double f = 40.0; f <= 400.0; f += 0.5)
    {
        grid.push_back (f);
        level.push_back (magnitude (wide, rate, f, from, wide.size()));
        strongest = std::max (strongest, level.back());
    }
    double worstF = 0.0, worstT = 0.0;
    int count = 0;
    std::string table;
    for (std::size_t i = 1; i + 1 < grid.size(); ++i)
    {
        if (! (level[i] > level[i - 1] && level[i] >= level[i + 1] && level[i] > strongest * 0.0316))
            continue;
        const double fw = spectralPeakFrequency (std::vector<double> (wide.begin() + static_cast<long> (from), wide.end()),
                                                 1.0 / rate, grid[i] - 1.0, grid[i] + 1.0);
        const double fb = spectralPeakFrequency (std::vector<double> (base.begin() + static_cast<long> (from), base.end()),
                                                 1.0 / rate, grid[i] - 1.5, grid[i] + 1.5);
        const double sw = decayRate (wide, rate, fw, 0.05, seconds), sb = decayRate (base, rate, fb, 0.05, seconds);
        const double df = fb / fw - 1.0, dt = sw / sb - 1.0; // T60 ratio = sw / sb
        worstF = std::max (worstF, std::abs (df));
        worstT = std::max (worstT, std::abs (dt));
        ++count;
        table += format ("%.2f Hz: f %+.3f %%, T60 %.3f s vs %.3f s (%+.2f %%); ", fw, 100.0 * df,
                         6.9078 / sb, 6.9078 / sw, 100.0 * dt);
    }
    report.note (table);
    report.check ("doubling the air margins: partial frequencies", count > 0 && worstF < 0.01,
                  format ("worst change %.3f %% over %d partials (40-400 Hz, within 30 dB), limit 1 %%", 100.0 * worstF, count));
    report.check ("doubling the air margins: partial decay times", count > 0 && worstT < 0.02,
                  format ("worst T60 change %.2f %%, limit 2 %% (estimation-limited, see above)", 100.0 * worstT));
}
