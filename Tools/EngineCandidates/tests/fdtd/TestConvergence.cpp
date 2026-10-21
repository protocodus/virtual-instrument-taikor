// Test 8: convergence of the production take. The first 40 ms after the
// strike at the left capsule, at the production internal rate and at twice
// it (which also refines both head grids by sqrt 2, their spacing being set
// by the stability bound), compared in octave bands and in 1 kHz bands to
// 5 kHz (the output FIR's pass edge, Decimator.h); and at the production
// rate with finer air grids (L / 100 and L / 125 against L / 75, Courant
// 0.382 and 0.477). L / 150 (Courant 0.573, inside the 3-D limit 0.577) is
// out of reach of this PML: its explicit k^2 (z_x z_y + z_y z_z + z_z z_x) u
// term adds ~3 (sigma k)^2 = 0.7 to the leapfrog's (k Omega)^2 = 12 lambda^2
// = 3.93 in the layer's corners, past the bound 4, and the take blew up
// within 4 ms (with the reflecting box the same 6 ms stay bounded). Heavy: runs only with
// TAIKOR_FDTD_HEAVY (the air refinements with TAIKOR_FDTD_HEAVY_AIR).
// The check covers the octaves below 2 kHz, the band where test 1 puts the
// heads' modes within 1.4 % of the analytic values; the bands above are
// reported (head dispersion grows past 2 kHz, and the air grid's reaches
// several per cent by 4 kHz, FdtdDrumConfig.h).

#include "TestSupport.h"

#include "FdtdEngine.h"

#include <cmath>
#include <complex>
#include <cstdlib>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::fdtd::tests;

struct Take
{
    std::vector<double> trace;
    double rate = 0.0;
    ContactStatistics contact;
};

Take render (int oversample, int airCells, double seconds)
{
    setenv ("TAIKOR_FDTD_OVERSAMPLE", std::to_string (oversample).c_str(), 1);
    setenv ("TAIKOR_FDTD_AIR_CELLS", std::to_string (airCells).c_str(), 1);
    FdtdOptions options;
    options.selfChecks = false;
    FdtdEngine<float> engine (makeBergerPotential(), options);
    std::string error;
    engine.prepare (nagado().drum, nagado().mics, 48000.0, error);
    unsetenv ("TAIKOR_FDTD_OVERSAMPLE");
    unsetenv ("TAIKOR_FDTD_AIR_CELLS");
    engine.strike (nagado().stroke);
    Take take;
    take.rate = engine.internalRate();
    const int steps = static_cast<int> (seconds * take.rate);
    for (int n = 0; n < steps; ++n)
    {
        engine.advance (1);
        take.trace.push_back (engine.micPressure (0));
    }
    take.contact = engine.contact();
    return take;
}

// Energy of the Hann-windowed take in the band [f0, f1), by a direct DFT on
// a 5 Hz grid (the 40 ms window resolves 25 Hz).
double bandEnergy (const Take& take, double f0, double f1)
{
    const std::size_t n = take.trace.size();
    double total = 0.0;
    for (double f = f0; f < f1; f += 5.0)
    {
        const std::complex<double> step = std::polar (1.0, -2.0 * 3.14159265358979323846 * f / take.rate);
        std::complex<double> phasor (1.0, 0.0), sum (0.0, 0.0);
        for (std::size_t i = 0; i < n; ++i)
        {
            const double w = 0.5 - 0.5 * std::cos (2.0 * 3.14159265358979323846 * (i + 0.5) / static_cast<double> (n));
            sum += w * take.trace[i] * phasor;
            phasor *= step;
        }
        total += std::norm (sum / take.rate);
    }
    return total;
}

void compare (Report& report, const char* label, const Take& a, const Take& b)
{
    double worstLow = 0.0;
    std::string table;
    for (double f0 = 62.5; f0 < 16000.0; f0 *= 2.0)
    {
        const double ea = bandEnergy (a, f0, 2.0 * f0), eb = bandEnergy (b, f0, 2.0 * f0);
        const double db = 10.0 * std::log10 (eb / ea);
        if (f0 < 2000.0)
            worstLow = std::isfinite (db) ? std::max (worstLow, std::abs (db)) : HUGE_VAL;
        table += format ("%.0f-%.0f Hz %+.2f dB; ", f0, 2.0 * f0, db);
    }
    report.note (format ("%s: %s", label, table.c_str()));
    std::string fine;
    for (double f0 = 1000.0; f0 < 5000.0; f0 += 1000.0)
        fine += format ("%.0f-%.0f kHz %+.2f dB; ", 1e-3 * f0, 1e-3 * (f0 + 1000.0),
                        10.0 * std::log10 (bandEnergy (b, f0, f0 + 1000.0) / bandEnergy (a, f0, f0 + 1000.0)));
    report.note (format ("%s, to the FIR's pass edge: %s", label, fine.c_str()));
    report.note (format ("%s: contact %.4f / %.4f ms, peak %.2f / %.2f N, impulse %.5f / %.5f N s", label,
                         1e3 * a.contact.lastContact, 1e3 * b.contact.lastContact, a.contact.peakForce,
                         b.contact.peakForce, a.contact.impulse, b.contact.impulse));
    report.check (format ("%s: octave bands below 2 kHz", label), worstLow < 1.0,
                  format ("worst |level change| %.2f dB, limit 1 dB", worstLow));
}
} // namespace

FDTD_TEST (takeConvergence)
{
    if (! std::getenv ("TAIKOR_FDTD_HEAVY"))
        return;
    const auto numerics = FdtdNumerics::fromEnvironment();
    const double seconds = 0.04;
    const auto production = render (numerics.oversample, numerics.airCellsPerDepth, seconds);
    const auto finerRate = render (2 * numerics.oversample, numerics.airCellsPerDepth, seconds);
    compare (report, "rate x2 (heads refined by sqrt 2)", production, finerRate);
    if (std::getenv ("TAIKOR_FDTD_HEAVY_AIR"))
    {
        compare (report, "air L/100", production, render (numerics.oversample, 100, seconds));
        compare (report, "air L/125", production, render (numerics.oversample, 125, seconds));
    }
}
