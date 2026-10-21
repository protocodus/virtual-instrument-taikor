// Test 9: float32 air against double air on the production drum and stroke
// (a null test of the production precision): the engine's output (the
// capsule pressure through its valid-band filter, 48 kHz) of the two over
// the first tens of milliseconds after the strike, and the raw internal-rate
// capsule pressure for reference. Heavy: runs only with TAIKOR_FDTD_HEAVY.
//
// What limits float32 here: at the air's Courant number (0.29) a 60 Hz wave
// changes Psi by only ~0.4 % per step, while each stored value is rounded to
// eps |Psi|, so every step adds a relative phase noise ~ eps / (omega k) ~
// 1.5e-5 to the drum's lowest partials; as a random walk that is ~ -61 dB
// after 40 ms (measured: -58.6 dB, 50-500 Hz) and ~ -40 dB after 4.5 s, a
// ~1 % perturbation of each partial rather than a noise floor. Storing the
// step increment (Psi^n - Psi^{n-1}) instead of the older level would cut it
// by ~250x; the interior's potential offset, the other float hazard, is
// removed by FdtdEngine's recentring. The limits below record that state.

#include "TestSupport.h"

#include "FdtdEngine.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::fdtd::tests;

// The left capsule: raw at the internal rate (advance) when raw, else the
// engine's own 48 kHz output (process, 256-sample blocks).
template <typename Real>
std::vector<double> micTrace (int steps, bool rigidBox, bool raw)
{
    FdtdOptions options;
    options.selfChecks = false;
    options.rigidBox = rigidBox;
    FdtdEngine<Real> engine (makeBergerPotential(), options);
    std::string error;
    engine.prepare (nagado().drum, nagado().mics, 48000.0, error);
    engine.strike (nagado().stroke);
    std::vector<double> trace;
    if (raw)
        for (int n = 0; n < steps; ++n)
        {
            engine.advance (1);
            trace.push_back (engine.micPressure (0));
        }
    else
    {
        const int samples = static_cast<int> (steps * 48000.0 / engine.internalRate());
        std::vector<float> left (256), right (256);
        for (int done = 0; done < samples; done += 256)
        {
            const int count = std::min (256, samples - done);
            engine.process (left.data(), right.data(), count);
            trace.insert (trace.end(), left.begin(), left.begin() + count);
        }
    }
    return trace;
}
} // namespace

FDTD_TEST (floatNullTest)
{
    if (! std::getenv ("TAIKOR_FDTD_HEAVY"))
        return;
    const int steps = static_cast<int> ((std::getenv ("TAIKOR_FDTD_NULL_MS") ? std::atof (std::getenv ("TAIKOR_FDTD_NULL_MS")) : 30.0)
                                        * 1e-3 * 96000.0);
    for (const int mode : { 0, 1, 2 })
    {
        const bool rigid = mode == 2, raw = mode == 1;
        if (rigid && ! std::getenv ("TAIKOR_FDTD_NULL_RIGID"))
            continue;
        const auto single = micTrace<float> (steps, rigid, raw);
        const auto full = micTrace<double> (steps, rigid, raw);
        double signal = 0.0, difference = 0.0;
        std::string growth;
        for (std::size_t n = 0; n < full.size(); ++n)
        {
            signal += full[n] * full[n];
            difference += (single[n] - full[n]) * (single[n] - full[n]);
            if ((n + 1) % (full.size() / 6) == 0)
                growth += format (" %.1f", 10.0 * std::log10 (difference / signal));
        }
        if (const char* dump = std::getenv ("TAIKOR_FDTD_DUMP"))
            if (FILE* file = std::fopen (format ("%s.%d", dump, mode).c_str(), "wb"))
            {
                std::fwrite (single.data(), sizeof (double), single.size(), file);
                std::fwrite (full.data(), sizeof (double), full.size(), file);
                std::fclose (file);
            }
        const double db = 10.0 * std::log10 (difference / signal);
        const char* what = mode == 0 ? "engine output (48 kHz, valid band)"
                                     : (mode == 1 ? "raw capsule pressure (96 kHz, full band)" : "raw, rigid box");
        const double limit = -55.0;
        report.check (format ("float32 air against double air, %.0f ms, %s", steps / 96.0, what), db < limit,
                      format ("difference %.1f dB re the signal (energy; cumulative by sixths:%s), limit %.0f dB",
                              db, growth.c_str(), limit));
    }
}
