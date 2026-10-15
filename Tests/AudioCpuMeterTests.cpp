#include "DSP/AudioCpuMeter.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <cstdint>

int main()
{
    taikor::AudioCpuMeter meter;
    const auto near = [] (float a, float b) { return std::abs (a - b) < 0.001f; };
    for (double rate : { 48000.0, 96000.0 })
    {
        meter.reset (rate);
        if (meter.snapshot().average != 0.0f || meter.snapshot().overruns != 0) return 1;
        const int samples = static_cast<int> (rate / 100.0);
        meter.record (0.0025, samples);
        if (! near (meter.snapshot().average, 25.0f) || ! near (meter.snapshot().peak, 25.0f)) return 1;
        meter.record (0.0125, samples);
        const auto busy = meter.snapshot();
        if (! near (busy.average, 29.0f) || ! near (busy.peak, 125.0f) || busy.overruns != 1) return 1;
        meter.record (0.01, 0);
        meter.record (0.01, -1);
        meter.record (-1.0, samples);
        meter.record (std::numeric_limits<double>::quiet_NaN(), samples);
        meter.record (std::numeric_limits<double>::infinity(), samples);
        if (meter.snapshot().average != busy.average || meter.snapshot().overruns != busy.overruns) return 1;
        for (int i = 0; i < 200; ++i) meter.record (0.0, samples);
        if (meter.snapshot().average >= 0.1f || meter.snapshot().peak >= busy.peak
            || meter.snapshot().overruns != 1) return 1;
    }
    meter.reset (std::numeric_limits<double>::quiet_NaN());
    meter.record (0.0025, 480);
    if (! near (meter.snapshot().average, 25.0f) || meter.snapshot().overruns != 0) return 1;

    std::uint32_t seed = 0x12345678u;
    const auto next = [&seed]
    {
        seed = seed * 1664525u + 1013904223u;
        return seed;
    };
    for (const double rate : { 0.0, 1.0e-300, 1.0e12,
                               std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN() })
    {
        meter.reset (rate);
        for (int iteration = 0; iteration < 4096; ++iteration)
        {
            const int samples = static_cast<int> (next() % 2048u) - 4;
            const double elapsed = iteration % 173 == 0
                ? std::numeric_limits<double>::quiet_NaN()
                : static_cast<double> (next() % 2000u) / 1000.0;
            meter.record (elapsed, samples);
            const auto state = meter.snapshot();
            if (! std::isfinite (state.average) || ! std::isfinite (state.peak)
                || state.average < 0.0f || state.peak < 0.0f)
                return 1;
        }
    }
    std::cout << "Audio CPU meter deadline, peak, reset and invalid-input checks passed\n";
}
