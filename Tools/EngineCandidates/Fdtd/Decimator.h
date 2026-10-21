#pragma once

// The microphone signal leaves the engine through one linear-phase FIR
// (Kaiser-windowed sinc; J. F. Kaiser, "Nonrecursive digital filter design
// using the I0-sinh window function", Proc. IEEE ISCAS 1974; A. V. Oppenheim
// & R. W. Schafer, Discrete-Time Signal Processing, 3rd ed., sec. 7.6), then
// every factor-th sample is kept for the 48 kHz output.
//
// The filter's stop edge is the air grid's on-axis cut-off (FdtdEngine
// passes asin(lambda) / (pi k), 8.87 kHz for the production grid): above it
// the leapfrog air carries nothing along the axes, and just below it its
// group velocity falls to zero, so energy the stroke puts there stays
// trapped in the sealed, lossless barrel and rings for many seconds (the
// dev take held an 8-9 kHz ring at 0.6 dB/s). What feeds it is the contact
// force's own oscillation near 9 kHz: the rigid uniform-disk footprint of
// shared decision 3 has its first coupling null at beta b = 3.83, so the
// hide patch under the tip (~0.4 g) rides the Hertz spring almost decoupled
// from the rest of the head (a continuum membrane/disk/Hertz reference puts
// it at 9.55 kHz; the production heads at 8.83 kHz, 9.45 kHz at twice the
// rate; the engine's diagnostics report it). Limiting the output to the
// band the grid resolves is the usual practice for FDTD auralisation, whose
// dispersion bounds a scheme's valid bandwidth (K. Kowalczyk & M. van
// Walstijn, "Room acoustics simulation using 3-D compact explicit FDTD
// schemes", IEEE TASLP 19(1), 2011, 34-46). Numerical parameters: group
// delay 0.5 ms (48 samples at 96 kHz; FdtdEngine::process runs the drum that
// far ahead of its output, so the delay is cancelled), Kaiser attenuation
// 65 dB at the stop edge (beta = 0.1102 (A - 8.7)), so the transition band
// is (A - 8) / (2.285 * 2 pi * group delay) = 4.0 kHz wide and the pass edge
// sits near 4.9 kHz, where the axis phase error is ~5 %. The DC gain is one.

#include <algorithm>
#include <cmath>
#include <vector>

namespace taikor::candidates::fdtd
{

class Decimator
{
public:
    void design (int decimation, double internalRate, double stopHz)
    {
        factor = std::max (1, decimation);
        phase = 0;
        constexpr double pi = 3.14159265358979323846;
        constexpr double groupDelaySeconds = 0.5e-3, attenuationDb = 65.0;
        const int half = static_cast<int> (std::lround (groupDelaySeconds * internalRate));
        const int length = 2 * half + 1;
        const double beta = 0.1102 * (attenuationDb - 8.7);
        const double transitionHz = (attenuationDb - 8.0) / (2.285 * 2.0 * pi * (length - 1) / internalRate);
        stopEdge = stopHz;
        passEdge = stopHz - transitionHz;
        const double cutoff = stopHz - 0.5 * transitionHz;
        const auto besselI0 = [] (double x)
        {
            double sum = 1.0, term = 1.0;
            for (int k = 1; k < 80; ++k)
            {
                term *= (x / (2.0 * k)) * (x / (2.0 * k));
                sum += term;
            }
            return sum;
        };
        taps.clear();
        double total = 0.0;
        const double wc = 2.0 * pi * cutoff / internalRate;
        for (int n = 0; n < length; ++n)
        {
            const double t = n - half;
            const double sinc = t == 0.0 ? wc / pi : std::sin (wc * t) / (pi * t);
            const double ratio = t / half;
            const double window = besselI0 (beta * std::sqrt (std::max (0.0, 1.0 - ratio * ratio))) / besselI0 (beta);
            taps.push_back (sinc * window);
            total += taps.back();
        }
        for (auto& tap : taps)
            tap /= total;
        history.assign (taps.size(), 0.0);
        position = 0;
    }

    // Pushes one internal sample; returns true when an output sample is due.
    bool push (double x, double& y) noexcept
    {
        history[position] = x;
        position = position + 1 == history.size() ? 0 : position + 1;
        if (++phase < factor)
            return false;
        phase = 0;
        double sum = 0.0;
        std::size_t index = position;
        for (std::size_t n = taps.size(); n-- > 0;)
        {
            sum += taps[n] * history[index];
            index = index + 1 == history.size() ? 0 : index + 1;
        }
        y = sum;
        return true;
    }

    int delayInternalSamples() const noexcept { return static_cast<int> (taps.size() / 2); }
    double passEdgeHz() const noexcept { return passEdge; }
    double stopEdgeHz() const noexcept { return stopEdge; }

private:
    int factor = 1, phase = 0;
    std::vector<double> taps, history;
    std::size_t position = 0;
    double passEdge = 0.0, stopEdge = 0.0;
};

} // namespace taikor::candidates::fdtd
