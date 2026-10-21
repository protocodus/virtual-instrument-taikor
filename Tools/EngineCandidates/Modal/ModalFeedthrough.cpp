// The contact-force term of the capsule pressure, and the retardation across
// the head while the stick drives it.
//
// The capsule hears p(t) = rho int h(tau) y''(t - tau) dtau per mode, h the
// delay density of the mode's transfer: the Rayleigh nodes of ModalMics.cpp,
// weight w at delay R/c (so G(s) = int h e^{-s tau} dtau), and for the rear head
// V delta(tau - R'/c) / (4 pi R'). A mode y'' + 2 alpha y' + w0^2 y = f with
// analytic amplitude A (y = Re A, y' = Re{s A}, s = -alpha + j w) obeys
//   dA/dt = s A - j f / w,   y'' = Re{s^2 A} + f,
// so A(t - tau) = e^{-s tau} A(t) + (j/w) int_0^tau e^{s (u - tau)} f(t - u) du and,
// exactly,
//   p(t) = Re{rho s^2 G(s) A(t)}
//        + rho int_0^inf f(t - u) [h(u) + Re{(j s^2/w) int_u^inf h(tau) e^{-s (tau - u)} dtau}] du.
// The first term is the pole-residue output of ModalMics.cpp, exact for the
// free decay; the second is what it misses while the mode is forced: the
// force's own share of the acceleration (the "feedthrough", h(u)) and the
// retardation of the forced change of A across the head (the second part,
// which also cancels the first term's pressure before sound can arrive, so the
// sum is causal). With the contact force f_k = -c_k F for every mode, the
// correction is -int F(t - u) K(u) du with one kernel per capsule,
//   K(u) = rho sum_k c_k [h_k(u) + Re{(j s_k^2/w_k) int_u^inf h_k(tau) e^{-s_k (tau - u)} dtau}],
// nonzero only up to the longest delay (about 3.3 ms here).
//
// The force is known per contact substep d, so K is kept integrated over lag
// bins [l d, (l + 1) d). With E_l = sum of node weights in bin l and
// S_l = sum w e^{-s R/c} over the same nodes (exact phases), the bin integral is
//   rho c [E_l + Re{(j s/w) (E_l - e^{s l d} S_l + (e^{s (l+1) d} - e^{s l d}) sum_{b>l} S_b)}],
// exact for the quadrature, and sum_b S_b is the same G(s) the pole-residue
// output uses, so the pre-arrival cancellation holds to rounding
// (tests/modal output_causality, output_pressure_integral).
//
// Two other forces act on the modes and carry no such term of their own: the
// Berger force (a modulation of each mode's restoring force by dT/T, 6 % at
// the attack peak here and falling with the square of the amplitude; it leaves
// 1.8 % rms of the first 10 ms, mostly a phase error of the upper partials, and
// 3 % of the uncorrected pressure integral: tests/modal output_retarded_reference,
// output_pressure_integral), and the damping couplings inside the clusters of
// close modes (ModalRuntime.h; 7e-4 of the peak before sound can arrive,
// output_causality). Both would need a kernel per mode.

#include "ModalBuild.h"
#include "ModalInternal.h"
#include "ModalRuntime.h"

#include <chrono>
#include <cmath>
#include <mutex>

namespace taikor::candidates::modal
{
namespace
{
// A constant force F over one substep d changes a mode's analytic amplitude by
// beta_c F, beta_c = -(j/w) (e^{s d} - 1)/s; the two-point update of
// ModalRuntime.h changes it by beta_d F, beta_d = d* xi (1 + j (mu_r - 1)/mu_i)
// (d* = 1/(1 + a* + b*), xi = d^2/2), equal to second order in w d. Scaling the
// retardation part of the kernel by beta_d / beta_c makes it cancel the
// pole-residue sum of the discrete trajectory exactly before sound can arrive
// (the ratio differs from 1 by 3e-3 at 21.6 kHz for d = 1/768 kHz).
std::complex<double> discreteRatio (const BankMode& m, double d)
{
    const std::complex<double> s (-m.alpha, m.omega);
    const double R = std::exp (-m.alpha * d), O = std::cos (m.omega * d);
    const double den = 1.0 + 2.0 * R * O + R * R;
    const double a = (1.0 - 2.0 * R * O + R * R) / den, b = 2.0 * (1.0 - R * R) / den;
    const std::complex<double> zeta = std::exp (s * d);
    const std::complex<double> mu = (zeta - 1.0) / (zeta + 1.0);
    const std::complex<double> discrete = (0.5 * d * d / (1.0 + a + b)) * std::complex<double> (1.0, (mu.real() - 1.0) / mu.imag());
    const std::complex<double> continuous = std::complex<double> (0.0, -1.0 / m.omega) * (zeta - 1.0) / s;
    return discrete / continuous;
}
} // namespace

void buildFeedthrough (Bank& bank, const ModalModel& model, const DrumSpec& drum, const MicSpec& mics,
                       const StrokeSpec& stroke, const std::vector<int>& original, double bin, int threads)
{
    const auto t0 = std::chrono::steady_clock::now();
    const double a = model.radius, rho = drum.airDensity, c = drum.soundSpeed, len = drum.depth;
    const double thetaS = std::atan2 (stroke.strikeY, stroke.strikeX);
    double longest = 0.0;
    double thetaMic[2], rhoMic[2];
    for (int j = 0; j < 2; ++j)
    {
        const auto& p = mics.capsule[static_cast<std::size_t> (j)];
        thetaMic[j] = std::atan2 (p[1], p[0]);
        rhoMic[j] = std::hypot (p[0], p[1]);
        longest = std::max ({ longest, std::hypot (rhoMic[j] + a, p[2]), rearDistance (len, p) });
    }
    const int bins = static_cast<int> (std::ceil (longest / c / bin)) + 2;
    for (auto& kernel : bank.feedthrough)
        kernel.assign (static_cast<std::size_t> (bins), 0.0);
    bank.feedthroughBin = bin;

    std::mutex lock;
    parallelFor (static_cast<int> (bank.modes.size()), threads, [&] (int b)
    {
        const BankMode& bm = bank.modes[static_cast<std::size_t> (b)];
        const CoupledMode& cm = model.modes[static_cast<std::size_t> (original[static_cast<std::size_t> (b)])];
        if (bm.contact == 0.0)
            return;
        const std::complex<double> s (-bm.alpha, bm.omega);
        double b1Max = 0.0;
        const auto profile = batterProfile (cm, model.batter, a, b1Max);
        const double volume = rearVolume (cm, model.rear);
        std::vector<double> e (static_cast<std::size_t> (bins));
        std::vector<std::complex<double>> sum (static_cast<std::size_t> (bins));
        std::array<std::vector<double>, 2> local;
        for (int j = 0; j < 2; ++j)
        {
            std::fill (e.begin(), e.end(), 0.0);
            std::fill (sum.begin(), sum.end(), 0.0);
            const auto add = [&] (double weight, double dist)
            {
                const double tau = dist / c;
                const auto l = static_cast<std::size_t> (std::floor (tau / bin));
                e[l] += weight;
                sum[l] += weight * std::exp (-s * tau);
            };
            if (! cm.batterIndex.empty())
            {
                const double angular = std::cos (cm.m * (thetaMic[j] - thetaS));
                const auto& p = mics.capsule[static_cast<std::size_t> (j)];
                rayleighNodes (cm.m, a, profile, b1Max, std::abs (s) / c, rhoMic[j], p[2], 1,
                               [&] (double weight, double dist) { add (angular * weight, dist); });
            }
            if (volume != 0.0)
            {
                const double path = rearDistance (len, mics.capsule[static_cast<std::size_t> (j)]);
                add (volume / (4.0 * pi * path), path);
            }
            // Bin integrals, from the last bin down (tail = sum_{b > l} S_b).
            auto& out = local[static_cast<std::size_t> (j)];
            out.assign (static_cast<std::size_t> (bins), 0.0);
            const std::complex<double> factor = std::complex<double> (0.0, 1.0) * s / bm.omega * discreteRatio (bm, bin);
            // e^{s l d} by a downward recurrence, refreshed every 256 bins.
            const std::complex<double> stepPhase = std::exp (s * bin), stepBack = std::exp (-s * bin);
            std::complex<double> tail = 0.0;
            std::complex<double> phase = std::exp (s * (bin * (bins - 1)));
            for (int l = bins - 1; l >= 0; --l)
            {
                const auto ul = static_cast<std::size_t> (l);
                const std::complex<double> next = phase * stepPhase;  // e^{s (l + 1) d}
                const std::complex<double> retarded = e[ul] - phase * sum[ul] + (next - phase) * tail;
                out[ul] = rho * bm.contact * (e[ul] + (factor * retarded).real());
                tail += sum[ul];
                phase = (l - 1) % 256 == 0 ? std::exp (s * (bin * (l - 1))) : phase * stepBack;
            }
        }
        const std::lock_guard<std::mutex> guard (lock);
        for (int j = 0; j < 2; ++j)
            for (int l = 0; l < bins; ++l)
                bank.feedthrough[static_cast<std::size_t> (j)][static_cast<std::size_t> (l)] += local[static_cast<std::size_t> (j)][static_cast<std::size_t> (l)];
    });
    bank.feedthroughSeconds = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
}

} // namespace taikor::candidates::modal
