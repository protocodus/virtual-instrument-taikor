// (8) The Rayleigh microphone transfer: the exact on-axis baffled piston at a
// complex frequency, an independent brute-force quadrature of a head mode off
// axis (which also checks the transfer's reciprocity: the same symmetric kernel
// e^{-sR/c}/(2 pi R) gives the blocked-pressure force of a capsule-position
// monopole on the mode), and convergence of the node counts.

#include "ModalBuild.h"
#include "ModalFixtures.h"
#include "ModalTestUtil.h"

#include <algorithm>
#include <complex>

using namespace taikor::candidates;
using namespace taikor::candidates::modal;

MODAL_TEST (mics_piston_on_axis)
{
    const double a = 0.39, c = 343.0, z = 0.16;
    for (double f : { 120.0, 2000.0, 15000.0 })
    {
        const std::complex<double> s (-40.0, 2.0 * pi * f);
        const auto g = rayleighIntegral (0, a, [] (double) { return 1.0; }, 0.0, s, c, 0.0, z, 1);
        // int_0^a e^{-s R/c} r/R dr with R = sqrt(z^2 + r^2): (c/s)(e^{-s z/c} - e^{-s sqrt(z^2 + a^2)/c}).
        const auto exact = (c / s) * (std::exp (-s * z / c) - std::exp (-s * std::sqrt (z * z + a * a) / c));
        char label[96];
        std::snprintf (label, sizeof label, "piston on axis at %.0f Hz (complex s): |G - exact| / |exact|", f);
        modaltest::below (label, std::abs (g - exact) / std::abs (exact), 1e-10);
    }
}

MODAL_TEST (mics_head_mode_bruteforce)
{
    HeadMaterial h;
    h.tension = 6146.67;
    h.arealDensity = 0.94;
    h.bending = 0.266213;
    h.radius = 0.39;
    const auto modes = solveHeadModes (h, 2.0 * pi * 6000.0);
    const double c = 343.0, rhoM = 0.185, z = 0.16;
    for (const auto* mode : { &modes[3][1], &modes[12][4] })
    {
        const std::complex<double> s (-25.0, mode->omega);
        const auto profile = [&] (double r) { return evaluateShape (*mode, h.radius, r).value; };
        const auto g = rayleighIntegral (mode->m, h.radius, profile, mode->b1, s, c, rhoM, z, 1);
        // Brute force: 2-D product rule (Gauss-Legendre 600 in r, 1440 uniform in theta).
        std::vector<double> x, w;
        special::gaussLegendre (600, x, w);
        std::complex<double> ref = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            const double r = 0.5 * h.radius * (x[i] + 1.0), wr = 0.5 * h.radius * w[i];
            const double u = profile (r);
            for (int t = 0; t < 1440; ++t)
            {
                const double th = 2.0 * pi * t / 1440.0;
                const double dx = r * std::cos (th) - rhoM, dy = r * std::sin (th);
                const double dist = std::sqrt (dx * dx + dy * dy + z * z);
                ref += wr * r * (2.0 * pi / 1440.0) * u * std::cos (mode->m * th) * std::exp (-s * dist / c) / (2.0 * pi * dist);
            }
        }
        char label[96];
        std::snprintf (label, sizeof label, "(%d,%d) mode off axis vs brute force: relative error", mode->m, mode->n);
        modaltest::below (label, std::abs (g - ref) / std::abs (ref), 1e-8);
    }
}

MODAL_TEST (mics_convergence)
{
    HeadMaterial h;
    h.tension = 6146.67;
    h.arealDensity = 0.94;
    h.bending = 0.266213;
    h.radius = 0.39;
    const auto modes = solveHeadModes (h, 2.0 * pi * 21600.0);
    const double c = 343.0, rhoM = 0.185, z = 0.16;
    // High orders far below their cut-on (m = 120 at 10-20 kHz) have transfers
    // of ~1e-18, at the round-off floor of the oscillating integrand; the error
    // that matters is absolute, against the transfers that are heard.
    double worstRel = 0.0, worstAbs = 0.0, largest = 0.0;
    for (int m : { 0, 5, 40, 120 })
    {
        const auto& set = modes[static_cast<std::size_t> (m)];
        for (std::size_t i = 0; i < set.size(); i += std::max<std::size_t> (1, set.size() / 4))
        {
            const HeadMode& mode = set[i];
            const std::complex<double> s (-(0.95 + 0.0036 * mode.b1 * mode.b1), mode.omega);
            const auto profile = [&] (double r) { return evaluateShape (mode, h.radius, r).value; };
            const auto g1 = rayleighIntegral (m, h.radius, profile, mode.b1, s, c, rhoM, z, 1);
            const auto g2 = rayleighIntegral (m, h.radius, profile, mode.b1, s, c, rhoM, z, 2);
            // Pressure per unit modal displacement is rho s^2 G: compare on that scale.
            const double p1 = std::abs (g1 - g2) * std::norm (s), p2 = std::abs (g2) * std::norm (s);
            largest = std::max (largest, p2);
            worstAbs = std::max (worstAbs, p1);
            if (p2 > 1e-6 * largest)
                worstRel = std::max (worstRel, p1 / p2);
        }
    }
    modaltest::below ("node counts doubled: max relative change of transfers above -120 dB", worstRel, 1e-10);
    modaltest::below ("node counts doubled: max change of any transfer relative to the largest", worstAbs / largest, 1e-12);
}
