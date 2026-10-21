#pragma once

// The bachi-head collision as an energy-conserving discrete gradient of a
// one-sided power-law potential (S. Bilbao, A. Torin & V. Chatziioannou,
// "Numerical modeling of collisions in musical instruments", Acta Acust.
// united Ac. 101(1), 2015, 155-173, sec. III; V. Chatziioannou & M. van
// Walstijn, "Energy conserving schemes for the simulation of musical
// instrument contact dynamics", JSV 339, 2015):
//
//   Phi(eta) = K / (alpha + 1) [eta]_+^(alpha + 1),
//   f^n = (Phi(eta^{n+1}) - Phi(eta^{n-1})) / (eta^{n+1} - eta^{n-1}),
//
// eta = (footprint average of w) - (stick tip height), positive in contact,
// f > 0 pushing the stick out and the head in. With every other term linear
// in the unknowns, eta^{n+1} = etaFree - m f, so the step reduces to one
// scalar equation in r = eta^{n+1} - eta^{n-1}, a0 = eta^{n-1}:
//
//   Gamma(r) = r + m (Phi(r + a0) - Phi(a0)) / r - (etaFree - a0) = 0.
//
// Phi is convex, so the secant slope is non-decreasing in r, Gamma' >= 1 and
// the root is unique; it lies in [b - m F(b), b], b = etaFree - a0, where
// safeguarded Newton converges in a few iterations.

#include <algorithm>
#include <cmath>

namespace taikor::candidates::fdtd
{

struct HertzContact
{
    double stiffness = 0.0; // K, N / m^alpha
    double exponent = 1.5;  // alpha

    double potential (double eta) const noexcept
    {
        return eta > 0.0 ? stiffness / (exponent + 1.0) * std::pow (eta, exponent + 1.0) : 0.0;
    }
    double force (double eta) const noexcept { return eta > 0.0 ? stiffness * std::pow (eta, exponent) : 0.0; }
    double stiffnessAt (double eta) const noexcept
    {
        return eta > 0.0 ? stiffness * exponent * std::pow (eta, exponent - 1.0) : 0.0;
    }

    // Secant slope F(r) = (Phi(a0 + r) - Phi(a0)) / r (the discrete force),
    // its r-derivative, and their r -> 0 limits.
    void secant (double a0, double r, double& value, double& slope) const noexcept
    {
        const double scale = std::max ({ std::abs (a0), std::abs (r), 1.0e-30 });
        if (std::abs (r) < 1.0e-9 * scale)
        {
            value = force (a0) + 0.5 * stiffnessAt (a0) * r;
            slope = 0.5 * stiffnessAt (a0);
            return;
        }
        const double p1 = potential (a0 + r), p0 = potential (a0);
        value = (p1 - p0) / r;
        slope = (force (a0 + r) - value) / r;
    }

    // Solves Gamma(r) = 0; returns r and the force f = F(r).
    double solve (double a0, double etaFree, double m, double& f) const noexcept
    {
        const double b = etaFree - a0;
        double fb, unused;
        secant (a0, b, fb, unused);
        double lo = b - m * fb, hi = b;
        if (hi - lo <= 0.0)
        {
            f = fb;
            return b;
        }
        double r = hi;
        for (int iteration = 0; iteration < 100; ++iteration)
        {
            double value, slope;
            secant (a0, r, value, slope);
            const double gamma = r + m * value - b;
            if (gamma > 0.0)
                hi = r;
            else
                lo = r;
            double next = r - gamma / (1.0 + m * slope);
            if (! (next > lo && next < hi))
                next = 0.5 * (lo + hi);
            const double tolerance = 1.0e-15 * std::max ({ std::abs (a0), std::abs (b), 1.0e-18 });
            if (std::abs (next - r) <= tolerance || hi - lo <= tolerance)
            {
                r = next;
                break;
            }
            r = next;
        }
        double value, slope;
        secant (a0, r, value, slope);
        f = value;
        return r;
    }
};

} // namespace taikor::candidates::fdtd
