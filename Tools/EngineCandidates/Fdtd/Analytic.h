#pragma once

// Closed-form references the finite-difference drum is verified against.

#include "CandidateEngine.h"
#include "SpecialFunctions.h"

#include <cmath>
#include <vector>

namespace taikor::candidates::fdtd
{

// A clamped stiff membrane (a prestressed thin plate) of radius a:
//   sigma w_tt = T Lap w - D Lap^2 w,  w = d_r w = 0 at r = a.
// Mode (m, n): w = [J_m(al r) + C I_m(ga r)] cos(m theta) with
//   al^2 = (-T + sqrt(T^2 + 4 D sigma w^2)) / 2D, ga^2 = al^2 + T/D,
// and the clamped frequency equation (Leissa, Vibration of Plates, NASA
// SP-160, 1969, sec. 2.1 with in-plane tension; Suzuki, Yamaguchi & Izumi,
// Acoust. Sci. Tech. 30(5), 2009, eq. 4):
//   J_m(al a) ga I_m'(ga a) - al J_m'(al a) I_m(ga a) = 0.
// I_m is used scaled by e^{-ga a}, which leaves the roots unchanged.
struct StiffMode
{
    int m = 0, n = 1;
    double frequency = 0.0;   // Hz
    double alpha = 0.0, gamma = 0.0;
    double jRim = 0.0;        // J_m(al a); the shape's I-term is -jRim I_m(ga r) / I_m(ga a)
};

namespace detail
{
inline void stiffWavenumbers (const HeadSpec& head, double omega, double& al, double& ga) noexcept
{
    const double T = head.tension, D = head.bendingStiffness, s = head.arealDensity;
    const double root = std::sqrt (T * T + 4.0 * D * s * omega * omega);
    al = std::sqrt ((-T + root) / (2.0 * D));
    ga = std::sqrt ((T + root) / (2.0 * D));
}

inline double clampedDeterminant (const HeadSpec& head, double a, int m, double omega)
{
    double al, ga;
    stiffWavenumbers (head, omega, al, ga);
    std::vector<double> j, i;
    special::besselJSequence (al * a, m + 1, j);
    special::scaledBesselISequence (ga * a, m + 1, i);
    const double jm = j[static_cast<std::size_t> (m)];
    const double jmPrime = m == 0 ? -j[1] : 0.5 * (j[static_cast<std::size_t> (m - 1)] - j[static_cast<std::size_t> (m + 1)]);
    const double im = i[static_cast<std::size_t> (m)];
    const double imPrime = m == 0 ? i[1] : 0.5 * (i[static_cast<std::size_t> (m - 1)] + i[static_cast<std::size_t> (m + 1)]);
    return jm * ga * imPrime - al * jmPrime * im;
}
} // namespace detail

// All clamped modes with frequency below maxHz for orders 0..maxOrder,
// sorted by frequency within each order (n = 1, 2, ...).
[[nodiscard]] inline std::vector<StiffMode> clampedStiffModes (const HeadSpec& head, double a,
                                                                int maxOrder, double maxHz)
{
    std::vector<StiffMode> modes;
    const double omegaMax = 2.0 * special::pi * maxHz;
    const double c = std::sqrt (head.tension / head.arealDensity);
    for (int m = 0; m <= maxOrder; ++m)
    {
        // Scan omega so al a advances by < 0.1 per step (zeros of the
        // determinant are about pi apart in al a).
        const double step = 0.1 * c / a;
        double previous = detail::clampedDeterminant (head, a, m, 1.0e-3);
        double wPrev = 1.0e-3;
        int n = 0;
        for (double w = step; w <= omegaMax + step; w += step)
        {
            const double value = detail::clampedDeterminant (head, a, m, w);
            if ((value < 0.0) != (previous < 0.0))
            {
                const auto f = [&] (double x) { return detail::clampedDeterminant (head, a, m, x); };
                const double root = special::detail::refineRoot (f, wPrev, w, previous, value);
                if (root <= omegaMax)
                {
                    StiffMode mode;
                    mode.m = m;
                    mode.n = ++n;
                    mode.frequency = root / (2.0 * special::pi);
                    detail::stiffWavenumbers (head, root, mode.alpha, mode.gamma);
                    mode.jRim = special::besselJ (m, mode.alpha * a);
                    modes.push_back (mode);
                }
            }
            previous = value;
            wPrev = w;
        }
    }
    return modes;
}

// Radial shape of a clamped stiff mode at radius r (r <= a), unnormalised:
// J_m(al r) - J_m(al a) I_m(ga r) / I_m(ga a), the I ratio formed from the
// scaled functions so a near-membrane head (ga a ~ 10^3) cannot overflow.
[[nodiscard]] inline double clampedStiffShape (const StiffMode& mode, double r, double a)
{
    const double ratio = special::scaledBesselI (mode.m, mode.gamma * r) / special::scaledBesselI (mode.m, mode.gamma * a)
                         * std::exp (mode.gamma * (r - a));
    return special::besselJ (mode.m, mode.alpha * r) - mode.jRim * ratio;
}

// Rigid closed cylinder (radius a, length L) acoustic modes:
// f = c/(2 pi) sqrt((j'_mq / a)^2 + (p pi / L)^2) (Morse & Ingard,
// Theoretical Acoustics, 1968, sec. 9.2).
struct CavityMode
{
    int m = 0, q = 0, p = 0;
    double frequency = 0.0;
};

[[nodiscard]] inline std::vector<CavityMode> rigidCylinderModes (double a, double L, double c, double maxHz)
{
    std::vector<CavityMode> modes;
    const double kMax = 2.0 * special::pi * maxHz / c;
    const auto zeros = special::besselJPrimeZeros (40, kMax * a + 1.0);
    for (int m = 0; m < static_cast<int> (zeros.size()); ++m)
        for (int q = 0; q < static_cast<int> (zeros[static_cast<std::size_t> (m)].size()); ++q)
            for (int p = 0;; ++p)
            {
                const double kr = zeros[static_cast<std::size_t> (m)][static_cast<std::size_t> (q)] / a;
                const double kz = p * special::pi / L;
                const double kk = std::sqrt (kr * kr + kz * kz);
                if (kk > kMax)
                    break;
                modes.push_back ({ m, q, p, c * kk / (2.0 * special::pi) });
            }
    return modes;
}

} // namespace taikor::candidates::fdtd
