#pragma once

// The finite-difference drum's purely numerical parameters (grid spacings,
// rates, absorbing-layer size, solver truncation). None of them is physical:
// each is set by a stability bound or justified by a convergence or accuracy
// test in tests/fdtd/, cited beside its value. Tests may override them through
// the environment (TAIKOR_FDTD_<NAME>), the production render never does.

#include "CandidateEngine.h"
#include "HeadGrid.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>

namespace taikor::candidates::fdtd
{

struct FdtdNumerics
{
    // Internal rate = oversample * 48 kHz. The head grid spacing sits at its
    // stability bound (headSpacing below), which the bending term dominates
    // (h ~ 2 sqrt(kappa k)), so the heads' dispersion error falls only as
    // 1/fs. At 96 kHz test 1 measures the nagado heads' discrete modes (after
    // the staircase correction, at the production spacing) within 0.39 % of
    // the analytic clamped values below 1 kHz and within 1.48 % up to 2 kHz;
    // the dispersion relation puts 48 kHz at about twice those errors (~1.1 %
    // at 1 kHz, ~3.8 % at 2 kHz before correction).
    int oversample = 2;

    // Air cells across the barrel length L (Delta = L / airCellsPerDepth, so
    // both heads fall midway between air planes; Bilbao, JASA 131(1), 2012).
    // At Delta = 12.48 mm the rigid-cylinder axial modes come out within
    // 0.06 % and the transverse ones 0.4-1.0 % flat (the staircase barrel;
    // test 2). The air grid's own cut-off is asin(lambda) / (pi k) = 8.87 kHz
    // on the axes and 15.8 kHz on the diagonal at 96 kHz (Courant lambda
    // 0.286), and its phase error reaches ~1 % at 2 kHz and ~3 % at 4 kHz on
    // the axes. Test 8 (heavy) compares the raw left-capsule pressure over the
    // first 40 ms of the Don against finer air: L / 100 moves the octaves
    // below 2 kHz by at most 0.15 dB and the 1 kHz bands 1-5 kHz by -0.03,
    // +0.11, -0.38, +1.02 dB; L / 125 by 0.22 dB and -0.04, +0.24, +0.10,
    // +1.29 dB. So in air spacing the take is converged to ~0.2 dB below
    // 2 kHz and to ~1.3 dB up to the output FIR's 4.9 kHz pass edge; above
    // 2 kHz the heads' dispersion dominates (twice the rate: +1.60, +0.73,
    // +2.71 dB in 2-5 kHz, see oversample). The contact's impulse moves by
    // 0.05 % (0.28873 / 0.28865 / 0.28860 N s) but its peak force by -3.7 %
    // (118.9 / 115.6 / 114.5 N): the peak rides the ~9 kHz footprint
    // oscillation (Decimator.h), whose air loading this grid cannot carry on
    // its axes. Finer air costs Delta^-3 (L / 150, the Courant limit at this
    // rate, 8 times the cells) and needs a different absorbing layer
    // (TestConvergence.cpp).
    int airCellsPerDepth = 75;

    // Air between the drum's outer surface and the absorbing layer, m, beside
    // and below the drum and above the batter (the capsules sit 0.16 m above
    // it, so the layer stays 7 cells clear of them), and the layer itself
    // (Grote & Sim PML): thickness in cells and its theoretical
    // normal-incidence round-trip reflection. Test 5 measures the reflection
    // and the effect of doubling the margins.
    // 12 cells with R0 = 1e-3 reflect a 2-cell (4-cell) Gaussian pulse at
    // -53.8 (-40.4) dB peak, worst over normal and 45/55-degree incidence;
    // 10 cells manage only -47.8 (-36.9) dB and 14 cells -59.3 (-43.8) dB
    // (test 5; the oblique corner path dominates, R0^cos(theta)). Doubling
    // both margins moves the partials below 400 Hz by at most 0.03 % and
    // their decay times by -0.7 ... +1.3 % (8 partials of the Don over its
    // first second, estimation-limited; test 5). The raw left-capsule
    // pressure itself moves by -44 ... -36 dB in 40-400 Hz over that second,
    // but by -23 dB in 2-8 kHz over 20-100 ms (-28 dB to 0.4 s, -36 dB to
    // 1 s): the layer's reflection of the barrel's radiated high band, 7
    // cells above the capsules.
    double airMargin = 0.15;
    double airMarginAbove = 0.25;
    int pmlCells = 12;
    double pmlReflection = 1.0e-3;

    // Terms of the Neumann series for the head-air coupling solve in the
    // air-slice space: (I + nu G)^-1 ~ sum_{p=0..order} (-nu G)^p, nu ~ 6e-4
    // at 96 kHz, so order 2 leaves a relative residual ~ nu^3 ~ 3e-10.
    int neumannOrder = 2;

    // Float32 air in production: against double air the output differs by
    // -82 dB over the first 5 ms, -78 dB over 15 ms and -58.6 dB over 40 ms
    // (energy; test 9, where the mechanism is described: rounding of the
    // stored potential at Courant 0.29, a random walk in the partials' phase),
    // once the sealed interior's potential offset is removed every 64 steps
    // (FdtdEngineStep.cpp); double air costs about five times as much.
    bool doubleAir = false;

    // SAV regulariser, J: psi = sqrt(2 V + eps), g = grad V / sqrt(2 V + eps).
    // It only keeps g defined at V = 0, where grad V vanishes too (g stays
    // finite as eps -> 0). Its footprint on the take is the SAV's lag at the
    // onset, while V passes eps within the first steps, and it scales as
    // sqrt(eps): against 1e-18 J, 1e-12 J moved the batter by 6e-8 (relative
    // L2 after 8 ms of the Don), and halving 1e-18 J moves it by 4e-10.
    // Numerical; test 3 halves it.
    double savEpsilon = 1.0e-18;

    static FdtdNumerics fromEnvironment()
    {
        FdtdNumerics n;
        const auto read = [] (const char* name, auto& value)
        {
            if (const char* text = std::getenv (name))
                value = static_cast<std::decay_t<decltype (value)>> (std::atof (text));
        };
        read ("TAIKOR_FDTD_OVERSAMPLE", n.oversample);
        read ("TAIKOR_FDTD_AIR_CELLS", n.airCellsPerDepth);
        read ("TAIKOR_FDTD_AIR_MARGIN", n.airMargin);
        read ("TAIKOR_FDTD_AIR_MARGIN_ABOVE", n.airMarginAbove);
        read ("TAIKOR_FDTD_PML_CELLS", n.pmlCells);
        read ("TAIKOR_FDTD_PML_REFLECTION", n.pmlReflection);
        read ("TAIKOR_FDTD_NEUMANN_ORDER", n.neumannOrder);
        read ("TAIKOR_FDTD_DOUBLE_AIR", n.doubleAir);
        read ("TAIKOR_FDTD_SAV_EPSILON", n.savEpsilon);
        return n;
    }
};

// The energy a stroke brings into the drum at rest: the tip mass's kinetic
// energy, 1/2 m v0^2 (0.4325 J for the nagado Don).
[[nodiscard]] inline double strokeEnergy (const StrokeSpec& stroke) noexcept
{
    return 0.5 * stroke.tipMass * stroke.impactSpeed * stroke.impactSpeed;
}

// The largest membrane tension the SAV step can apply to a head when the
// whole drum holds at most `energy` joules. The SAV force is -g psi with
// g = grad V / sqrt(2 V + eps); for Berger (BergerPotential.cpp) that is the
// explicit tension T_B = sqrt(2 xi) psi acting on L w^n, xi = E h /
// (2 (1 - nu^2) pi a^2). The scheme conserves heads + air + stick + contact +
// psi^2 / 2, every term non-negative, so psi^2 / 2 <= energy and
//   |T_B| <= sqrt(2 xi) sqrt(2 energy) = 2 sqrt(xi energy)
// (2558 N/m on the nagado batter, 2453 N/m on the rear, against the ~330 N/m
// the Don reaches; the exact Berger tension xi int |grad w|^2 obeys the same
// bound through V = xi / 4 (int |grad w|^2)^2 <= energy).
[[nodiscard]] inline double bergerTensionBound (const HeadSpec& head, double radius, double energy) noexcept
{
    const double nu = head.poissonRatio;
    const double xi = head.youngsModulus * head.thickness
                      / (2.0 * (1.0 - nu * nu) * 3.14159265358979323846 * radius * radius);
    return 2.0 * std::sqrt (xi * std::max (0.0, energy));
}

// The head grid spacing: the stability bound at the stiffened tension
// T + extraTension, widened by one part in a million so rounding cannot put
// the scheme on the wrong side of it. The engine passes the bound above for
// the stroke's energy. Why the margin is needed: at the linear bound the
// grid's top (checkerboard, time-Nyquist) mode carries almost no discrete
// energy, so the SAV step's explicit tension - which acts on that mode like
// extra tension beyond the bound - can grow it at almost no energy cost.
// With h at the linear bound the Don grew such a component to ~40 um within
// 2.5 ms, holding 40-50 % of int |grad w|^2, and psi fell to 0.04-0.4 of
// sqrt(2 V(w)): the batter felt an SAV tension of 15-40 N/m over 3-12 ms
// where the correctly tracked Berger tension is 55-360 N/m. At the
// stiffened bound (+0.66 % in h) psi tracks sqrt(2 V) to ~1e-3 and the
// time-Nyquist share stays below 1e-4 (test savTracking, TestSav.cpp). By
// the frozen-coefficient argument the scheme with explicit tension T + T_B
// is within its bound for every |T_B| <= extraTension.
[[nodiscard]] inline double headSpacing (const HeadSpec& head, double k, double extraTension = 0.0) noexcept
{
    return stableHeadSpacing (head, k, extraTension) * (1.0 + 1.0e-6);
}

} // namespace taikor::candidates::fdtd
