// (2) The King integrals: the exact baffled piston (added mass 8 rho a^3/3 at
// low ka and the Rayleigh/Struve impedance at finite ka), the monopole limit of
// the (0,n) modes, and the truncation of the wavenumber integral.

#include "ModalExterior.h"
#include "ModalTestUtil.h"

#include <algorithm>

using namespace taikor::candidates;
using namespace taikor::candidates::modal;

namespace
{
// Struve H_1 by its power series (Abramowitz & Stegun 12.1.3), fine for x <= 8.
double struveH1 (double x)
{
    double sum = 0.0;
    for (int k = 0; k < 60; ++k)
    {
        const double term = std::pow (-1.0, k) * std::pow (0.5 * x, 2 * k + 2)
                          / (std::tgamma (k + 1.5) * std::tgamma (k + 2.5));
        sum += term;
        if (std::abs (term) < 1e-18 * std::abs (sum))
            break;
    }
    return sum;
}

HeadMaterial batter()
{
    HeadMaterial h;
    h.tension = 6146.67;
    h.arealDensity = 0.94;
    h.bending = 0.266213;
    h.radius = 0.39;
    return h;
}
} // namespace

MODAL_TEST (exterior_piston)
{
    const double a = 0.39, rho = 1.2041, c = 343.0;
    const auto piston = [a] (double kappa) { return kappa == 0.0 ? 0.5 * a * a : a * special::besselJ (1, kappa * a) / kappa; };
    double mass = 0.0, res = 0.0;
    kingIntegrals (0, a, 0.001 / a, 3000.0 / a, rho, c, piston, mass, res);
    modaltest::nearRel ("piston added mass at ka = 0.001 vs 8 rho a^3/3", mass, 8.0 * rho * a * a * a / 3.0, 1e-5);
    for (double ka : { 0.25, 1.0, 3.0 })
    {
        const double k = ka / a;
        kingIntegrals (0, a, k, 3000.0 / a, rho, c, piston, mass, res);
        const double x = 2.0 * ka;
        const double exactR = rho * c * pi * a * a * (1.0 - 2.0 * special::besselJ (1, x) / x);
        const double exactM = rho * pi * a * a * (2.0 * struveH1 (x) / x) / k;
        std::printf ("  ka %.2f\n", ka);
        modaltest::nearRel ("  piston resistance vs rho c S (1 - 2 J1(2ka)/2ka)", res, exactR, 1e-9);
        modaltest::nearRel ("  piston added mass vs rho S X1(2ka)/k (Struve)", mass, exactM, 2e-6);
    }
}

MODAL_TEST (exterior_monopole_limit)
{
    // At ka -> 0 an (0,n) mode radiates as a baffled monopole of its volume:
    // R -> rho c k^2 V^2 / (2 pi) (equivalently radiation efficiency 2 k^2/b1^2
    // for the ideal membrane's (0,n) with V = 2 pi a^2 J_1(j)/j scaling).
    const HeadMaterial h = batter();
    const auto modes = solveHeadModes (h, 2.0 * pi * 800.0);
    const double rho = 1.2041, c = 343.0;
    ExteriorTable table (h.radius, 12.0 * modes[0].back().b1 + 60.0 / h.radius, 4);
    for (int n = 0; n < 3; ++n)
    {
        const HeadMode& mode = modes[0][static_cast<std::size_t> (n)];
        const double k = 0.002 / h.radius;
        const double r = table.selfResistance (mode, k, rho, c);
        const double monopole = rho * c * k * k * mode.volume * mode.volume / (2.0 * pi);
        char label[96];
        std::snprintf (label, sizeof label, "(0,%d) resistance / monopole limit at ka = 0.002", n + 1);
        modaltest::nearRel (label, r / monopole, 1.0, 1e-4);
    }
}

MODAL_TEST (exterior_truncation)
{
    const HeadMaterial h = batter();
    const auto modes = solveHeadModes (h, 2.0 * pi * 4000.0);
    double b1Max = 0.0;
    for (const auto& set : modes)
        for (const auto& mode : set)
            b1Max = std::max (b1Max, mode.b1);
    const double rho = 1.2041;
    const double kappaMax = 12.0 * b1Max + 60.0 / h.radius;
    ExteriorTable table (h.radius, kappaMax, 20);
    ExteriorTable wide (h.radius, 2.0 * kappaMax, 20);
    double worst = 0.0;
    for (int m : { 0, 1, 7, 20 })
        for (std::size_t i = 0; i < modes[static_cast<std::size_t> (m)].size(); i += 4)
        {
            const HeadMode& mode = modes[static_cast<std::size_t> (m)][i];
            const double k = mode.omega / 343.0;
            const double a1 = table.selfAddedMass (mode, k, rho), a2 = wide.selfAddedMass (mode, k, rho);
            worst = std::max (worst, std::abs (a1 - a2) / a2);
        }
    modaltest::below ("self added mass, relative change when kappaMax doubles", worst, 1e-8);
    // Air-loading ratio of (0,1) and (1,1) for orientation (research note: 0.36 / 0.18).
    const double k01 = modes[0][0].omega / 343.0, k11 = modes[1][0].omega / 343.0;
    modaltest::note ("(0,1) exterior added mass / hide mass = %.4f",
                     table.selfAddedMass (modes[0][0], k01, rho) / (h.arealDensity * modes[0][0].normSq));
    modaltest::note ("(1,1) exterior added mass / hide mass = %.4f",
                     table.selfAddedMass (modes[1][0], k11, rho) / (h.arealDensity * modes[1][0].normSq));
    // Mutual matrix: symmetric by construction; resistance must be PSD.
    std::vector<const HeadMode*> set;
    for (std::size_t i = 0; i < 6; ++i)
        set.push_back (&modes[0][i]);
    std::vector<double> res;
    table.resistance (set, 2.0 * pi * 300.0 / 343.0, rho, 343.0, res);
    double minDiag = 1e300, maxOffRatio = 0.0;
    for (std::size_t i = 0; i < 6; ++i)
    {
        minDiag = std::min (minDiag, res[i * 6 + i]);
        for (std::size_t j = 0; j < 6; ++j)
            if (i != j)
                maxOffRatio = std::max (maxOffRatio, res[i * 6 + j] * res[i * 6 + j] / (res[i * 6 + i] * res[j * 6 + j]));
    }
    modaltest::check (minDiag > 0.0 && maxOffRatio <= 1.0 + 1e-12, "resistance matrix: positive diagonal, |R_ij|^2 <= R_ii R_jj",
                      maxOffRatio, 1.0, 0.0);
}
