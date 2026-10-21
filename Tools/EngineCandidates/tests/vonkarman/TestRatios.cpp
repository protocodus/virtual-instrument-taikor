// (5) Von Karman stretching energy against exact Berger (immovable edge),
// V_B = C (int |grad w|^2)^2 / (8 pi a^2), for mode shapes and a local dent.
//
// Why the ratio differs from 1. Per unit area the plane-stress energy is
//   C/2 [ (tr eps)^2 - 2 (1 - nu) det eps ],
// and with an immovable edge int tr eps dA = 1/2 int |grad w|^2 dA exactly
// (the div u term integrates to the rim flux, zero). Berger keeps only the
// first invariant and spreads it uniformly, which by Jensen's inequality is
// the least trace energy compatible with that integral: V_B = C (I/2)^2 /
// (2 A). Von Karman differs in two ways:
//  - the stretching strain stays where the slopes are, so a strain field
//    concentrated in part of the head carries more trace energy than its
//    average (ratio > 1, largest for a local dent);
//  - the second invariant det eps is positive where the stretch is nearly
//    isotropic (the centre of a bowl), and that lowers the energy below the
//    trace-only value (ratio < 1 for smooth global shapes, e.g. the exact
//    parabola 0.726 and quartic 0.943 of TestAnalytic.cpp).

#include "SpecialFunctions.h"
#include "TestUtil.h"

#include <algorithm>
#include <tuple>

namespace vktest
{

void testRatios()
{
    namespace sf = taikor::candidates::special;
    const double a = headRadius();
    const double nu = batter().poissonRatio;
    const auto zeros = sf::besselJZeros (2, 12.0);
    const double j01 = zeros[0][0], j02 = zeros[0][1], j11 = zeros[1][0], j21 = zeros[2][0];
    const auto exact01 = axisymmetricEnergy ([&] (double r) { return -j01 / a * sf::besselJ (1, j01 * r / a); }, a, nu, c11());
    const auto exact02 = axisymmetricEnergy ([&] (double r) { return -j02 / a * sf::besselJ (1, j02 * r / a); }, a, nu, c11());
    const double r01 = exact01.first / exact01.second, r02 = exact02.first / exact02.second;
    note (format ("continuum (radial quadrature): (0,1) V/V_B = %.5f, (0,2) V/V_B = %.5f", r01, r02));
    // Non-axisymmetric modes: the Fourier-decomposed continuum problem
    // (TestFourier.cpp; numerical: 800 linear elements, within 2e-6 of their
    // Richardson limit), cross-checked on (0,1) against the quadrature.
    const auto dJ = [] (int m, double x) { return 0.5 * (sf::besselJ (m - 1, x) - sf::besselJ (m + 1, x)); };
    const auto fourier = [&] (int m, double z)
    {
        const auto e = fourierEnergy ([&] (double r) { return sf::besselJ (m, z * r / a); },
                                      [&] (double r) { return m == 0 ? -z / a * sf::besselJ (1, z * r / a) : z / a * dJ (m, z * r / a); },
                                      m, a, nu, c11(), 800);
        return e.first / e.second;
    };
    const double f01 = fourier (0, j01), f11 = fourier (1, j11), f21 = fourier (2, j21);
    note (format ("continuum (Fourier-decomposed FE): (0,1) %.5f, (1,1) %.5f, (2,1) %.5f", f01, f11, f21));
    check (std::abs (f01 / r01 - 1) < 1e-5,
           format ("the two continuum methods agree on (0,1): %.6f against %.6f (< 1e-5)", f01, r01));
    const auto& s = shipping().stroke;

    double dentRatio[2] {}, ratio01[2] {}, ratio11[2] {}, ratio21[2] {};
    int index = 0;
    for (double h : { 0.0069, 0.0048 })
    {
        Fixture f = makeFixture (a, h, 1.0e-12);
        if (! f.vk)
            return check (false, "prepare: " + f.error);
        const auto& g = f.grid;
        struct Shape { const char* name; std::vector<double> w; };
        std::vector<Shape> shapes;
        shapes.push_back ({ "(0,1)", sampleField (g, [&] (double x, double y) { return sf::besselJ (0, j01 * std::hypot (x, y) / a); }) });
        shapes.push_back ({ "(1,1)", sampleField (g, [&] (double x, double y) { return sf::besselJ (1, j11 * std::hypot (x, y) / a) * std::cos (std::atan2 (y, x)); }) });
        shapes.push_back ({ "(2,1)", sampleField (g, [&] (double x, double y) { return sf::besselJ (2, j21 * std::hypot (x, y) / a) * std::cos (2 * std::atan2 (y, x)); }) });
        shapes.push_back ({ "(0,2)", sampleField (g, [&] (double x, double y) { return sf::besselJ (0, j02 * std::hypot (x, y) / a); }) });
        // Local dent: the static deflection of the linear membrane under the
        // test stroke's footprint (uniform disk of footprintRadius at the
        // strike point), -T Lap w = load.
        const auto load = sampleField (g, [&] (double x, double y)
        {
            return std::hypot (x - s.strikeX, y - s.strikeY) < s.footprintRadius ? 1.0 : 0.0;
        });
        shapes.push_back ({ "dent", solvePoisson (g, load) });
        for (auto& shape : shapes)
        {
            double peak = 0.0;
            for (double v : shape.w)
                peak = std::max (peak, std::abs (v));
            for (double& v : shape.w)
                v *= 0.001 / peak; // 1 mm peak; both energies are quartic, the ratio is amplitude-free
            const double V = f.vk->evaluate (shape.w.data(), nullptr);
            const double VB = bergerEnergy (g, shape.w.data(), a);
            // Peak of the von Karman mean membrane stress against Berger's
            // uniform tension N_B = C I / (2 pi a^2) (dV_B/dI times 2).
            const double NB = c11() * gradientSquared (g, shape.w.data()) / (2.0 * std::numbers::pi * a * a);
            double peakMean = 0.0;
            for (int k = 0; k < g.size(); ++k)
                if (g.interior[static_cast<std::size_t> (k)] != 0)
                    peakMean = std::max (peakMean, 0.5 * (f.vk->stressXX()[k] + f.vk->stressYY()[k]));
            note (format ("h = %.1f mm %-6s V_vK / V_B = %8.4f   peak (N_xx + N_yy)/2 / N_B = %7.3f", h * 1e3,
                          shape.name, V / VB, peakMean / NB));
            if (std::string (shape.name) == "dent") dentRatio[index] = V / VB;
            if (std::string (shape.name) == "(0,1)") ratio01[index] = V / VB;
            if (std::string (shape.name) == "(1,1)") ratio11[index] = V / VB;
            if (std::string (shape.name) == "(2,1)") ratio21[index] = V / VB;
        }
        ++index;
    }
    check (std::abs (ratio01[1] / r01 - 1) < 0.01,
           format ("(0,1): discrete V/V_B %.5f (6.9 mm), %.5f (4.8 mm) against the continuum %.5f (< 1 %%)", ratio01[0], ratio01[1], r01));
    // Converging (the finer grid closer) and within 1 % at 4.8 mm, as (0,1).
    for (const auto& [name, ratio, exact] : { std::tuple { "(1,1)", ratio11, f11 }, std::tuple { "(2,1)", ratio21, f21 } })
        check (std::abs (ratio[1] / exact - 1) < 0.01 && std::abs (ratio[1] - exact) < std::abs (ratio[0] - exact),
               format ("%s: discrete V/V_B %.5f (6.9 mm), %.5f (4.8 mm) converging on the continuum %.5f (< 1 %%)", name,
                       ratio[0], ratio[1], exact));
    check (dentRatio[0] > 1.0 && dentRatio[1] > 1.0,
           format ("local dent stores more stretching energy than Berger spreads over the head: V/V_B = %.3f, %.3f (> 1)",
                   dentRatio[0], dentRatio[1]));
}

} // namespace vktest
