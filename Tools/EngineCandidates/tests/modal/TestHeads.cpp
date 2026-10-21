// Head-mode checks: the clamped stiff-membrane roots in both limits, and the
// closed-form integrals against brute-force quadrature.

#include "ModalMath.h"
#include "ModalTestUtil.h"

#include <algorithm>

using namespace taikor::candidates;
using namespace taikor::candidates::modal;

namespace
{
HeadMaterial nagadoBatter()
{
    // The values the harness resolves for the nagado batter (only used as a
    // representative stiff hide here; the engine reads them from DrumSpec).
    HeadMaterial h;
    h.tension = 6146.67;
    h.arealDensity = 0.94;
    h.bending = 0.266213;
    h.radius = 0.39;
    return h;
}

// Gauss-Legendre on [0, a] with n nodes.
void radialRule (double a, int n, std::vector<double>& r, std::vector<double>& w)
{
    special::gaussLegendre (n, r, w);
    for (std::size_t i = 0; i < r.size(); ++i)
    {
        r[i] = 0.5 * a * (r[i] + 1.0);
        w[i] *= 0.5 * a;
    }
}
} // namespace

// (1a) D -> 0: roots tend to the zeros of J_m, with an error of order 1/(b2 a).
MODAL_TEST (heads_membrane_limit)
{
    HeadMaterial h = nagadoBatter();
    h.bending = 1.0e-9 * h.tension * h.radius * h.radius; // b2 a ~ 3e4
    const double wMax = 2.0 * pi * 800.0;
    const auto modes = solveHeadModes (h, wMax);
    const auto zeros = special::besselJZeros (6, 40.0);
    double worst = 0.0;
    for (int m = 0; m <= 5; ++m)
        for (int n = 0; n < 4; ++n)
        {
            const double x = modes[static_cast<std::size_t> (m)][static_cast<std::size_t> (n)].b1 * h.radius;
            const double z = zeros[static_cast<std::size_t> (m)][static_cast<std::size_t> (n)];
            worst = std::max (worst, std::abs (x - z) / z);
        }
    const double b2a = std::sqrt (h.tension / h.bending) * h.radius;
    modaltest::below ("max relative root error vs J_m zeros (m<=5, n<=4), times b2 a", worst * b2a, 2.0);
}

// (1b) T -> 0: the clamped circular plate, frequency parameter lambda^2 =
// omega a^2 sqrt(sigma/D) from Leissa, Vibration of Plates, NASA SP-160 (1969),
// Table 2.1: (0,1) 10.2158, (1,1) 21.260, (2,1) 34.877, (0,2) 39.771, (3,1) 51.030.
MODAL_TEST (heads_plate_limit)
{
    HeadMaterial h = nagadoBatter();
    h.tension = 1.0e-8 * h.bending / (h.radius * h.radius);
    const auto modes = solveHeadModes (h, 2.0 * pi * 4000.0);
    const auto lambda2 = [&] (int m, int n)
    {
        const auto& mode = modes[static_cast<std::size_t> (m)][static_cast<std::size_t> (n - 1)];
        return mode.omega * h.radius * h.radius * std::sqrt (h.arealDensity / h.bending);
    };
    modaltest::nearRel ("clamped plate (0,1)", lambda2 (0, 1), 10.2158, 2e-5);
    modaltest::nearRel ("clamped plate (1,1)", lambda2 (1, 1), 21.260, 5e-5);
    modaltest::nearRel ("clamped plate (2,1)", lambda2 (2, 1), 34.877, 5e-5);
    modaltest::nearRel ("clamped plate (0,2)", lambda2 (0, 2), 39.771, 5e-5);
    modaltest::nearRel ("clamped plate (3,1)", lambda2 (3, 1), 51.030, 5e-5);
}

// (1c) Suzuki, Yamaguchi & Izumi (2009): a 2 mm cowhide (E 3.5 GPa, 2000 kg/m^3,
// nu 0.3) clamped on r = 0.23 m at T = 20,800 N/m has T a^2 / D = 429 and its
// (0,1) mode at 127 Hz (their stretched-plate calculation, quoted to 3 digits).
MODAL_TEST (heads_suzuki2009)
{
    HeadMaterial h;
    h.radius = 0.23;
    h.arealDensity = 2000.0 * 0.002;
    h.bending = 3.5e9 * 0.002 * 0.002 * 0.002 / (12.0 * (1.0 - 0.09));
    h.tension = 20800.0;
    modaltest::nearRel ("T a^2 / D", h.tension * h.radius * h.radius / h.bending, 429.0, 2e-3);
    const auto modes = solveHeadModes (h, 2.0 * pi * 400.0);
    modaltest::near ("(0,1) frequency, Hz", modes[0][0].omega / (2.0 * pi), 127.0, 0.5);
}

// Closed-form norms, volume, Hankel transform, disk average and gradient form
// against Gauss-Legendre quadrature of the shapes themselves.
MODAL_TEST (heads_closed_forms)
{
    const HeadMaterial h = nagadoBatter();
    const double a = h.radius;
    const auto modes = solveHeadModes (h, 2.0 * pi * 6000.0);
    std::vector<double> r, w;
    radialRule (a, 1600, r, w);
    double normErr = 0.0, volErr = 0.0, hankelErr = 0.0, diskErr = 0.0;
    for (int m : { 0, 1, 4, 17 })
    {
        const auto& set = modes[static_cast<std::size_t> (m)];
        for (std::size_t i = 0; i < set.size(); i += 3)
        {
            const HeadMode& mode = set[i];
            double norm = 0.0, vol = 0.0;
            for (std::size_t k = 0; k < r.size(); ++k)
            {
                const double v = evaluateShape (mode, a, r[k]).value;
                norm += v * v * r[k] * w[k];
                vol += v * r[k] * w[k];
            }
            norm *= angularNorm (m);
            normErr = std::max (normErr, std::abs (norm - mode.normSq) / mode.normSq);
            if (m == 0)
                volErr = std::max (volErr, std::abs (2.0 * pi * vol - mode.volume) / std::sqrt (mode.normSq));
            for (double kappa : { 3.0, 0.7 * mode.b1, mode.b1 * (1.0 + 3e-8), 1.4 * mode.b1, 900.0 })
            {
                double quad = 0.0;
                for (std::size_t k = 0; k < r.size(); ++k)
                    quad += evaluateShape (mode, a, r[k]).value * special::besselJ (m, kappa * r[k]) * r[k] * w[k];
                const double jm = special::besselJ (m, kappa * a), jp = special::besselJPrime (m, kappa * a);
                const double closed = hankelTransform (mode, a, kappa, jm, jp);
                hankelErr = std::max (hankelErr, std::abs (closed - quad) / std::sqrt (mode.normSq / angularNorm (m)));
            }
            // Disk average by 2-D Gauss quadrature over the footprint.
            const double rs = 0.0585, b = 0.01125;
            std::vector<double> u, uw;
            special::gaussLegendre (48, u, uw);
            double sum = 0.0;
            for (std::size_t p = 0; p < u.size(); ++p)
                for (int t = 0; t < 96; ++t)
                {
                    const double rho = 0.5 * b * (u[p] + 1.0), wr = 0.5 * b * uw[p];
                    const double th = 2.0 * pi * t / 96.0;
                    const double x = rs + rho * std::cos (th), y = rho * std::sin (th);
                    const double rr = std::hypot (x, y), ang = std::atan2 (y, x);
                    sum += evaluateShape (mode, a, rr).value * std::cos (m * ang) * rho * wr * (2.0 * pi / 96.0);
                }
            const double quadMean = sum / (pi * b * b);
            diskErr = std::max (diskErr, std::abs (quadMean - diskAverage (mode, a, rs, b)));
        }
    }
    modaltest::below ("normSq closed form, max relative error", normErr, 1e-10);
    modaltest::below ("volume closed form, max error / sqrt(norm)", volErr, 1e-10);
    modaltest::below ("Hankel transform closed form, max error / sqrt(norm)", hankelErr, 1e-9);
    modaltest::below ("disk average (mean-value theorem), max abs error", diskErr, 1e-10);
}

MODAL_TEST (heads_gradient_form)
{
    const HeadMaterial h = nagadoBatter();
    const double a = h.radius;
    const auto modes = solveHeadModes (h, 2.0 * pi * 12000.0);
    std::vector<double> r, w;
    radialRule (a, 1600, r, w);
    for (int m : { 0, 9 })
    {
        const auto& set = modes[static_cast<std::size_t> (m)];
        const auto form = gradientForm (set, h, 1e-12);
        const std::size_t n = set.size();
        std::vector<std::vector<ShapeValue>> s (n);
        for (std::size_t i = 0; i < n; ++i)
            for (double rr : r)
                s[i].push_back (evaluateShape (set[i], a, rr));
        double worst = 0.0;
        for (std::size_t i = 0; i < n; i += 2)
            for (std::size_t j = i; j < n; j += 3)
            {
                double g = 0.0;
                for (std::size_t k = 0; k < r.size(); ++k)
                    g += (s[i][k].slope * s[j][k].slope + m * m * s[i][k].value * s[j][k].value / (r[k] * r[k]))
                       * r[k] * w[k];
                g *= angularNorm (m);
                double model = i == j ? form.diagonal[i] : 0.0;
                for (const auto& l : form.lowRank)
                    model -= l[i] * l[j];
                worst = std::max (worst, std::abs (g - model) / std::sqrt (form.diagonal[i] * form.diagonal[j]));
            }
        std::printf ("  m %d: %zu modes, rank %zu, truncation %.2e\n", m, n, form.lowRank.size(), form.truncation);
        modaltest::below ("gradient form vs quadrature, max normalised error", worst, 1e-9);
    }
}
