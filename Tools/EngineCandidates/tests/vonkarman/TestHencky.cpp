// (3) Hencky's problem: a clamped circular membrane with no pre-tension and an
// immovable edge under uniform lateral pressure p. The exact solution has
//   w(0) = c a (p a / (E h))^(1/3),   N_r(0) = (b0 / 4) (E h)^(1/3) (p a)^(2/3),
// and for nu = 0.3 the shooting solution of Hencky's radial equations
// (d(r N_r)/dr = N_t, u' + w'^2/2 = (N_r - nu N_t)/(E h), u/r = (N_t - nu N_r)
// /(E h), N_r w' = -p r/2, u(0) = u(a) = 0; solved to 1e-12 with DOP853 and
// Brent) gives c = 0.653443 and b0 = 1.72440. The b0 agrees with Fichter's
// corrected power-series value 1.7244 (W. B. Fichter, "Some solutions for the
// large deflections of uniformly loaded circular membranes", NASA TP-3658,
// 1997); the 0.662 often quoted from Timoshenko & Woinowsky-Krieger inherits
// Hencky's erroneous b0 = 1.713.
//
// The discrete equilibrium dV_h/dw = p h^2 is found by a damped fixed-point
// iteration on the shape with the exact optimal amplitude: phi solves
// G(N(w)) phi = p h^2 (G the frozen-stress operator, SPD for a tensile N), the
// amplitude of a shape follows from the quartic energy, alpha^3 =
// (f . phi) / (4 V(phi)), and the shape is relaxed by 1/2 per step.

#include "TestUtil.h"

#include <algorithm>

namespace vktest
{

namespace
{
constexpr double henckyC = 0.6534433447;  // w(0) / (a (p a / E h)^(1/3)), nu = 0.3
constexpr double henckyB0 = 1.7244040646; // 4 N_r(0) / ((E h)^(1/3) (p a)^(2/3)), nu = 0.3

// CG on the frozen-stress operator, to a relative residual of 1e-12.
std::vector<double> solveStressOperator (VonKarmanPotential& vk, const std::vector<double>& f)
{
    const auto size = f.size();
    std::vector<double> x (size, 0.0), r = f, p = f, q (size, 0.0);
    double rr = 0.0;
    for (double v : r)
        rr += v * v;
    const double bb = rr;
    for (int it = 0; it < 50000 && rr > 1e-24 * bb; ++it)
    {
        vk.applyStressOperator (p.data(), q.data());
        double pq = 0.0;
        for (std::size_t k = 0; k < size; ++k)
            pq += p[k] * q[k];
        const double alpha = rr / pq;
        double rrNew = 0.0;
        for (std::size_t k = 0; k < size; ++k)
        {
            x[k] += alpha * p[k];
            r[k] -= alpha * q[k];
            rrNew += r[k] * r[k];
        }
        for (std::size_t k = 0; k < size; ++k)
            p[k] = r[k] + rrNew / rr * p[k];
        rr = rrNew;
    }
    return x;
}

double dotv (const std::vector<double>& a, const std::vector<double>& b)
{
    double s = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k)
        s += a[k] * b[k];
    return s;
}
} // namespace

void testHencky()
{
    const double a = headRadius();
    const auto& head = batter();
    const double Eh = head.youngsModulus * head.thickness;
    // The discrete problem is exactly self-similar in p (w ~ p^(1/3)), so p
    // only sets the scale: 200 Pa gives w(0)/a ~ 2 %, the moderate-rotation
    // range of the test stroke.
    const double p = 200.0;
    // The shipping skin stores nu in single precision (0.3f).
    check (std::abs (head.poissonRatio - 0.3) < 1e-6, "batter nu = 0.3 (the Hencky coefficients are for nu = 0.3)");
    double errEff[5] {}, errRaw[5] {}, errStress[5] {};
    int index = 0;
    for (int ah : { 20, 28, 40, 56, 80 })
    {
        Fixture fx = makeFixture (a, a / ah, 1.0e-12);
        if (! fx.vk)
            return check (false, "prepare: " + fx.error);
        auto& vk = *fx.vk;
        const auto& g = fx.grid;
        const auto f = sampleField (g, [&] (double, double) { return p * g.h * g.h; });
        const double ref = henckyC * a * std::cbrt (p * a / Eh);
        auto w = sampleField (g, [&] (double x, double y) { return ref * (1 - (x * x + y * y) / (a * a)); });
        std::vector<double> grad (w.size());
        auto rescale = [&] (std::vector<double>& shape)
        {
            const double V = vk.evaluate (shape.data(), nullptr);
            const double alpha = std::cbrt (dotv (f, shape) / (4.0 * V));
            for (double& v : shape)
                v *= alpha;
        };
        double residual = 1.0;
        int it = 0;
        for (; it < 200; ++it)
        {
            vk.evaluate (w.data(), grad.data());
            double num = 0.0;
            for (std::size_t k = 0; k < w.size(); ++k)
                num += (grad[k] - f[k]) * (grad[k] - f[k]);
            residual = std::sqrt (num / dotv (f, f));
            if (residual < 1e-9)
                break;
            auto phi = solveStressOperator (vk, f);
            rescale (phi);
            for (std::size_t k = 0; k < w.size(); ++k)
                w[k] = 0.5 * (w[k] + phi[k]);
            rescale (w);
        }
        vk.evaluate (w.data(), grad.data());
        const int centre = g.index ((g.nx - 1) / 2, (g.ny - 1) / 2);
        const double w0 = w[static_cast<std::size_t> (centre)];
        const double n0 = 0.5 * (vk.stressXX()[centre] + vk.stressYY()[centre]);
        // Effective radius of the staircase: -Lap_h w = 1 is solved exactly by
        // (a^2 - r^2)/4 away from the rim, so its centre value is a_eff^2 / 4.
        const auto ones = sampleField (g, [] (double, double) { return 1.0; });
        const double aEff = std::sqrt (4.0 * solvePoisson (g, ones)[static_cast<std::size_t> (centre)]);
        const double refEff = henckyC * aEff * std::cbrt (p * aEff / Eh);
        const double n0Ref = 0.25 * henckyB0 * std::cbrt (Eh) * std::pow (p * aEff, 2.0 / 3.0);
        errRaw[index] = w0 / ref - 1;
        errEff[index] = w0 / refEff - 1;
        errStress[index] = n0 / n0Ref - 1;
        note (format ("a/h = %2d: %3d its (residual %.1e), w(0) = %.4f mm, Hencky %.4f mm -> err %+.3e; staircase a_eff/a - 1 = %+.2e -> err %+.3e; N(0) err %+.3e",
                      ah, it, residual, w0 * 1e3, ref * 1e3, errRaw[index], aEff / a - 1, errEff[index], errStress[index]));
        check (residual < 1e-9, format ("a/h = %d: discrete equilibrium dV/dw = p h^2 reached (%.1e)", ah, residual));
        ++index;
    }
    // The error is first order in h (a nonzero edge slope meets the O(h)
    // staircase rim, as for the parabola of TestAnalytic.cpp): err * (a/h)
    // is constant, and Richardson extrapolation of a/h = 56 and 80 removes it.
    const double extrapolatedW = (80.0 * errEff[4] - 56.0 * errEff[3]) / 24.0;
    const double extrapolatedN = (80.0 * errStress[4] - 56.0 * errStress[3]) / 24.0;
    note (format ("err x a/h, w(0) at a_eff: %.3f %.3f %.3f %.3f %.3f; N(0): %.3f %.3f %.3f %.3f %.3f",
                  20 * errEff[0], 28 * errEff[1], 40 * errEff[2], 56 * errEff[3], 80 * errEff[4],
                  20 * errStress[0], 28 * errStress[1], 40 * errStress[2], 56 * errStress[3], 80 * errStress[4]));
    check (std::abs (errEff[3]) < 0.01 && std::abs (errEff[4]) < 0.01 && std::abs (errEff[4]) < std::abs (errEff[3]),
           format ("w(0) against Hencky (c = 0.65344) at the staircase radius: %+.2e (a/h = 56), %+.2e (a/h = 80), converging (< 1 %%)",
                   errEff[3], errEff[4]));
    check (std::abs (extrapolatedW) < 1e-3, format ("w(0), Richardson-extrapolated to h -> 0: error %+.1e (< 1e-3)", extrapolatedW));
    check (std::abs (errStress[3]) < 0.01 && std::abs (errStress[4]) < 0.01,
           format ("centre stress against Hencky/Fichter (b0 = 1.7244): %+.2e (a/h = 56), %+.2e (a/h = 80) (< 1 %%)",
                   errStress[3], errStress[4]));
    check (std::abs (extrapolatedN) < 1e-3, format ("N(0), Richardson-extrapolated to h -> 0: error %+.1e (< 1e-3)", extrapolatedN));
    (void) errRaw;
}

} // namespace vktest
