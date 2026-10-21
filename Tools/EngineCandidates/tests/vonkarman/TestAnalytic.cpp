// (4) The discrete energy converges to exact continuum values with an
// immovable edge, and (5) the von Karman energy against exact Berger for
// modal and local shapes.
//
// Exact continuum energies for axisymmetric w(r) with u(a) = 0. Radial
// equilibrium d(r N_r)/dr = N_theta with N_r = C (eps_r + nu eps_t),
// N_t = C (eps_t + nu eps_r), eps_r = u' + w'^2 / 2, eps_t = u / r reduces to
//   d/dr [ (1/r) d(r u)/dr ] = -1/2 [ (1 - nu) w'^2 / r + (w'^2)' ],
// integrated twice with u(0) = 0 and u(a) = 0 (sympy, checked by
// substitution). With C = E h / (1 - nu^2):
//   w = w0 (1 - r^2/a^2):   u = (3 - nu) w0^2 r (a^2 - r^2) / (4 a^4),
//                           V = pi C w0^4 (1 + nu)(7 - nu) / (24 a^2),
//                           V / V_Berger = (1 + nu)(7 - nu) / 12;
//   w = w0 (1 - r^2/a^2)^2: V = pi C w0^4 (1 + nu)(23 - 9 nu) / (126 a^2),
//                           V / V_Berger = (1 + nu)(23 - 9 nu) / 28;
// with V_Berger = C (int |grad w|^2)^2 / (8 pi a^2). For any other
// axisymmetric w the same double integral is evaluated by quadrature below.

#include "TestUtil.h"

#include <algorithm>

namespace vktest
{

// Continuum V and V_Berger of an axisymmetric field from its slope w'(r), by
// the double integral above on a fine uniform radial grid (trapezoid).
std::pair<double, double> axisymmetricEnergy (const std::function<double (double)>& slope, double a, double nu, double C)
{
    const int M = 200000;
    const double dr = a / M;
    std::vector<double> r (M + 1), s2 (M + 1), F (M + 1), G (M + 1), H (M + 1);
    for (int i = 0; i <= M; ++i)
    {
        r[static_cast<std::size_t> (i)] = i * dr;
        const double s = slope (i * dr);
        s2[static_cast<std::size_t> (i)] = s * s;
    }
    for (int i = 0; i <= M; ++i)
    {
        const auto ii = static_cast<std::size_t> (i);
        const double ds2 = i == 0 ? (s2[1] - s2[0]) / dr
                         : i == M ? (s2[ii] - s2[ii - 1]) / dr
                                  : (s2[ii + 1] - s2[ii - 1]) / (2 * dr);
        const double s2r = i == 0 ? 0.0 : s2[ii] / r[ii]; // w'(0) = 0
        F[ii] = -0.5 * ((1.0 - nu) * s2r + ds2);
    }
    G[0] = 0.0; // G = int_0^r F
    for (int i = 1; i <= M; ++i)
        G[static_cast<std::size_t> (i)] = G[static_cast<std::size_t> (i - 1)] + 0.5 * dr * (F[static_cast<std::size_t> (i - 1)] + F[static_cast<std::size_t> (i)]);
    // r u = int_0^r r' (G + c1) dr' = H(r) + c1 r^2 / 2, u(a) = 0 fixes c1.
    H[0] = 0.0;
    for (int i = 1; i <= M; ++i)
    {
        const auto ii = static_cast<std::size_t> (i);
        H[ii] = H[ii - 1] + 0.5 * dr * (r[ii - 1] * G[ii - 1] + r[ii] * G[ii]);
    }
    const double c1 = -2.0 * H[static_cast<std::size_t> (M)] / (a * a);
    double V = 0.0, I = 0.0;
    for (int i = 1; i <= M; ++i)
    {
        const auto ii = static_cast<std::size_t> (i);
        // u = H / r + c1 r / 2, u' = G + c1 - u / r (from (r u)' = r (G + c1)).
        const double u = H[ii] / r[ii] + 0.5 * c1 * r[ii];
        const double du = G[ii] + c1 - u / r[ii];
        const double er = du + 0.5 * s2[ii];
        const double et = u / r[ii];
        const double density = 0.5 * C * (er * er + et * et + 2.0 * nu * er * et);
        const double weight = i == M ? 0.5 : 1.0;
        V += weight * density * 2.0 * std::numbers::pi * r[ii] * dr;
        I += weight * s2[ii] * 2.0 * std::numbers::pi * r[ii] * dr;
    }
    const double VB = C * I * I / (8.0 * std::numbers::pi * a * a);
    return { V, VB };
}

void testAnalytic()
{
    const double a = headRadius();
    const double nu = batter().poissonRatio;
    const double C = c11();
    const double w0 = 0.01;
    const double pi = std::numbers::pi;
    const double Vp = pi * C * std::pow (w0, 4) * (1 + nu) * (7 - nu) / (24 * a * a);
    const double Vq = pi * C * std::pow (w0, 4) * (1 + nu) * (23 - 9 * nu) / (126 * a * a);
    const double Rp = (1 + nu) * (7 - nu) / 12.0, Rq = (1 + nu) * (23 - 9 * nu) / 28.0;
    // The quadrature must reproduce the closed forms before it is trusted
    // for other shapes.
    const auto qp = axisymmetricEnergy ([&] (double r) { return -2 * w0 * r / (a * a); }, a, nu, C);
    const auto qq = axisymmetricEnergy ([&] (double r) { return -4 * w0 * r / (a * a) * (1 - r * r / (a * a)); }, a, nu, C);
    check (std::abs (qp.first / Vp - 1) < 1e-8 && std::abs (qq.first / Vq - 1) < 1e-8,
           format ("radial quadrature reproduces the closed forms (%.1e, %.1e)", qp.first / Vp - 1, qq.first / Vq - 1));
    note (format ("exact: parabola V = %.6e J (V/V_B = %.6f), quartic V = %.6e J (V/V_B = %.6f)", Vp, Rp, Vq, Rq));

    double parab56 = 0, parab80 = 0, parab40 = 0, quartWorst28 = 0;
    for (int ah : { 10, 14, 20, 28, 40, 56, 80 })
    {
        Fixture f = makeFixture (a, a / ah, 1.0e-12);
        if (! f.vk)
            return check (false, "prepare: " + f.error);
        const auto wp = sampleField (f.grid, [&] (double x, double y) { return w0 * (1 - (x * x + y * y) / (a * a)); });
        const auto wq = sampleField (f.grid, [&] (double x, double y) { const double b = 1 - (x * x + y * y) / (a * a); return w0 * b * b; });
        const double vp = f.vk->evaluate (wp.data(), nullptr);
        const double vq = f.vk->evaluate (wq.data(), nullptr);
        const double ep = vp / Vp - 1, eq = vq / Vq - 1;
        const double rp = vp / bergerEnergy (f.grid, wp.data(), a) / Rp - 1;
        const double rq = vq / bergerEnergy (f.grid, wq.data(), a) / Rq - 1;
        note (format ("a/h = %2d (h = %5.2f mm): parabola V err %+.3e (ratio err %+.3e) | quartic V err %+.3e (ratio err %+.3e)",
                      ah, a / ah * 1e3, ep, rp, eq, rq));
        if (ah == 40) parab40 = std::abs (ep);
        if (ah == 56) parab56 = std::abs (ep);
        if (ah == 80) parab80 = std::abs (ep);
        if (ah >= 28) quartWorst28 = std::max (quartWorst28, std::abs (eq));
    }
    check (quartWorst28 < 2e-3, format ("quartic (smooth edge): |V_h / V - 1| <= %.2e for a/h >= 28 (< 2e-3)", quartWorst28));
    check (parab56 < 0.035 && parab80 < 0.025 && parab80 < 0.6 * parab40,
           format ("parabola (kinked edge): first-order convergence, |err| %.3e at a/h = 40, %.3e at 56, %.3e at 80", parab40, parab56, parab80));
}

} // namespace vktest
