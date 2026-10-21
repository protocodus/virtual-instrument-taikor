// Continuum von Karman stretching energy, immovable edge, of a
// non-axisymmetric field w = f(r) cos(m theta), for the ratios test.
//
// With w_r = f' cos(m th) and w_th / r = -g sin(m th), g = m f / r, the
// quadratic strain 1/2 grad w (x) grad w has angular orders 0 and n = 2m:
//   order 0: q_rr = f'^2 / 4, q_thth = g^2 / 4 (no shear);
//   order n: q_rr = f'^2 / 4 cos n th, q_thth = -g^2 / 4 cos n th,
//            gamma = -f' g / 2 sin n th.
// The plane-stress energy is quadratic in eps = sym grad u + q and the
// orders are orthogonal over theta, so it splits into an order-0 problem in
// u_r = U0(r) (weight 2 pi) and an order-n problem in u_r = U cos n th,
// u_th = W sin n th (weight pi) with
//   eps_rr = U' + q_rr, eps_thth = (U + n W) / r + q_thth,
//   gamma  = W' - W / r - n U / r + gamma_q.
// Each is minimised over u by linear finite elements on r in [0, a] with
// u = 0 at r = a (immovable edge) and at r = 0 (regularity; U0(0) = 0 by
// symmetry), 4-point Gauss per element; the minimum of 1/2 d^T K d + b^T d
// + c is c + b^T d / 2 at K d = -b (banded Cholesky). m = 0 reduces to the
// axisymmetric problem of TestAnalytic.cpp and reproduces it.

#include "TestUtil.h"

#include <array>

namespace vktest
{
namespace
{
// Minimum over the in-plane field of one angular order: n = 0 (U only) or
// n > 0 (U and W interleaved per node). q holds {q_rr, q_thth, gamma_q} at
// every quadrature point.
double orderEnergy (int n, double weight, const std::vector<double>& rs, const std::vector<double>& ws,
                    const std::vector<int>& element, const std::vector<std::array<double, 3>>& q, int nel, double a,
                    double nu, double C)
{
    const int fields = n == 0 ? 1 : 2, dofs = fields * (nel + 1), band = 2 * fields - 1;
    std::vector<double> K (static_cast<std::size_t> (dofs * (band + 1)), 0.0), b (static_cast<std::size_t> (dofs), 0.0);
    const auto at = [&] (int i, int j) -> double& { return K[static_cast<std::size_t> (i * (band + 1) + (j - i))]; }; // j >= i
    const double Cm[3][3] = { { C, C * nu, 0 }, { C * nu, C, 0 }, { 0, 0, C * (1 - nu) / 2 } };
    double c = 0.0;
    const double dr = a / nel;
    for (std::size_t k = 0; k < rs.size(); ++k)
    {
        const int e = element[k];
        const double r = rs[k], r0 = e * dr;
        const double N[2] = { 1.0 - (r - r0) / dr, (r - r0) / dr }, dN[2] = { -1.0 / dr, 1.0 / dr };
        // B maps the element's dofs {U_e, U_e+1, W_e, W_e+1} to the strain.
        double B[3][4] = {};
        for (int s = 0; s < 2; ++s)
        {
            B[0][s] = dN[s];
            B[1][s] = N[s] / r;
            if (n > 0)
            {
                B[1][2 + s] = n * N[s] / r;
                B[2][s] = -n * N[s] / r;
                B[2][2 + s] = dN[s] - N[s] / r;
            }
        }
        const double wt = weight * r * ws[k];
        const auto& qk = q[k];
        double Cq[3] = {};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                Cq[i] += Cm[i][j] * qk[static_cast<std::size_t> (j)];
        c += 0.5 * wt * (qk[0] * Cq[0] + qk[1] * Cq[1] + qk[2] * Cq[2]);
        const int local = 2 * fields;
        int global[4];
        for (int s = 0; s < 2; ++s)
        {
            global[s] = fields * (e + s);
            if (n > 0)
                global[2 + s] = fields * (e + s) + 1;
        }
        for (int i = 0; i < local; ++i)
        {
            double BCq = 0.0;
            for (int p = 0; p < 3; ++p)
                BCq += B[p][i] * Cq[p];
            b[static_cast<std::size_t> (global[i])] += wt * BCq;
            for (int j = 0; j < local; ++j)
            {
                if (global[j] < global[i])
                    continue;
                double v = 0.0;
                for (int p = 0; p < 3; ++p)
                    for (int s = 0; s < 3; ++s)
                        v += B[p][i] * Cm[p][s] * B[s][j];
                at (global[i], global[j]) += wt * v;
            }
        }
    }
    // u = 0 at r = 0 and r = a: identity rows, zero right-hand side.
    for (int f = 0; f < fields; ++f)
        for (const int node : { 0, nel })
        {
            const int i = fields * node + f;
            for (int j = std::max (0, i - band); j <= std::min (dofs - 1, i + band); ++j)
                (j >= i ? at (i, j) : at (j, i)) = 0.0;
            at (i, i) = 1.0;
            b[static_cast<std::size_t> (i)] = 0.0;
        }
    // Banded Cholesky K = L L^T (upper storage holds L^T), then K d = -b.
    for (int i = 0; i < dofs; ++i)
    {
        for (int k = std::max (0, i - band); k < i; ++k)
        {
            const double lki = at (k, i);
            for (int j = i; j <= std::min (dofs - 1, k + band); ++j)
                at (i, j) -= lki * at (k, j);
        }
        const double d = std::sqrt (at (i, i));
        for (int j = i; j <= std::min (dofs - 1, i + band); ++j)
            at (i, j) /= d;
    }
    std::vector<double> d (static_cast<std::size_t> (dofs));
    for (int i = 0; i < dofs; ++i) // L y = -b
    {
        double s = -b[static_cast<std::size_t> (i)];
        for (int k = std::max (0, i - band); k < i; ++k)
            s -= at (k, i) * d[static_cast<std::size_t> (k)];
        d[static_cast<std::size_t> (i)] = s / at (i, i);
    }
    for (int i = dofs - 1; i >= 0; --i) // L^T d = y
    {
        double s = d[static_cast<std::size_t> (i)];
        for (int j = i + 1; j <= std::min (dofs - 1, i + band); ++j)
            s -= at (i, j) * d[static_cast<std::size_t> (j)];
        d[static_cast<std::size_t> (i)] = s / at (i, i);
    }
    double bd = 0.0;
    for (int i = 0; i < dofs; ++i)
        bd += b[static_cast<std::size_t> (i)] * d[static_cast<std::size_t> (i)];
    return c + 0.5 * bd;
}
} // namespace

std::pair<double, double> fourierEnergy (const std::function<double (double)>& f, const std::function<double (double)>& slope,
                                         int m, double a, double nu, double C, int nel)
{
    static constexpr double gx[4] = { -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526 };
    static constexpr double gw[4] = { 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538 };
    std::vector<double> rs, ws;
    std::vector<int> element;
    const double dr = a / nel;
    for (int e = 0; e < nel; ++e)
        for (int p = 0; p < 4; ++p)
        {
            rs.push_back ((e + 0.5) * dr + 0.5 * dr * gx[p]);
            ws.push_back (0.5 * dr * gw[p]);
            element.push_back (e);
        }
    std::vector<std::array<double, 3>> q0 (rs.size()), qn (rs.size());
    double I = 0.0;
    for (std::size_t k = 0; k < rs.size(); ++k)
    {
        const double fp = slope (rs[k]), g = m * f (rs[k]) / rs[k];
        I += (m == 0 ? 2.0 : 1.0) * std::numbers::pi * ws[k] * (fp * fp + g * g) * rs[k];
        if (m == 0)
            q0[k] = { 0.5 * fp * fp, 0.0, 0.0 };
        else
        {
            q0[k] = { 0.25 * fp * fp, 0.25 * g * g, 0.0 };
            qn[k] = { 0.25 * fp * fp, -0.25 * g * g, -0.5 * fp * g };
        }
    }
    double V = orderEnergy (0, 2.0 * std::numbers::pi, rs, ws, element, q0, nel, a, nu, C);
    if (m > 0)
        V += orderEnergy (2 * m, std::numbers::pi, rs, ws, element, qn, nel, a, nu, C);
    return { V, C * I * I / (8.0 * std::numbers::pi * a * a) };
}

} // namespace vktest
