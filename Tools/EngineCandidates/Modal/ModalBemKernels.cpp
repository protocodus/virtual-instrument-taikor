// Geometry and ring kernels of the axisymmetric exterior solver (ModalBem.h).

#include "ModalBem.h"
#include "ModalMath.h"

#include <cmath>
#include <map>

namespace taikor::candidates::modal
{
namespace
{
struct Rule
{
    std::vector<double> x, w; // on [-1, 1]
};

// Gauss-Legendre rules, built once per size and thread.
const Rule& gaussRule (int n)
{
    thread_local std::map<int, Rule> rules;
    auto& rule = rules[n];
    if (rule.x.empty())
        special::gaussLegendre (n, rule.x, rule.w);
    return rule;
}

// Radical inverse in base b (Halton sequence).
double radicalInverse (int i, int base)
{
    double f = 1.0, result = 0.0;
    for (int n = i; n > 0; n /= base)
    {
        f /= base;
        result += f * (n % base);
    }
    return result;
}
} // namespace

BemGeometry polylineGeometry (const std::vector<std::array<double, 2>>& nodes, const std::vector<int>& parts,
                              const std::vector<double>& steps, const std::vector<std::array<double, 2>>& chief)
{
    BemGeometry geo;
    geo.chief = chief;
    for (std::size_t s = 0; s + 1 < nodes.size(); ++s)
    {
        const double dr = nodes[s + 1][0] - nodes[s][0], dz = nodes[s + 1][1] - nodes[s][1];
        const double len = std::hypot (dr, dz);
        if (len <= 0.0)
            continue;
        const int count = std::max (1, static_cast<int> (std::ceil (len / steps[s] - 1e-9)));
        for (int p = 0; p < count; ++p)
        {
            BemPanel panel;
            panel.r0 = nodes[s][0] + dr * p / count;
            panel.z0 = nodes[s][1] + dz * p / count;
            panel.r1 = nodes[s][0] + dr * (p + 1) / count;
            panel.z1 = nodes[s][1] + dz * (p + 1) / count;
            panel.length = len / count;
            panel.tr = dr / len;
            panel.tz = dz / len;
            panel.nr = -panel.tz;
            panel.nz = panel.tr;
            panel.part = parts[s];
            geo.panels.push_back (panel);
        }
    }
    return geo;
}

BemGeometry drumGeometry (double a, double wall, double length, double headStep, double bodyStep, int chiefPoints)
{
    const double b = a + wall;
    const std::vector<std::array<double, 2>> nodes { { 0.0, 0.0 }, { a, 0.0 }, { b, 0.0 }, { b, -length },
                                                     { a, -length }, { 0.0, -length } };
    std::vector<std::array<double, 2>> chief;
    // Spread through the interior away from the walls and the axis (Halton, uniform in volume).
    for (int i = 1; i <= chiefPoints; ++i)
    {
        const double u = radicalInverse (i, 2), v = radicalInverse (i, 3);
        chief.push_back ({ b * std::sqrt (0.02 + 0.70 * u), -length * (0.12 + 0.76 * v) });
    }
    return polylineGeometry (nodes, { 0, 2, 2, 2, 1 }, { headStep, headStep, bodyStep, headStep, headStep }, chief);
}

std::array<double, bemVelocityNodes> bemNodeOffsets (const BemPanel& panel)
{
    const Rule& g = gaussRule (bemVelocityNodes);
    std::array<double, bemVelocityNodes> t {};
    for (int i = 0; i < bemVelocityNodes; ++i)
        t[static_cast<std::size_t> (i)] = 0.5 * panel.length * (g.x[static_cast<std::size_t> (i)] + 1.0);
    return t;
}

void ringKernels (double rx, double zx, double ry, double zy, double nr, double nz, double k, int maxOrder, int refine,
                  std::complex<double>* g, std::complex<double>* h)
{
    for (int m = 0; m <= maxOrder; ++m)
        g[m] = h[m] = 0.0;
    const double dz = zy - zx;
    const double d2 = (rx - ry) * (rx - ry) + dz * dz;
    const double prod = rx * ry;
    const double eps = prod > 0.0 ? std::sqrt (d2 / prod) : 1.0e300;
    // The static ring kernel of order m falls as e^{-m eta}, cosh eta = 1 + d^2 / (2 r_x r_y)
    // (Legendre Q_{m-1/2}); e^{-j k R} spreads it over about k sqrt(r_x r_y) more orders.
    // Orders beyond both are below 1e-17 of order 0 and left at zero.
    if (prod > 0.0 && maxOrder > 0)
    {
        const double eta = std::acosh (1.0 + 0.5 * eps * eps);
        const double reach = 40.0 / std::max (eta, 1.0e-12) + k * std::sqrt (prod) + 2.0;
        if (reach < maxOrder)
            maxOrder = static_cast<int> (std::ceil (reach));
    }
    // Gauss-Legendre on [0, pi] integrates cos(M phi) to 1e-13 with about M + 12 nodes.
    const int oscillation = static_cast<int> (std::ceil (maxOrder + k * std::sqrt (prod)));
    // Kernel values at every node first (weights folded in), then all orders by
    // the Chebyshev recurrence of cos(m phi) as dot products over the nodes.
    thread_local std::vector<double> cosine, gr, gi, hr, hi, tPrev, tCur;
    cosine.clear();
    gr.clear();
    gi.clear();
    hr.clear();
    hi.clear();
    const auto accumulate = [&] (double phi, double weight)
    {
        const double half = std::sin (0.5 * phi);
        const double oneMinusCos = 2.0 * half * half;
        const double c = 1.0 - oneMinusCos;
        const double r = std::sqrt (d2 + 2.0 * prod * oneMinusCos);
        const double cr = std::cos (k * r), sr = -std::sin (k * r);   // e^{-j k R}
        const double dRdn = ((ry - rx * c) * nr + dz * nz) / r;
        const double gs = weight / r;
        const double hs = -weight * dRdn / (r * r);                    // times (1 + j k R) e^{-j k R}
        cosine.push_back (c);
        gr.push_back (gs * cr);
        gi.push_back (gs * sr);
        hr.push_back (hs * (cr - k * r * sr));
        hi.push_back (hs * (sr + k * r * cr));
    };
    const auto gaussOn = [&] (double lo, double hi, int n)
    {
        const Rule& rule = gaussRule (n);
        for (int i = 0; i < n; ++i)
            accumulate (lo + 0.5 * (hi - lo) * (rule.x[static_cast<std::size_t> (i)] + 1.0),
                        0.5 * (hi - lo) * rule.w[static_cast<std::size_t> (i)]);
    };
    // phi1: the angle over which cos(m phi) e^{-j k R} turns by about 2 rad.
    const double phi1 = std::min (pi, 2.0 / (1.0 + maxOrder + k * std::sqrt (prod)));
    if (eps < phi1)
    {
        // Near: the 1/sqrt(eps^2 + phi^2) peak is narrower than the oscillation.
        // sinh map phi = eps sinh(u) on [0, phi1], regular Gauss beyond.
        const double uMax = std::asinh (phi1 / eps);
        const Rule& rule = gaussRule (refine * (8 + 2 * static_cast<int> (std::ceil (uMax))));
        for (std::size_t i = 0; i < rule.x.size(); ++i)
        {
            const double u = 0.5 * uMax * (rule.x[i] + 1.0);
            accumulate (eps * std::sinh (u), 0.5 * uMax * rule.w[i] * eps * std::cosh (u));
        }
        if (phi1 < pi)
            gaussOn (phi1, pi, refine * (8 + oscillation + static_cast<int> (std::ceil (3.0 * std::log (pi / phi1)))));
    }
    else
        gaussOn (0.0, pi, refine * (16 + oscillation + static_cast<int> (std::ceil (std::min (pi / eps, 1.0e4)))));
    const std::size_t count = cosine.size();
    tPrev.assign (count, 1.0);
    tCur = cosine;
    {
        double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
        for (std::size_t i = 0; i < count; ++i)
        {
            a0 += gr[i];
            a1 += gi[i];
            a2 += hr[i];
            a3 += hi[i];
        }
        g[0] = { a0, a1 };
        h[0] = { a2, a3 };
    }
    for (int m = 1; m <= maxOrder; ++m)
    {
        double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
        double* tp = tPrev.data();
        double* tc = tCur.data();
        const double* cs = cosine.data();
        for (std::size_t i = 0; i < count; ++i)
        {
            const double t = tc[i];
            a0 += t * gr[i];
            a1 += t * gi[i];
            a2 += t * hr[i];
            a3 += t * hi[i];
            const double next = 2.0 * cs[i] * t - tp[i];
            tp[i] = t;
            tc[i] = next;
        }
        g[m] = { a0, a1 };
        h[m] = { a2, a3 };
    }
    const double scale = 1.0 / (2.0 * pi);
    for (int m = 0; m <= maxOrder; ++m)
    {
        g[m] *= scale;
        h[m] *= scale;
    }
}

} // namespace taikor::candidates::modal
