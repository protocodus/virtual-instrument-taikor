#include "ModalExterior.h"

#include <algorithm>
#include <cmath>

namespace taikor::candidates::modal
{
namespace
{
constexpr int panelNodes = 8;       // Gauss nodes per panel of width pi/a
constexpr int nearNodes = 32;       // Gauss nodes in u on [k, kappa_p]

struct Rule
{
    std::vector<double> x, w;
};

const Rule& gauss (int n)
{
    // Rules are built once per size and reused (thread-safe static init).
    static const Rule r8 = [] { Rule r; special::gaussLegendre (panelNodes, r.x, r.w); return r; }();
    static const Rule r32 = [] { Rule r; special::gaussLegendre (nearNodes, r.x, r.w); return r; }();
    if (n == panelNodes)
        return r8;
    return r32;
}

void jAndPrime (int m, double x, double& j, double& jp)
{
    thread_local std::vector<double> seq;
    special::besselJSequence (x, m + 1, seq);
    j = seq[static_cast<std::size_t> (m)];
    jp = 0.5 * ((m == 0 ? -seq[1] : seq[static_cast<std::size_t> (m - 1)]) - seq[static_cast<std::size_t> (m + 1)]);
}

int resistanceCount (double k, double a) { return 24 + static_cast<int> (std::ceil (1.5 * k * a)); }

// First panel boundary at or beyond k + 2 pi / a.
int firstPanel (double k, double a) { return static_cast<int> (std::ceil ((k + 2.0 * pi / a) / (pi / a))); }
} // namespace

ExteriorTable::ExteriorTable (double radius, double kappaMaxIn, int maxOrder)
    : a (radius), kappaMax (kappaMaxIn), orders (maxOrder)
{
    const double h = pi / a;
    const int panels = static_cast<int> (std::ceil (kappaMax / h));
    kappaMax = panels * h;
    const Rule& g = gauss (panelNodes);
    for (int p = 0; p < panels; ++p)
        for (int i = 0; i < panelNodes; ++i)
        {
            panelKappa.push_back ((p + 0.5 * (g.x[static_cast<std::size_t> (i)] + 1.0)) * h);
            panelWeight.push_back (0.5 * h * g.w[static_cast<std::size_t> (i)]);
        }
    const std::size_t stride = static_cast<std::size_t> (orders + 1);
    panelJ.resize (panelKappa.size() * stride);
    panelJp.resize (panelKappa.size() * stride);
    thread_local std::vector<double> seq;
    for (std::size_t n = 0; n < panelKappa.size(); ++n)
    {
        special::besselJSequence (panelKappa[n] * a, orders + 1, seq);
        for (int m = 0; m <= orders; ++m)
        {
            const auto k = static_cast<std::size_t> (m);
            panelJ[n * stride + k] = seq[k];
            panelJp[n * stride + k] = 0.5 * ((m == 0 ? -seq[1] : seq[k - 1]) - seq[k + 1]);
        }
    }
}

void ExteriorTable::reactanceNodes (int m, double k, std::vector<Node>& nodes) const
{
    nodes.clear();
    const int p0 = firstPanel (k, a);
    const double kappaP = p0 * pi / a;
    const Rule& g = gauss (nearNodes);
    if (k > 0.0)
    {
        const double uMax = std::acosh (kappaP / k);
        for (int i = 0; i < nearNodes; ++i)
        {
            const double u = 0.5 * uMax * (g.x[static_cast<std::size_t> (i)] + 1.0);
            const double kappa = k * std::cosh (u);
            Node node { kappa, 0.5 * uMax * g.w[static_cast<std::size_t> (i)] * kappa, 0.0, 0.0 };
            jAndPrime (m, kappa * a, node.j, node.jp);
            nodes.push_back (node);
        }
    }
    else
    {
        for (int i = 0; i < nearNodes; ++i)
        {
            const double kappa = 0.5 * kappaP * (g.x[static_cast<std::size_t> (i)] + 1.0);
            Node node { kappa, 0.5 * kappaP * g.w[static_cast<std::size_t> (i)], 0.0, 0.0 };
            jAndPrime (m, kappa * a, node.j, node.jp);
            nodes.push_back (node);
        }
    }
    const std::size_t stride = static_cast<std::size_t> (orders + 1);
    for (std::size_t n = static_cast<std::size_t> (p0) * panelNodes; n < panelKappa.size(); ++n)
    {
        const double kappa = panelKappa[n];
        nodes.push_back ({ kappa, panelWeight[n] * kappa / std::sqrt (kappa * kappa - k * k),
                           panelJ[n * stride + static_cast<std::size_t> (m)], panelJp[n * stride + static_cast<std::size_t> (m)] });
    }
}

void ExteriorTable::resistanceNodes (int m, double k, std::vector<Node>& nodes) const
{
    nodes.clear();
    const int count = resistanceCount (k, a);
    std::vector<double> x, w;
    special::gaussLegendre (count, x, w);
    for (int i = 0; i < count; ++i)
    {
        const double t = 0.25 * pi * (x[static_cast<std::size_t> (i)] + 1.0);
        const double kappa = k * std::sin (t);
        Node node { kappa, 0.25 * pi * w[static_cast<std::size_t> (i)] * kappa, 0.0, 0.0 };
        jAndPrime (m, kappa * a, node.j, node.jp);
        nodes.push_back (node);
    }
}

namespace
{
void gram (const std::vector<const HeadMode*>& modes, double radius, const std::vector<double>& kappa,
           const std::vector<double>& weight, const std::vector<double>& j, const std::vector<double>& jp,
           double scale, std::vector<double>& out)
{
    const std::size_t n = modes.size(), nodes = kappa.size();
    std::vector<double> w (n * nodes);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t p = 0; p < nodes; ++p)
            w[i * nodes + p] = hankelTransform (*modes[i], radius, kappa[p], j[p], jp[p]);
    out.assign (n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = i; k < n; ++k)
        {
            double sum = 0.0;
            const double* wi = &w[i * nodes];
            const double* wk = &w[k * nodes];
            for (std::size_t p = 0; p < nodes; ++p)
                sum += weight[p] * wi[p] * wk[p];
            out[i * n + k] = out[k * n + i] = scale * sum;
        }
}
} // namespace

void ExteriorTable::addedMass (const std::vector<const HeadMode*>& modes, double k, double density,
                               std::vector<double>& out) const
{
    if (modes.empty())
    {
        out.clear();
        return;
    }
    const int m = modes.front()->m;
    std::vector<Node> nodes;
    reactanceNodes (m, k, nodes);
    std::vector<double> kap, wt, j, jp;
    for (const auto& node : nodes)
    {
        kap.push_back (node.kappa);
        wt.push_back (node.weight);
        j.push_back (node.j);
        jp.push_back (node.jp);
    }
    gram (modes, a, kap, wt, j, jp, density * angularNorm (m), out);
}

void ExteriorTable::resistance (const std::vector<const HeadMode*>& modes, double k, double density,
                                double soundSpeed, std::vector<double>& out) const
{
    if (modes.empty())
    {
        out.clear();
        return;
    }
    const int m = modes.front()->m;
    std::vector<Node> nodes;
    resistanceNodes (m, k, nodes);
    std::vector<double> kap, wt, j, jp;
    for (const auto& node : nodes)
    {
        kap.push_back (node.kappa);
        wt.push_back (node.weight);
        j.push_back (node.j);
        jp.push_back (node.jp);
    }
    gram (modes, a, kap, wt, j, jp, density * soundSpeed * k * angularNorm (m), out);
}

double ExteriorTable::selfAddedMass (const HeadMode& mode, double k, double density) const
{
    std::vector<double> out;
    addedMass ({ &mode }, k, density, out);
    return out[0];
}

double ExteriorTable::selfResistance (const HeadMode& mode, double k, double density, double soundSpeed) const
{
    std::vector<double> out;
    resistance ({ &mode }, k, density, soundSpeed, out);
    return out[0];
}

void kingIntegrals (int m, double radius, double k, double kappaMax, double density, double soundSpeed,
                    const std::function<double (double)>& transform, double& addedMassOut, double& resistanceOut)
{
    const double a = radius;
    const Rule& near = gauss (nearNodes);
    const Rule& panel = gauss (panelNodes);
    const int p0 = firstPanel (k, a);
    const double kappaP = p0 * pi / a;
    double sum = 0.0;
    const double uMax = std::acosh (kappaP / k);
    for (int i = 0; i < nearNodes; ++i)
    {
        const double u = 0.5 * uMax * (near.x[static_cast<std::size_t> (i)] + 1.0);
        const double kappa = k * std::cosh (u);
        const double f = transform (kappa);
        sum += 0.5 * uMax * near.w[static_cast<std::size_t> (i)] * kappa * f * f;
    }
    const double h = pi / a;
    for (double left = kappaP; left < kappaMax; left += h)
        for (int i = 0; i < panelNodes; ++i)
        {
            const double kappa = left + 0.5 * h * (panel.x[static_cast<std::size_t> (i)] + 1.0);
            const double f = transform (kappa);
            sum += 0.5 * h * panel.w[static_cast<std::size_t> (i)] * kappa / std::sqrt (kappa * kappa - k * k) * f * f;
        }
    addedMassOut = density * angularNorm (m) * sum;
    const int count = resistanceCount (k, a);
    std::vector<double> x, w;
    special::gaussLegendre (count, x, w);
    double res = 0.0;
    for (int i = 0; i < count; ++i)
    {
        const double t = 0.25 * pi * (x[static_cast<std::size_t> (i)] + 1.0);
        const double kappa = k * std::sin (t);
        const double f = transform (kappa);
        res += 0.25 * pi * w[static_cast<std::size_t> (i)] * kappa * f * f;
    }
    resistanceOut = density * soundSpeed * k * angularNorm (m) * res;
}

} // namespace taikor::candidates::modal
