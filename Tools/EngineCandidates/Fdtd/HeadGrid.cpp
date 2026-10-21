#include "HeadGrid.h"

#include <algorithm>
#include <cmath>

namespace taikor::candidates::fdtd
{

double stableHeadSpacing (const HeadSpec& head, double k, double extraTension) noexcept
{
    const double c2 = (head.tension + extraTension) / head.arealDensity;
    const double kappa2 = head.bendingStiffness / head.arealDensity;
    const double a = c2 * k * k + 4.0 * head.sigma1 * k;
    return std::sqrt (a + std::sqrt (a * a + 16.0 * kappa2 * k * k));
}

void HeadGrid::build (const HeadSpec& spec, double radius, double timeStep, double spacing,
                      double inclusionRadius, bool lossless)
{
    head = spec;
    headRadius = radius;
    k = timeStep;
    h = spacing;
    s0 = lossless ? 0.0 : spec.sigma0;
    s1 = lossless ? 0.0 : spec.sigma1;
    mesh = makeDiskGrid (inclusionRadius, h);
    const int nx = mesh.nx, ny = mesh.ny;
    interiorBegin.assign (static_cast<std::size_t> (ny), 0);
    interiorEnd.assign (static_cast<std::size_t> (ny), 0);
    interior = 0;
    for (int j = 0; j < ny; ++j)
    {
        int first = -1, last = -1;
        for (int i = 0; i < nx; ++i)
            if (isInterior (mesh.index (i, j)))
            {
                if (first < 0)
                    first = i;
                last = i;
            }
        if (first >= 0)
        {
            // The staircase disk is convex along rows: one contiguous run.
            interiorBegin[static_cast<std::size_t> (j)] = first;
            interiorEnd[static_cast<std::size_t> (j)] = last + 1;
            interior += last + 1 - first;
        }
    }
    // Interior plus first ring: the row's own run widened by one, joined with
    // the runs of the rows above and below (each of their interior nodes has
    // this row's node as a vertical neighbour).
    extendedBegin.assign (static_cast<std::size_t> (ny), 0);
    extendedEnd.assign (static_cast<std::size_t> (ny), 0);
    for (int j = 1; j < ny - 1; ++j)
    {
        int b = nx, e = 0;
        const auto join = [&] (int row, int widen)
        {
            const auto r = static_cast<std::size_t> (row);
            if (interiorEnd[r] > interiorBegin[r])
            {
                b = std::min (b, interiorBegin[r] - widen);
                e = std::max (e, interiorEnd[r] + widen);
            }
        };
        join (j, 1);
        join (j - 1, 0);
        join (j + 1, 0);
        if (e > b)
        {
            extendedBegin[static_cast<std::size_t> (j)] = b;
            extendedEnd[static_cast<std::size_t> (j)] = e;
        }
    }
}

void HeadGrid::laplacian (const double* w, double* lap, int j0, int j1) const noexcept
{
    const int nx = mesh.nx;
    const double scale = 1.0 / (h * h);
    j0 = std::max (j0, 1);
    j1 = std::min (j1, mesh.ny - 1);
    for (int j = j0; j < j1; ++j)
    {
        const int b = extBegin (j), e = extEnd (j);
        const double* row = w + j * nx;
        double* out = lap + j * nx;
        for (int i = b; i < e; ++i)
            out[i] = scale * (row[i - 1] + row[i + 1] + row[i - nx] + row[i + nx] - 4.0 * row[i]);
    }
}

void HeadGrid::explicitPart (const double* w, const double* wPrev, const double* lap,
                             const double* lapPrev, double* out, int j0, int j1) const noexcept
{
    const int nx = mesh.nx;
    const double invH2 = 1.0 / (h * h);
    const double kk = k * k / head.arealDensity;
    const double cT = kk * head.tension;
    const double cD = kk * head.bendingStiffness * invH2;
    const double cPrev = 1.0 - s0 * k;
    const double cLoss = 2.0 * s1 * k;
    for (int j = j0; j < j1; ++j)
    {
        const int b = rowBegin (j), e = rowEnd (j);
        const int base = j * nx;
        for (int i = b; i < e; ++i)
        {
            const int l = base + i;
            const double bilap = lap[l - 1] + lap[l + 1] + lap[l - nx] + lap[l + nx] - 4.0 * lap[l];
            out[l] = 2.0 * w[l] - cPrev * wPrev[l] + cT * lap[l] - cD * bilap
                     + cLoss * (lap[l] - lapPrev[l]);
        }
    }
}

double HeadGrid::energy (const double* next, const double* now) const noexcept
{
    const int nx = mesh.nx, ny = mesh.ny;
    const double sigma = head.arealDensity;
    double kinetic = 0.0, massCorrection = 0.0, tension = 0.0, bending = 0.0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            const int l = mesh.index (i, j);
            if (isInterior (l))
            {
                const double d = next[l] - now[l];
                kinetic += d * d;
            }
            // Edges to the right and up; w = 0 on clamped nodes, so edges
            // with one clamped end are counted and clamped-clamped ones vanish.
            if (i + 1 < nx)
            {
                const double d1 = next[l + 1] - next[l], d0 = now[l + 1] - now[l];
                tension += d1 * d0;
                const double dd = d1 - d0;
                massCorrection += dd * dd;
            }
            if (j + 1 < ny)
            {
                const double d1 = next[l + nx] - next[l], d0 = now[l + nx] - now[l];
                tension += d1 * d0;
                const double dd = d1 - d0;
                massCorrection += dd * dd;
            }
        }
    // Bending: (D h^2 / 2) sum over interior + first ring of (L w^{n+1})(L w^n).
    std::vector<double> lapNext (static_cast<std::size_t> (mesh.size()), 0.0), lapNow (lapNext);
    laplacian (next, lapNext.data(), 0, ny);
    laplacian (now, lapNow.data(), 0, ny);
    for (int j = 1; j < ny - 1; ++j)
        for (int i = extBegin (j); i < extEnd (j); ++i)
        {
            const int l = mesh.index (i, j);
            bending += lapNext[static_cast<std::size_t> (l)] * lapNow[static_cast<std::size_t> (l)];
        }
    const double h2 = h * h;
    return 0.5 * sigma * h2 / (k * k) * (kinetic - s1 * k / h2 * massCorrection)
           + 0.5 * head.tension * tension + 0.5 * head.bendingStiffness * h2 * bending;
}

double HeadGrid::lossPower (const double* next, const double* prev) const noexcept
{
    const int nx = mesh.nx, ny = mesh.ny;
    const double sigma = head.arealDensity;
    const double inv2k = 1.0 / (2.0 * k);
    double sum0 = 0.0, sum1 = 0.0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
            const int l = mesh.index (i, j);
            const double v = (next[l] - prev[l]) * inv2k;
            if (isInterior (l))
                sum0 += v * v;
            if (i + 1 < nx)
            {
                const double d = ((next[l + 1] - prev[l + 1]) - (next[l] - prev[l])) * inv2k;
                sum1 += d * d;
            }
            if (j + 1 < ny)
            {
                const double d = ((next[l + nx] - prev[l + nx]) - (next[l] - prev[l])) * inv2k;
                sum1 += d * d;
            }
        }
    return 2.0 * sigma * s0 * h * h * sum0 + 2.0 * sigma * s1 * sum1;
}

double HeadGrid::dot (const double* a, const double* b) const noexcept
{
    double sum = 0.0;
    for (int j = 0; j < mesh.ny; ++j)
    {
        const int base = j * mesh.nx;
        for (int i = rowBegin (j); i < rowEnd (j); ++i)
            sum += a[base + i] * b[base + i];
    }
    return sum;
}

std::vector<std::vector<double>> simulateInVacuo (const HeadGrid& head, const std::vector<double>& shape,
                                                  const std::vector<std::vector<double>>& probes, int steps)
{
    const auto n = static_cast<std::size_t> (head.size());
    std::vector<double> prev (shape), now (shape), next (n, 0.0);
    std::vector<double> lapNow (n, 0.0), lapPrev (n, 0.0);
    for (std::size_t l = 0; l < n; ++l)
        if (! head.isInterior (static_cast<int> (l)))
            prev[l] = now[l] = 0.0;
    head.laplacian (prev.data(), lapPrev.data(), 0, head.rows());
    std::vector<std::vector<double>> traces (probes.size());
    const double invAlpha = 1.0 / head.alpha();
    for (int step = 0; step < steps; ++step)
    {
        for (std::size_t p = 0; p < probes.size(); ++p)
            traces[p].push_back (head.dot (now.data(), probes[p].data()));
        head.laplacian (now.data(), lapNow.data(), 0, head.rows());
        head.explicitPart (now.data(), prev.data(), lapNow.data(), lapPrev.data(), next.data(), 0, head.rows());
        for (int j = 0; j < head.rows(); ++j)
        {
            const int base = j * head.grid().nx;
            for (int i = head.rowBegin (j); i < head.rowEnd (j); ++i)
                next[static_cast<std::size_t> (base + i)] *= invAlpha;
        }
        std::swap (prev, now);
        std::swap (now, next);
        std::swap (lapPrev, lapNow);
    }
    return traces;
}

double pronyFrequency (const std::vector<double>& s, double k) noexcept
{
    double num = 0.0, den = 0.0;
    for (std::size_t n = 1; n + 1 < s.size(); ++n)
    {
        num += s[n] * (s[n + 1] + s[n - 1]);
        den += 2.0 * s[n] * s[n];
    }
    const double c = std::clamp (num / den, -1.0, 1.0);
    return std::acos (c) / (2.0 * 3.14159265358979323846 * k);
}

} // namespace taikor::candidates::fdtd
