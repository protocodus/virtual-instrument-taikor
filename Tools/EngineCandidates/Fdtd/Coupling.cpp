#include "Coupling.h"

#include <algorithm>
#include <cmath>
#include <map>

namespace taikor::candidates::fdtd
{

SparseMatrix sparseProduct (const SparseMatrix& a, const SparseMatrix& b)
{
    SparseMatrix result;
    result.rows = a.rows;
    result.cols = b.cols;
    result.start.assign (static_cast<std::size_t> (a.rows + 1), 0);
    std::vector<double> accumulator (static_cast<std::size_t> (b.cols), 0.0);
    std::vector<int> touched;
    std::vector<char> seen (static_cast<std::size_t> (b.cols), 0);
    for (int r = 0; r < a.rows; ++r)
    {
        touched.clear();
        for (int e = a.start[static_cast<std::size_t> (r)]; e < a.start[static_cast<std::size_t> (r) + 1]; ++e)
        {
            const int mid = a.column[static_cast<std::size_t> (e)];
            const double av = a.value[static_cast<std::size_t> (e)];
            for (int f = b.start[static_cast<std::size_t> (mid)]; f < b.start[static_cast<std::size_t> (mid) + 1]; ++f)
            {
                const int col = b.column[static_cast<std::size_t> (f)];
                if (! seen[static_cast<std::size_t> (col)])
                {
                    seen[static_cast<std::size_t> (col)] = 1;
                    touched.push_back (col);
                }
                accumulator[static_cast<std::size_t> (col)] += av * b.value[static_cast<std::size_t> (f)];
            }
        }
        std::sort (touched.begin(), touched.end());
        for (const int col : touched)
        {
            result.column.push_back (col);
            result.value.push_back (accumulator[static_cast<std::size_t> (col)]);
            accumulator[static_cast<std::size_t> (col)] = 0.0;
            seen[static_cast<std::size_t> (col)] = 0;
        }
        result.start[static_cast<std::size_t> (r) + 1] = static_cast<int> (result.value.size());
    }
    return result;
}

SparseMatrix sparseTranspose (const SparseMatrix& a, double scale)
{
    SparseMatrix t;
    t.rows = a.cols;
    t.cols = a.rows;
    t.start.assign (static_cast<std::size_t> (t.rows + 1), 0);
    for (const int col : a.column)
        ++t.start[static_cast<std::size_t> (col) + 1];
    for (int r = 0; r < t.rows; ++r)
        t.start[static_cast<std::size_t> (r) + 1] += t.start[static_cast<std::size_t> (r)];
    t.column.assign (a.column.size(), 0);
    t.value.assign (a.value.size(), 0.0);
    std::vector<int> fill (t.start.begin(), t.start.end() - 1);
    for (int r = 0; r < a.rows; ++r)
        for (int e = a.start[static_cast<std::size_t> (r)]; e < a.start[static_cast<std::size_t> (r) + 1]; ++e)
        {
            const int col = a.column[static_cast<std::size_t> (e)];
            const auto slot = static_cast<std::size_t> (fill[static_cast<std::size_t> (col)]++);
            t.column[slot] = r;
            t.value[slot] = scale * a.value[static_cast<std::size_t> (e)];
        }
    return t;
}

void HeadCoupling::build (const HeadGrid& head, const AirPlaneLayout& layout,
                          const std::vector<std::array<int, 2>>& columns, double nuValue, int neumannOrder)
{
    columnList = columns;
    nu = nuValue;
    const auto& g = head.grid();
    std::vector<int> lookup (static_cast<std::size_t> (layout.nx * layout.ny), -1);
    for (std::size_t c = 0; c < columns.size(); ++c)
        lookup[static_cast<std::size_t> (columns[c][1] * layout.nx + columns[c][0])] = static_cast<int> (c);
    // I: quadratic B-spline weights of the 3 x 3 cell centres nearest the
    // node (w(-1) = (1/2 - t)^2 / 2, w(0) = 3/4 - t^2, w(+1) = (1/2 + t)^2 / 2
    // per axis, t the offset from the nearest centre in cells), restricted to
    // columns inside the head and renormalised. Its transpose J averages the
    // head's velocity with a kernel whose spectrum falls as sinc^3 rather
    // than the bilinear hat's sinc^2: head bending waves shorter than two air
    // cells (subsonic, 10-16 kHz on the nagado batter) would otherwise alias
    // on the air grid into supersonic, radiating components. At 1.5 pi/Delta
    // the kernel passes -31 dB where the bilinear one passes -21 dB.
    interpolate = SparseMatrix(); // build() may run again when the heads are rebuilt
    interpolate.rows = g.size();
    interpolate.cols = planeSize();
    interpolate.start.assign (static_cast<std::size_t> (g.size() + 1), 0);
    for (int l = 0; l < g.size(); ++l)
    {
        if (head.isInterior (l))
        {
            const double fx = (g.x (l % g.nx) - layout.x0) / layout.delta - 0.5;
            const double fy = (g.y (l / g.nx) - layout.y0) / layout.delta - 0.5;
            const int ic = static_cast<int> (std::lround (fx)), jc = static_cast<int> (std::lround (fy));
            const double tx = fx - ic, ty = fy - jc;
            const double wx[3] = { 0.5 * (0.5 - tx) * (0.5 - tx), 0.75 - tx * tx, 0.5 * (0.5 + tx) * (0.5 + tx) };
            const double wy[3] = { 0.5 * (0.5 - ty) * (0.5 - ty), 0.75 - ty * ty, 0.5 * (0.5 + ty) * (0.5 + ty) };
            int cols[9];
            double weights[9];
            int count = 0;
            double total = 0.0;
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx)
                {
                    const int i = ic + dx, j = jc + dy;
                    if (i < 0 || j < 0 || i >= layout.nx || j >= layout.ny)
                        continue;
                    const int c = lookup[static_cast<std::size_t> (j * layout.nx + i)];
                    const double weight = wx[dx + 1] * wy[dy + 1];
                    if (c < 0 || weight <= 0.0)
                        continue;
                    cols[count] = c;
                    weights[count++] = weight;
                    total += weight;
                }
            for (int e = 0; e < count; ++e)
            {
                interpolate.column.push_back (cols[e]);
                interpolate.value.push_back (weights[e] / total);
            }
        }
        interpolate.start[static_cast<std::size_t> (l) + 1] = static_cast<int> (interpolate.value.size());
    }
    const double h = head.spacing();
    transposeScale = layout.delta * layout.delta / (h * h);
    spread = sparseTranspose (interpolate, 1.0 / transposeScale);
    gram = sparseProduct (spread, interpolate);
    buildInverse (nullptr, neumannOrder);
}

void HeadCoupling::setSliceScale (const std::vector<double>& scale, int neumannOrder)
{
    buildInverse (&scale, neumannOrder);
}

void HeadCoupling::buildInverse (const std::vector<double>* scale, int neumannOrder)
{
    // Neumann series: sum_{p = 0..order} (-nu G S)^p (S = I without a scale).
    SparseMatrix identity;
    identity.rows = identity.cols = planeSize();
    for (int c = 0; c <= planeSize(); ++c)
        identity.start.push_back (c);
    identity.column.resize (static_cast<std::size_t> (planeSize()));
    for (int c = 0; c < planeSize(); ++c)
        identity.column[static_cast<std::size_t> (c)] = c;
    identity.value.assign (static_cast<std::size_t> (planeSize()), 1.0);
    SparseMatrix minusNuG = gram;
    for (auto& v : minusNuG.value)
        v *= -nu;
    if (scale)
        for (std::size_t e = 0; e < minusNuG.value.size(); ++e)
            minusNuG.value[e] *= (*scale)[static_cast<std::size_t> (minusNuG.column[e])];
    SparseMatrix term = identity;
    std::map<std::pair<int, int>, double> sum;
    for (int p = 0; p <= neumannOrder; ++p)
    {
        for (int r = 0; r < term.rows; ++r)
            for (int e = term.start[static_cast<std::size_t> (r)]; e < term.start[static_cast<std::size_t> (r) + 1]; ++e)
                sum[{ r, term.column[static_cast<std::size_t> (e)] }] += term.value[static_cast<std::size_t> (e)];
        if (p < neumannOrder)
            term = sparseProduct (term, minusNuG);
    }
    inverse = SparseMatrix();
    inverse.rows = inverse.cols = planeSize();
    inverse.start.assign (static_cast<std::size_t> (planeSize() + 1), 0);
    for (const auto& [key, value] : sum)
    {
        inverse.column.push_back (key.second);
        inverse.value.push_back (value);
        ++inverse.start[static_cast<std::size_t> (key.first) + 1];
    }
    for (int r = 0; r < planeSize(); ++r)
        inverse.start[static_cast<std::size_t> (r) + 1] += inverse.start[static_cast<std::size_t> (r)];
}

} // namespace taikor::candidates::fdtd
