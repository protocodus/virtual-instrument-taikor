#include "AirDomain.h"

#include <algorithm>
#include <cmath>

namespace taikor::candidates::fdtd
{

template <typename Real>
void AirDomain<Real>::build (const AirGeometry& g, double soundSpeed, double rho, double timeStep)
{
    c = soundSpeed;
    density = rho;
    k = timeStep;
    delta = g.delta;
    lambda = c * k / delta;
    pml = g.pmlCells;
    const int border = pml + 1; // PML cells plus the fixed (or rigid) outer layer
    const int nh = static_cast<int> (std::ceil ((g.radius + g.wall + g.margin) / delta - 1.0e-9));
    const int nDepth = static_cast<int> (std::lround (g.depth / delta));
    const int nBelow = static_cast<int> (std::ceil (g.margin / delta - 1.0e-9));
    const int nAbove = static_cast<int> (std::ceil ((g.margin + g.marginAbove) / delta - 1.0e-9));
    nx = ny = 2 * (nh + border);
    nz = nDepth + nBelow + nAbove + 2 * border;
    x0 = y0 = -(nh + border) * delta;
    z0 = -(nDepth + nBelow + border) * delta;
    const auto total = static_cast<std::size_t> (cells());
    {
        const int fields = g.pmlCells > 0 ? 11 : 3;
        // Field m starts at byte m * 1408 (11 x 128 B) modulo 16 KB, the
        // set period of the M1's 8-way 128 KB L1 data cache, so the fields'
        // same-index lines spread over the sets.
        constexpr std::size_t period = 16384 / sizeof (Real), step = 1408 / sizeof (Real);
        std::size_t offset = 0;
        std::vector<std::size_t> starts;
        for (int m = 0; m < fields; ++m)
        {
            const std::size_t want = (static_cast<std::size_t> (m) * step) % period;
            offset += (want + period - offset % period) % period;
            starts.push_back (offset);
            offset += total;
        }
        storage.assign (offset + 64, Real (0));
        for (int m = 0; m < 3; ++m)
            psiBase[m] = storage.data() + starts[static_cast<std::size_t> (m)];
        if (fields == 11)
        {
            auxBase[0] = storage.data() + starts[3];
            auxBase[1] = storage.data() + starts[4];
            for (int p = 0; p < 2; ++p)
                for (int d = 0; d < 3; ++d)
                    phiBase[p][d] = storage.data() + starts[static_cast<std::size_t> (5 + 3 * p + d)];
        }
    }
    types.assign (total, CellType::Regular);
    masks.assign (total, 0);
    const double a = g.radius, outer = g.radius + g.wall, depth = nDepth * delta;
    for (int kz = 0; kz < nz; ++kz)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const auto cell = static_cast<std::size_t> (index (i, j, kz));
                const bool edge = i == 0 || j == 0 || kz == 0 || i == nx - 1 || j == ny - 1 || kz == nz - 1;
                const bool layer = i < border || j < border || kz < border
                                   || i >= nx - border || j >= ny - border || kz >= nz - border;
                const double r = std::hypot (x (i), y (j)), zc = z (kz);
                if (edge)
                    types[cell] = pml > 0 ? CellType::Fixed : CellType::Solid;
                else if (pml > 0 && layer)
                    types[cell] = CellType::Pml;
                else if (zc > -depth && zc < 0.0 && r >= a && r < outer)
                    types[cell] = CellType::Solid;
            }
    // Open faces: both cells air (anything but Solid and Fixed).
    const auto air = [this] (std::size_t cell)
    { return types[cell] != CellType::Solid && types[cell] != CellType::Fixed; };
    const int stride[3] = { 1, nx, nx * ny };
    for (int kz = 1; kz < nz - 1; ++kz)
        for (int j = 1; j < ny - 1; ++j)
            for (int i = 1; i < nx - 1; ++i)
            {
                const auto cell = static_cast<std::size_t> (index (i, j, kz));
                if (! air (cell))
                    continue;
                std::uint8_t mask = 0;
                for (int d = 0; d < 3; ++d)
                {
                    if (air (cell - static_cast<std::size_t> (stride[d])))
                        mask |= static_cast<std::uint8_t> (1u << (2 * d));
                    if (air (cell + static_cast<std::size_t> (stride[d])))
                        mask |= static_cast<std::uint8_t> (1u << (2 * d + 1));
                }
                masks[cell] = mask;
            }
    // PML profiles z(d) = zmax (d / thickness)^3 with zmax = 2 c ln(1/R) / thickness,
    // so the normal-incidence round trip through the layer is exp(-zmax thickness / 2c) = R.
    for (int d = 0; d < 3; ++d)
    {
        const int n = d == 0 ? nx : (d == 1 ? ny : nz);
        zetaCell[d].assign (static_cast<std::size_t> (n), 0.0);
        zetaFace[d].assign (static_cast<std::size_t> (n), 0.0);
        if (pml == 0)
            continue;
        const double thickness = pml * delta;
        const double zmax = 2.0 * c * std::log (1.0 / g.pmlReflection) / thickness;
        const auto profile = [&] (double depthIn)
        {
            const double s = std::clamp (depthIn / thickness, 0.0, 1.0);
            return zmax * s * s * s;
        };
        for (int i = 0; i < n; ++i)
        {
            // Low side: interface face at index border (between cells pml and pml+1).
            const double lowCell = (border - (i + 0.5)) * delta, lowFace = (border - (i + 1.0)) * delta;
            const double highCell = ((i + 0.5) - (n - border)) * delta, highFace = ((i + 1.0) - (n - border)) * delta;
            zetaCell[d][static_cast<std::size_t> (i)] = profile (std::max (lowCell, highCell));
            zetaFace[d][static_cast<std::size_t> (i)] = profile (std::max (lowFace, highFace));
        }
    }
    for (int d = 0; d < 3; ++d)
    {
        zetaCellR[d].assign (zetaCell[d].begin(), zetaCell[d].end());
        zetaFaceR[d].assign (zetaFace[d].begin(), zetaFace[d].end());
        faceDecay[d].assign (zetaFace[d].size(), Real (1));
        faceGain[d].assign (zetaFace[d].size(), Real (0));
        for (std::size_t i = 0; i < zetaFace[d].size(); ++i)
        {
            const double half = 0.5 * k * zetaFace[d][i];
            faceDecay[d][i] = static_cast<Real> ((1.0 - half) / (1.0 + half));
            faceGain[d][i] = static_cast<Real> (k * c * c / (delta * (1.0 + half)));
        }
        // The face between the fixed outer cell 0 and the first layer cell
        // carries no auxiliary field (phi = 0 there, as on the inner face).
        faceDecay[d][0] = faceGain[d][0] = Real (0);
    }
    finalise();
}

template <typename Real>
void AirDomain<Real>::closeHeadFaces (int kzAbove, const std::vector<std::array<int, 2>>& columns, bool markSlices)
{
    for (const auto& column : columns)
    {
        const auto above = static_cast<std::size_t> (index (column[0], column[1], kzAbove));
        const auto below = static_cast<std::size_t> (index (column[0], column[1], kzAbove - 1));
        masks[above] = static_cast<std::uint8_t> (masks[above] & ~FaceZm);
        masks[below] = static_cast<std::uint8_t> (masks[below] & ~FaceZp);
        // Only air cells become Slice cells: a drum without exterior air has
        // the rigid boundary layer on the head's outer side.
        const auto air = [this] (std::size_t cell)
        { return types[cell] != CellType::Solid && types[cell] != CellType::Fixed; };
        if (markSlices && air (above))
            types[above] = CellType::Slice;
        if (markSlices && air (below))
            types[below] = CellType::Slice;
    }
}

template <typename Real>
std::vector<std::array<int, 2>> AirDomain<Real>::headColumns (double radius) const
{
    std::vector<std::array<int, 2>> columns;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            if (std::hypot (x (i), y (j)) < radius)
                columns.push_back ({ i, j });
    return columns;
}

template <typename Real>
void AirDomain<Real>::finalise()
{
    runs.clear();
    irregular.clear();
    runsOfPlane.assign (static_cast<std::size_t> (nz + 1), 0);
    for (int kz = 0; kz < nz; ++kz)
    {
        runsOfPlane[static_cast<std::size_t> (kz)] = static_cast<int> (runs.size());
        for (int j = 0; j < ny; ++j)
        {
            int start = -1;
            for (int i = 0; i <= nx; ++i)
            {
                bool regular = false;
                if (i < nx)
                {
                    const auto cell = static_cast<std::size_t> (index (i, j, kz));
                    auto& t = types[cell];
                    if (t == CellType::Regular || t == CellType::Irregular)
                    {
                        t = masks[cell] == 63 ? CellType::Regular : CellType::Irregular;
                        if (t == CellType::Irregular)
                            irregular.push_back (static_cast<int> (cell));
                    }
                    regular = t == CellType::Regular;
                }
                if (regular && start < 0)
                    start = i;
                if (! regular && start >= 0)
                {
                    runs.push_back ({ index (start, j, kz), i - start });
                    start = -1;
                }
            }
        }
    }
    runsOfPlane[static_cast<std::size_t> (nz)] = static_cast<int> (runs.size());
}

template <typename Real>
void AirDomain<Real>::listInterior (double radius, double depth)
{
    interiorCells.clear();
    for (int kz = 0; kz < nz; ++kz)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                const int cell = index (i, j, kz);
                const CellType t = types[static_cast<std::size_t> (cell)];
                if ((t == CellType::Regular || t == CellType::Irregular || t == CellType::Slice)
                    && std::hypot (x (i), y (j)) < radius && z (kz) < 0.0 && z (kz) > -depth)
                    interiorCells.push_back (cell);
            }
}

template <typename Real>
void AirDomain<Real>::shiftInterior (long long n, double shift, int begin, int end) noexcept
{
    Real* a = level (n);
    Real* b = level (n - 1);
    const Real s = static_cast<Real> (shift);
    for (int c = begin; c < end; ++c)
    {
        const int cell = interiorCells[static_cast<std::size_t> (c)];
        a[cell] += s;
        b[cell] += s;
    }
}

template <typename Real>
double AirDomain<Real>::interiorSample (long long n) const noexcept
{
    // The mean over a sparse, fixed subset is enough: any constant will do.
    const Real* a = level (n);
    double sum = 0.0;
    int count = 0;
    for (std::size_t c = 0; c < interiorCells.size(); c += 97)
    {
        sum += a[interiorCells[c]];
        ++count;
    }
    return count > 0 ? sum / count : 0.0;
}

template <typename Real>
long long AirDomain<Real>::pmlCellCount() const
{
    long long count = 0;
    for (const auto t : types)
        count += t == CellType::Pml ? 1 : 0;
    return count;
}

template <typename Real>
std::vector<int> AirDomain<Real>::pmlPlaneSplit (int parts) const
{
    std::vector<long long> perPlane (static_cast<std::size_t> (nz), 0);
    long long total = 0;
    for (int kz = 0; kz < nz; ++kz)
    {
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                perPlane[static_cast<std::size_t> (kz)] += types[static_cast<std::size_t> (index (i, j, kz))] == CellType::Pml ? 1 : 0;
        total += perPlane[static_cast<std::size_t> (kz)];
    }
    std::vector<int> bounds { 0 };
    long long running = 0;
    for (int kz = 0; kz < nz; ++kz)
    {
        running += perPlane[static_cast<std::size_t> (kz)];
        if (static_cast<int> (bounds.size()) < parts && running * parts >= total * static_cast<long long> (bounds.size()))
            bounds.push_back (kz + 1);
    }
    while (static_cast<int> (bounds.size()) < parts + 1)
        bounds.push_back (nz);
    bounds.back() = nz;
    return bounds;
}

template <typename Real>
int AirDomain<Real>::planeAbove (double zFace) const noexcept
{
    return static_cast<int> (std::lround ((zFace - z0) / delta));
}

template <typename Real>
void AirDomain<Real>::stepRegular (long long n, int k0, int k1) noexcept
{
    const Real* __restrict cur = level (n);
    const Real* __restrict prev = level (n - 1);
    Real* __restrict next = level (n + 1);
    const Real l2 = static_cast<Real> (lambda * lambda);
    const Real c0 = static_cast<Real> (2.0 - 6.0 * lambda * lambda);
    const int sy = nx, sz = nx * ny;
    const int r0 = runsOfPlane[static_cast<std::size_t> (k0)], r1 = runsOfPlane[static_cast<std::size_t> (k1)];
    for (int r = r0; r < r1; ++r)
    {
        const int start = runs[static_cast<std::size_t> (r)].offset;
        const int end = start + runs[static_cast<std::size_t> (r)].length;
        for (int x = start; x < end; ++x)
            next[x] = c0 * cur[x] + l2 * (cur[x - 1] + cur[x + 1] + cur[x - sy] + cur[x + sy] + cur[x - sz] + cur[x + sz])
                      - prev[x];
    }
}

template <typename Real>
void AirDomain<Real>::stepMasked (long long n, int cell) noexcept
{
    const Real* cur = level (n);
    const Real* prev = level (n - 1);
    Real* next = level (n + 1);
    const int offsets[6] = { -1, 1, -nx, nx, -nx * ny, nx * ny };
    const std::uint8_t mask = masks[static_cast<std::size_t> (cell)];
    double sum = 0.0;
    int open = 0;
    for (int f = 0; f < 6; ++f)
        if (mask & (1u << f))
        {
            sum += static_cast<double> (cur[cell + offsets[f]]);
            ++open;
        }
    const double l2 = lambda * lambda;
    next[cell] = static_cast<Real> ((2.0 - open * l2) * static_cast<double> (cur[cell]) + l2 * sum
                                    - static_cast<double> (prev[cell]));
}

template <typename Real>
void AirDomain<Real>::stepIrregular (long long n, int begin, int end) noexcept
{
    for (int i = begin; i < end; ++i)
        stepMasked (n, irregular[static_cast<std::size_t> (i)]);
}

template class AirDomain<float>;
template class AirDomain<double>;

} // namespace taikor::candidates::fdtd
