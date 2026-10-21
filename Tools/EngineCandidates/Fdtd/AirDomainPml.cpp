// The absorbing layer, the microphones and the energy of the air domain
// (see AirDomain.h for the equations).

#include "AirDomain.h"

#include <algorithm>
#include <cmath>

namespace taikor::candidates::fdtd
{
namespace
{
// The PML cells of row (j, kz) as up to two i-ranges.
int pmlRanges (int nx, int ny, int nz, int border, int j, int kz, int ranges[4]) noexcept
{
    if (kz < border || kz >= nz - border || j < border || j >= ny - border)
    {
        ranges[0] = 1;
        ranges[1] = nx - 1;
        return 1;
    }
    ranges[0] = 1;
    ranges[1] = border;
    ranges[2] = nx - border;
    ranges[3] = nx - 1;
    return 2;
}
} // namespace

// Both kernels walk the layer row by row with the row's profile values
// hoisted, in the air's own precision, so the inner loops vectorise. Time
// staggering: u at integer steps, phi at half steps, psi (= int u dt) at
// integer steps, and the u update reads phi^{n+1/2} only (a leapfrog pair,
// as the velocity-pressure Yee scheme: phi^{n+1/2} is computed from u^n and
// psi^n after step n-1 completes).
namespace
{
// One contiguous run [i0, i1) of layer cells of row (j, kz); profiles are
// given as cell values and face values ((i + 1/2) at index i), with the face
// coefficients a = (1 - k z/2) / (1 + k z/2), b = k c^2 / (Delta (1 + k z/2)).
// Every cell
// forms phi^{n+1/2} on its six faces from the old buffer (phi^{n-1/2}) and
// the complete levels u^n, psi^n, stores its own three (+x, +y, +z) faces,
// and updates u and psi, so the pass needs no ordering between cells.
template <typename Real>
inline void pmlRow (const Real* __restrict cur, const Real* __restrict prev, Real* __restrict next,
                    const Real* __restrict psiNow, Real* __restrict psiNext,
                    const Real* __restrict oldX, const Real* __restrict oldY, const Real* __restrict oldZ,
                    Real* __restrict newX, Real* __restrict newY, Real* __restrict newZ,
                    const Real* __restrict zxc, const Real* __restrict zxf,
                    const Real* __restrict ax, const Real* __restrict bx,
                    int row, int i0, int i1, int sy, int sz,
                    Real zy, Real zz, Real zyP, Real zyM, Real zzP, Real zzM,
                    Real ayP, Real byP, Real ayM, Real byM, Real azP, Real bzP, Real azM, Real bzM,
                    Real cc, Real kk, Real kd, Real halfK) noexcept
{
    const Real yz = zy * zz;
    for (int i = i0; i < i1; ++i)
    {
        const int x = row + i;
        const Real zx = zxc[i];
        const Real u = cur[x], p0 = psiNow[x];
        // phi_j,t = -z_j phi_j + c^2 (z_i + z_k - z_j) d_j u + c^2 z_i z_k d_j psi, faces +j and -j.
        const Real fxP = ax[i] * oldX[x] + bx[i] * ((zy + zz - zxf[i]) * (cur[x + 1] - u) + yz * (psiNow[x + 1] - p0));
        const Real fxM = ax[i - 1] * oldX[x - 1] + bx[i - 1] * ((zy + zz - zxf[i - 1]) * (u - cur[x - 1]) + yz * (p0 - psiNow[x - 1]));
        const Real fyP = ayP * oldY[x] + byP * ((zz + zx - zyP) * (cur[x + sy] - u) + zz * zx * (psiNow[x + sy] - p0));
        const Real fyM = ayM * oldY[x - sy] + byM * ((zz + zx - zyM) * (u - cur[x - sy]) + zz * zx * (p0 - psiNow[x - sy]));
        const Real fzP = azP * oldZ[x] + bzP * ((zx + zy - zzP) * (cur[x + sz] - u) + zx * zy * (psiNow[x + sz] - p0));
        const Real fzM = azM * oldZ[x - sz] + bzM * ((zx + zy - zzM) * (u - cur[x - sz]) + zx * zy * (p0 - psiNow[x - sz]));
        newX[x] = fxP;
        newY[x] = fyP;
        newZ[x] = fzP;
        const Real s = zx + zy + zz, p = zx * (zy + zz) + yz, q = zx * yz;
        const Real lap = cur[x - 1] + cur[x + 1] + cur[x - sy] + cur[x + sy] + cur[x - sz] + cur[x + sz] - Real (6) * u;
        const Real div = (fxP - fxM) + (fyP - fyM) + (fzP - fzM);
        const Real value = (Real (2) * u - (Real (1) - halfK * s) * prev[x] + cc * lap + kd * div - kk * (p * u + q * p0))
                           / (Real (1) + halfK * s);
        next[x] = value;
        psiNext[x] = p0 + halfK * (value + u);
    }
}

// The same for cells where only the x-profile can be non-zero (zy = zz = 0 on
// the row, or zx = 0 with only one of zy, zz non-zero handled by the caller
// swapping nothing: it is used where exactly one of the three profiles is
// non-zero). There P = Q = 0 and every face's z_i z_k vanishes, so psi enters
// no term and nobody reads it: it is neither read nor integrated.
template <typename Real>
inline void pmlRowFace (const Real* __restrict cur, const Real* __restrict prev, Real* __restrict next,
                        const Real* __restrict oldX, const Real* __restrict oldY, const Real* __restrict oldZ,
                        Real* __restrict newX, Real* __restrict newY, Real* __restrict newZ,
                        const Real* __restrict zxc, const Real* __restrict zxf,
                        const Real* __restrict ax, const Real* __restrict bx,
                        int row, int i0, int i1, int sy, int sz,
                        Real zy, Real zz, Real zyP, Real zyM, Real zzP, Real zzM,
                        Real ayP, Real byP, Real ayM, Real byM, Real azP, Real bzP, Real azM, Real bzM,
                        Real cc, Real kd, Real halfK) noexcept
{
    for (int i = i0; i < i1; ++i)
    {
        const int x = row + i;
        const Real zx = zxc[i];
        const Real u = cur[x];
        const Real fxP = ax[i] * oldX[x] + bx[i] * (zy + zz - zxf[i]) * (cur[x + 1] - u);
        const Real fxM = ax[i - 1] * oldX[x - 1] + bx[i - 1] * (zy + zz - zxf[i - 1]) * (u - cur[x - 1]);
        const Real fyP = ayP * oldY[x] + byP * (zz + zx - zyP) * (cur[x + sy] - u);
        const Real fyM = ayM * oldY[x - sy] + byM * (zz + zx - zyM) * (u - cur[x - sy]);
        const Real fzP = azP * oldZ[x] + bzP * (zx + zy - zzP) * (cur[x + sz] - u);
        const Real fzM = azM * oldZ[x - sz] + bzM * (zx + zy - zzM) * (u - cur[x - sz]);
        newX[x] = fxP;
        newY[x] = fyP;
        newZ[x] = fzP;
        const Real s = zx + zy + zz;
        const Real lap = cur[x - 1] + cur[x + 1] + cur[x - sy] + cur[x + sy] + cur[x - sz] + cur[x + sz] - Real (6) * u;
        const Real div = (fxP - fxM) + (fyP - fyM) + (fzP - fzM);
        next[x] = (Real (2) * u - (Real (1) - halfK * s) * prev[x] + cc * lap + kd * div) / (Real (1) + halfK * s);
    }
}
} // namespace

// Time staggering: u and psi (= int u dt) at integer steps, phi at half
// steps; step n reads u^n, u^{n-1}, psi^n and phi^{n-1/2}, and writes u^{n+1},
// psi^{n+1} and phi^{n+1/2}, the leapfrog pair of the velocity-pressure Yee
// scheme. phi on a face shared by two cells is computed by both from the same
// old values and stored by its owner (the cell on its low side), so the cells
// can be updated in any order and by any thread.
template <typename Real>
void AirDomain<Real>::stepPml (long long n, int k0, int k1) noexcept
{
    if (pml == 0)
        return;
    const int border = pml + 1, sy = nx, sz = nx * ny;
    const Real cc = static_cast<Real> (c * c * k * k / (delta * delta));
    const Real kk = static_cast<Real> (k * k), kd = static_cast<Real> (k * k / delta), halfK = static_cast<Real> (0.5 * k);
    const auto oldP = static_cast<std::size_t> (n & 1), newP = static_cast<std::size_t> ((n + 1) & 1);
    for (int kz = std::max (k0, 1); kz < std::min (k1, nz - 1); ++kz)
        for (int j = 1; j < ny - 1; ++j)
        {
            int ranges[4];
            const int count = pmlRanges (nx, ny, nz, border, j, kz, ranges);
            const auto jj = static_cast<std::size_t> (j), kp = static_cast<std::size_t> (kz);
            const Real zy = zetaCellR[1][jj], zz = zetaCellR[2][kp];
            const auto full = [&] (int i0, int i1)
            {
                pmlRow<Real> (level (n), level (n - 1), level (n + 1), auxBase[oldP], auxBase[newP],
                              phiBase[oldP][0], phiBase[oldP][1], phiBase[oldP][2],
                              phiBase[newP][0], phiBase[newP][1], phiBase[newP][2],
                              zetaCellR[0].data(), zetaFaceR[0].data(), faceDecay[0].data(), faceGain[0].data(),
                              index (0, j, kz), i0, i1, sy, sz,
                              zy, zz, zetaFaceR[1][jj], zetaFaceR[1][jj - 1],
                              zetaFaceR[2][kp], zetaFaceR[2][kp - 1], faceDecay[1][jj], faceGain[1][jj],
                              faceDecay[1][jj - 1], faceGain[1][jj - 1], faceDecay[2][kp], faceGain[2][kp],
                              faceDecay[2][kp - 1], faceGain[2][kp - 1], cc, kk, kd, halfK);
            };
            const auto face = [&] (int i0, int i1)
            {
                pmlRowFace<Real> (level (n), level (n - 1), level (n + 1),
                                  phiBase[oldP][0], phiBase[oldP][1], phiBase[oldP][2],
                                  phiBase[newP][0], phiBase[newP][1], phiBase[newP][2],
                                  zetaCellR[0].data(), zetaFaceR[0].data(), faceDecay[0].data(), faceGain[0].data(),
                                  index (0, j, kz), i0, i1, sy, sz,
                                  zy, zz, zetaFaceR[1][jj], zetaFaceR[1][jj - 1],
                                  zetaFaceR[2][kp], zetaFaceR[2][kp - 1], faceDecay[1][jj], faceGain[1][jj],
                                  faceDecay[1][jj - 1], faceGain[1][jj - 1], faceDecay[2][kp], faceGain[2][kp],
                                  faceDecay[2][kp - 1], faceGain[2][kp - 1], cc, kd, halfK);
            };
            // Cells with at most one non-zero profile skip psi. Along a row
            // zx is non-zero on [1, border) and [nx - border, nx - 1) only.
            const bool rowHasY = zy != Real (0), rowHasZ = zz != Real (0);
            for (int r = 0; r < count; ++r)
            {
                const int i0 = ranges[2 * r], i1 = ranges[2 * r + 1];
                if (! rowHasY && ! rowHasZ)
                {
                    face (i0, i1);
                    continue;
                }
                if (rowHasY && rowHasZ)
                {
                    full (i0, i1);
                    continue;
                }
                // One of zy, zz: the x-layer parts are edges, the rest a face.
                const int lowEnd = std::min (i1, border), highStart = std::max (i0, nx - border);
                if (lowEnd > i0)
                    full (i0, lowEnd);
                if (highStart > std::max (i0, lowEnd))
                    face (std::max (i0, lowEnd), highStart);
                if (i1 > highStart)
                    full (std::max (highStart, i0), i1);
            }
        }
}

namespace
{
// Trilinear interpolation of a cell-centred field (or of a difference of two).
template <typename Real, typename Domain>
double trilinear (const Domain& air, const Real* a, const Real* b, const std::array<double, 3>& point) noexcept
{
    const double fx = (point[0] - air.x0) / air.delta - 0.5, fy = (point[1] - air.y0) / air.delta - 0.5;
    const double fz = (point[2] - air.z0) / air.delta - 0.5;
    const int i = static_cast<int> (std::floor (fx)), j = static_cast<int> (std::floor (fy)), kz = static_cast<int> (std::floor (fz));
    // A point outside the domain (a test's small box) reads as silence.
    if (i < 0 || j < 0 || kz < 0 || i + 1 >= air.nx || j + 1 >= air.ny || kz + 1 >= air.nz)
        return 0.0;
    const double tx = fx - i, ty = fy - j, tz = fz - kz;
    double sum = 0.0;
    for (int dz = 0; dz < 2; ++dz)
        for (int dy = 0; dy < 2; ++dy)
            for (int dx = 0; dx < 2; ++dx)
            {
                const double weight = (dx ? tx : 1.0 - tx) * (dy ? ty : 1.0 - ty) * (dz ? tz : 1.0 - tz);
                const int cell = air.index (i + dx, j + dy, kz + dz);
                sum += weight * (static_cast<double> (a[cell]) - (b ? static_cast<double> (b[cell]) : 0.0));
            }
    return sum;
}
} // namespace

template <typename Real>
double AirDomain<Real>::pressure (long long n, const std::array<double, 3>& point) const noexcept
{
    return density * trilinear<Real> (*this, level (n + 1), level (n - 1), point) / (2.0 * k);
}

template <typename Real>
double AirDomain<Real>::potential (long long n, const std::array<double, 3>& point) const noexcept
{
    return trilinear<Real> (*this, level (n), static_cast<const Real*> (nullptr), point);
}

template <typename Real>
double AirDomain<Real>::energy (long long n) const
{
    const Real* next = level (n + 1);
    const Real* now = level (n);
    const auto counted = [this] (int cell)
    {
        const CellType t = types[static_cast<std::size_t> (cell)];
        return t == CellType::Regular || t == CellType::Irregular || t == CellType::Slice;
    };
    const int stride[3] = { 1, nx, nx * ny };
    double kinetic = 0.0, potential = 0.0;
    for (int cell = 0; cell < static_cast<int> (cells()); ++cell)
    {
        if (! counted (cell))
            continue;
        const double d = static_cast<double> (next[cell]) - now[cell];
        kinetic += d * d;
        const std::uint8_t mask = masks[static_cast<std::size_t> (cell)];
        for (int dir = 0; dir < 3; ++dir)
        {
            const int nb = cell + stride[dir];
            if ((mask & (1u << (2 * dir + 1))) && counted (nb))
                potential += (static_cast<double> (next[nb]) - next[cell]) * (static_cast<double> (now[nb]) - now[cell]);
        }
    }
    return density * delta * delta * delta / (2.0 * c * c * k * k) * kinetic + 0.5 * density * delta * potential;
}

#define TAIKOR_AIR_INSTANTIATE(Real)                                                                \
    template void AirDomain<Real>::stepPml (long long, int, int) noexcept;                          \
    template double AirDomain<Real>::pressure (long long, const std::array<double, 3>&) const noexcept; \
    template double AirDomain<Real>::potential (long long, const std::array<double, 3>&) const noexcept; \
    template double AirDomain<Real>::energy (long long) const;
TAIKOR_AIR_INSTANTIATE (float)
TAIKOR_AIR_INSTANTIATE (double)
#undef TAIKOR_AIR_INSTANTIATE

} // namespace taikor::candidates::fdtd
