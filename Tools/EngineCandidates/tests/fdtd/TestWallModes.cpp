// Test (wall a): decay of the sealed barrel's duct modes under the Kirchhoff
// side-wall boundary layer (WallLoss.h), rigid ends, production air grid,
// against the continuum rate of each rigid-cylinder mode
//   alpha = (c^2/2) sqrt(w nu/2) [oint_side |grad_t psi|^2 / w^2
//           + (gamma - 1) / (c^2 sqrt Pr) oint_side psi^2] / int psi^2 dV
// (Pierce, Acoustics, sec. 10.5; the modal engine's wallLossRate), evaluated
// at the measured discrete frequency.
//
// Method: modes (m, q, p) with distinct p are seeded together as the analytic
// shapes J_m(j'_mq r / a) cos(m theta) cos(p pi (z + L) / L); the barrel is
// exactly separable in z on this grid (the ends are cell faces, the wall
// terms are the same on every plane), so each shape's projection holds its
// own p alone. Five runs - no loss, thermal only, viscous only, both, and
// both at a tenth of their strength - and the ratio of each mode's energy
// envelope lossy / lossless, whose log falls at 2 alpha (the ratio cancels
// the beating of the staircase's slightly mixed transverse modes, common to
// both runs).
//
// Kirchhoff's rate is first order in the wall admittance. At full strength a
// plane wave (0, 0, p) also couples through the wall to its transverse
// neighbours (0, q, p), which lie only c^2 k_r^2 / 2w away (140 rad/s at
// 6.5 kHz against rates of ~2 1/s), and that second-order term adds +2.5 %
// at 4.2 kHz and +4.1 % at 6.6 kHz; it scales with the square of the strength
// (physics of any wall this lossy, not of the grid), so the first-order rate is
// read from the run at a tenth of the strength. Measured (96 kHz, 0.2 s):
// axial modes 1.1-6.6 kHz within 0.4 % to first order; single transverse
// modes (0,2,1) -0.5 %, (0,5,2) +4.9 % (the staircase gives a radial flow a
// little shear it does not have), (5,0,9) -1.3 %, (1,0,10) +0.2 %, (4,0,8)
// -2.7 %.

#include "WallSupport.h"

#include "HeadCalibration.h"
#include "SpecialFunctions.h"

#include <algorithm>
#include <cmath>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::fdtd::tests;
constexpr double pi = 3.14159265358979323846;

struct DuctMode
{
    int m, q, p;
    const char* family;
    double jPrime = 0.0, frequency = 0.0;
    std::vector<double> column;  // the transverse shape on the interior columns
};

// The continuum Kirchhoff rate of mode (m, q, p), viscous and thermal parts,
// with the boundary layer at angular frequency w. The viscous part is written
// with the mode's own wavenumber, alpha_v = (1/2) sqrt(w nu/2) oint |grad_t psi|^2
// / (k^2 int psi^2) (for the continuum mode k = w / c): a discrete mode keeps
// its shape's ratio of wall shear to kinetic energy while its frequency moves
// with the air grid's dispersion (-11 % for the plane wave p = 40 on axis), so
// its reference is that ratio at the frequency it rings at.
void kirchhoff (const DuctMode& mode, const DrumSpec& drum, double w, double& viscous, double& thermal)
{
    const Viscothermal gas;
    const double a = drum.radius, len = drum.depth;
    const double nu = gas.viscosity / drum.airDensity;
    const double j = mode.m == 0 && mode.q == 0 ? 1.0 : special::besselJ (mode.m, mode.jPrime);
    const double j2 = j * j;
    const double thetaCos = mode.m == 0 ? 2.0 * pi : pi, thetaSin = mode.m == 0 ? 0.0 : pi;
    const double zCos = mode.p == 0 ? len : 0.5 * len, zSin = mode.p == 0 ? 0.0 : 0.5 * len;
    const double kz = mode.p * pi / len, kr = mode.jPrime / a;
    const double tangential = a * j2 * (double (mode.m * mode.m) / (a * a) * thetaSin * zCos + kz * kz * thetaCos * zSin);
    const double pressure = a * j2 * thetaCos * zCos;
    const double ratio = mode.jPrime > 0.0 ? double (mode.m * mode.m) / (mode.jPrime * mode.jPrime) : 0.0;
    const double volume = thetaCos * zCos * 0.5 * a * a * (1.0 - ratio) * j2;
    const double root = std::sqrt (w * nu / 2.0);
    viscous = 0.5 * root * tangential / ((kr * kr + kz * kz) * volume);
    thermal = 0.5 * root * (gas.gamma - 1.0) / std::sqrt (gas.prandtl) * pressure / volume;
}
} // namespace

FDTD_TEST (wallDuctModes)
{
    const auto& drum = nagado().drum;
    const auto numerics = FdtdNumerics::fromEnvironment();
    const double k = 1.0 / (48000.0 * numerics.oversample);
    // Distinct p throughout. On the square grid the staircase couples m to
    // m +- 4, 8, ...; where such a neighbour lies within tens of hertz the
    // analytic shape excites two discrete modes (the check below reports
    // those as mixed and leaves them out).
    std::vector<DuctMode> modes { { 0, 0, 6, "axial" },       { 0, 0, 12, "axial" },     { 0, 0, 24, "axial" },
                                  { 0, 0, 40, "axial" },      { 0, 2, 1, "radial" },     { 0, 5, 2, "radial" },
                                  { 0, 9, 3, "radial" },      { 5, 0, 9, "azimuthal" },  { 15, 0, 13, "azimuthal" },
                                  { 8, 0, 4, "azimuthal" },   { 16, 0, 5, "azimuthal" }, { 32, 0, 7, "azimuthal" },
                                  { 1, 0, 10, "mixed" },      { 4, 0, 8, "mixed" } };
    const auto zeros = special::besselJPrimeZeros (40, 60.0);
    for (auto& mode : modes)
    {
        mode.jPrime = zeros[static_cast<std::size_t> (mode.m)][static_cast<std::size_t> (mode.q)];
        const double kr = mode.jPrime / drum.radius, kz = mode.p * pi / drum.depth;
        mode.frequency = drum.soundSpeed * std::sqrt (kr * kr + kz * kz) / (2.0 * pi);
    }
    const double seconds = 0.2, skip = 0.02;
    const long long steps = static_cast<long long> (seconds / k);
    // traces[run][mode]: projections after each step.
    std::vector<std::vector<std::vector<double>>> traces (5);
    std::string layout;
    for (int run = 0; run < 5; ++run)
    {
        SealedBarrel barrel;
        barrel.build (drum, k, run > 0);
        if (run == 1)
            barrel.wall.scale (1.0, 0.0);
        if (run == 2)
            barrel.wall.scale (0.0, 1.0);
        if (run == 4)
            barrel.wall.scale (0.1, 0.1);
        auto& air = barrel.air;
        // Interior columns and planes.
        std::vector<int> columns;
        for (int j = 0; j < air.ny; ++j)
            for (int i = 0; i < air.nx; ++i)
                if (std::hypot (air.x (i), air.y (j)) < drum.radius && air.type (air.index (i, j, air.nz / 2)) != CellType::Solid)
                    columns.push_back (j * air.nx + i);
        std::vector<int> planes;
        for (int kz = 0; kz < air.nz; ++kz)
            if (air.z (kz) < 0.0 && air.z (kz) > -drum.depth)
                planes.push_back (kz);
        if (run == 3)
            layout = format ("%zu columns x %zu planes, Delta %.3f mm, %d wall cells, %d friction elements", columns.size(),
                             planes.size(), 1e3 * air.delta, barrel.wall.wallCellCount(), barrel.wall.frictionCount());
        for (auto& mode : modes)
        {
            mode.column.resize (columns.size());
            for (std::size_t c = 0; c < columns.size(); ++c)
            {
                const double x = air.x (columns[c] % air.nx), y = air.y (columns[c] / air.nx);
                mode.column[c] = special::besselJ (mode.m, mode.jPrime * std::hypot (x, y) / drum.radius)
                                 * std::cos (mode.m * std::atan2 (y, x));
            }
        }
        const auto zShape = [&] (const DuctMode& mode, int kz) { return std::cos (mode.p * pi * (air.z (kz) + drum.depth) / drum.depth); };
        for (const int kz : planes)
            for (std::size_t c = 0; c < columns.size(); ++c)
            {
                double value = 0.0;
                for (const auto& mode : modes)
                    value += mode.column[c] * zShape (mode, kz);
                const int cell = kz * air.nx * air.ny + columns[c];
                air.level (0)[cell] = air.level (-1)[cell] = value;
            }
        const int threads = std::min (8, candidateThreadCount());
        std::vector<std::vector<double>> partial (static_cast<std::size_t> (threads), std::vector<double> (modes.size()));
        auto& out = traces[static_cast<std::size_t> (run)];
        out.assign (modes.size(), {});
        barrel.run (steps, [&] (int t, int count, long long n)
        {
            int p0, p1;
            splitRange (static_cast<int> (planes.size()), count, t, p0, p1);
            auto& sums = partial[static_cast<std::size_t> (t)];
            std::fill (sums.begin(), sums.end(), 0.0);
            const double* level = air.level (n + 1);
            for (int pl = p0; pl < p1; ++pl)
            {
                const int kz = planes[static_cast<std::size_t> (pl)];
                const double* plane = level + static_cast<std::size_t> (kz) * air.nx * air.ny;
                for (std::size_t m = 0; m < modes.size(); ++m)
                {
                    double dot = 0.0;
                    for (std::size_t c = 0; c < columns.size(); ++c)
                        dot += modes[m].column[c] * plane[columns[c]];
                    sums[m] += dot * zShape (modes[m], kz);
                }
            }
        }, [&] (long long)
        {
            for (std::size_t m = 0; m < modes.size(); ++m)
            {
                double sum = 0.0;
                for (const auto& p : partial)
                    sum += p[m];
                out[m].push_back (sum);
            }
        }, threads);
    }
    report.note (layout);
    const long long first = static_cast<long long> (skip / k);
    // Slope of ln(E_run / E_lossless) against t: -2 alpha.
    const auto rate = [&] (int run, std::size_t m, double w, long long from, long long to)
    {
        const auto envelope = [&] (const std::vector<double>& a, long long n)
        {
            // Constant for a discrete sinusoid at w: a_n = cos (w n k) gives 1.
            const double d = (a[static_cast<std::size_t> (n + 1)] - a[static_cast<std::size_t> (n - 1)]) / (2.0 * std::sin (w * k));
            return d * d + a[static_cast<std::size_t> (n)] * a[static_cast<std::size_t> (n)];
        };
        double st = 0.0, sy = 0.0, stt = 0.0, sty = 0.0, count = 0.0;
        for (long long n = from; n < to; ++n)
        {
            const double t = n * k;
            const double y = std::log (envelope (traces[static_cast<std::size_t> (run)][m], n) / envelope (traces[0][m], n));
            st += t;
            sy += y;
            stt += t * t;
            sty += t * y;
            count += 1.0;
        }
        return -0.5 * (count * sty - st * sy) / (count * stt - st * st);
    };
    // A mode counts as one discrete mode when its lossless projection's
    // envelope stays flat (no beating with a staircase-mixed neighbour).
    const auto beating = [&] (std::size_t m, double w)
    {
        double lo = 1e300, hi = 0.0, sum = 0.0, count = 0.0;
        for (long long n = first; n < steps - 1; ++n)
        {
            const auto& a = traces[0][m];
            const double d = (a[static_cast<std::size_t> (n + 1)] - a[static_cast<std::size_t> (n - 1)]) / (2.0 * std::sin (w * k));
            const double e = d * d + a[static_cast<std::size_t> (n)] * a[static_cast<std::size_t> (n)];
            lo = std::min (lo, e);
            hi = std::max (hi, e);
            sum += e;
            count += 1.0;
        }
        return (hi - lo) / (sum / count);
    };
    double worst[3] = {}, worstWeak[3] = {}; // axial, transverse, and all single modes: worst |error|
    int single = 0;
    std::string table;
    for (std::size_t m = 0; m < modes.size(); ++m)
    {
        const double f = spectralPeakFrequency (traces[0][m], k, 0.8 * modes[m].frequency, 1.03 * modes[m].frequency);
        const double w = 2.0 * pi * f;
        double viscous, thermal;
        kirchhoff (modes[m], drum, w, viscous, thermal);
        const double full = rate (3, m, w, first, steps - 1), th = rate (1, m, w, first, steps - 1),
                     vi = rate (2, m, w, first, steps - 1), weak = 10.0 * rate (4, m, w, first, steps - 1);
        const double depth = beating (m, w);
        const double error = full / (viscous + thermal) - 1.0, errorWeak = weak / (viscous + thermal) - 1.0;
        const bool axial = modes[m].m == 0 && modes[m].q == 0, pure = depth < 0.1;
        if (pure)
        {
            ++single;
            const int family = axial ? 0 : 1;
            worst[family] = std::max (worst[family], std::abs (error));
            worstWeak[family] = std::max (worstWeak[family], std::abs (errorWeak));
        }
        table += format ("(%d,%d,%d) %s %.1f Hz, beat %.2f%s: %.4f 1/s (%.2f dB/s) vs Kirchhoff %.4f (%+.1f %%), first order %+.1f %%; "
                         "thermal %.4f vs %.4f (%+.1f %%), viscous %.4f vs %.4f (%+.1f %%)\n  ",
                         modes[m].m, modes[m].q, modes[m].p, modes[m].family, f, depth, pure ? "" : " (mixed)", full,
                         20.0 / std::log (10.0) * full, viscous + thermal, 100.0 * error, 100.0 * errorWeak, th, thermal,
                         100.0 * (th / thermal - 1.0), vi, viscous, viscous > 0.0 ? 100.0 * (vi / viscous - 1.0) : 0.0);
    }
    report.note (table);
    report.check ("wall loss: axial duct modes 1-6.6 kHz at the Kirchhoff rate (first order)", worstWeak[0] < 0.015,
                  format ("worst |error| %.2f %% at 1/10 strength (%.2f %% at full strength, where the wall's coupling of "
                          "each plane wave to its transverse neighbours adds a second-order term), limit 1.5 %%",
                          100.0 * worstWeak[0], 100.0 * worst[0]));
    report.check ("wall loss: single transverse duct modes at the Kirchhoff rate", worst[1] < 0.1 && single >= 8,
                  format ("worst |error| %.1f %% over the single modes (first order %.1f %%), limit 10 %%; %d of %zu modes single",
                          100.0 * worst[1], 100.0 * worstWeak[1], single, modes.size()));
}
