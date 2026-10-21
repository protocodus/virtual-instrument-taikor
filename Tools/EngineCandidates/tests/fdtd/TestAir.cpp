// Test 2: the air alone in the barrel closed by two rigid heads, on the
// production air grid: resonances against the rigid closed cylinder
// f = c/(2 pi) sqrt((j'_mq / a)^2 + (p pi / L)^2).

#include "TestSupport.h"

#include "AirDomain.h"
#include "Analytic.h"
#include "FdtdDrumConfig.h"
#include "HeadCalibration.h"
#include "ThreadPool.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::fdtd::tests;
} // namespace

FDTD_TEST (cavityModes)
{
    const auto& drum = nagado().drum;
    const auto numerics = FdtdNumerics::fromEnvironment();
    const double k = 1.0 / (48000.0 * numerics.oversample);
    AirGeometry geometry;
    geometry.radius = drum.radius;
    geometry.depth = drum.depth;
    geometry.wall = drum.shellWallThickness;
    geometry.delta = drum.depth / numerics.airCellsPerDepth;
    geometry.margin = geometry.delta;
    geometry.pmlCells = 0;
    AirDomain<double> air;
    air.build (geometry, drum.soundSpeed, drum.airDensity, k);
    const auto columns = air.headColumns (drum.radius);
    air.closeHeadFaces (air.planeAbove (0.0), columns, false);
    air.closeHeadFaces (air.planeAbove (-drum.depth), columns, false);
    air.finalise();
    // Rest at a smooth off-axis blob (sigma = 2 cells), so every family is excited.
    const std::array<double, 3> centre { 0.17, 0.09, -0.31 };
    const double s = 2.0 * geometry.delta;
    for (int kz = 0; kz < air.nz; ++kz)
        for (int j = 0; j < air.ny; ++j)
            for (int i = 0; i < air.nx; ++i)
            {
                const int cell = air.index (i, j, kz);
                if (air.type (cell) == CellType::Solid || std::hypot (air.x (i), air.y (j)) >= drum.radius
                    || air.z (kz) > 0.0 || air.z (kz) < -drum.depth)
                    continue;
                const double dx = air.x (i) - centre[0], dy = air.y (j) - centre[1], dz = air.z (kz) - centre[2];
                const double value = std::exp (-(dx * dx + dy * dy + dz * dz) / (2.0 * s * s));
                air.level (0)[cell] = air.level (-1)[cell] = value;
            }
    const std::vector<std::array<double, 3>> receivers { { 0.21, -0.13, -0.08 }, { -0.05, 0.3, -0.62 }, { 0.33, 0.02, -0.47 } };
    const int steps = static_cast<int> (0.3 / k);
    std::vector<std::vector<double>> traces (receivers.size());
    ThreadPool pool (std::min (8, candidateThreadCount()));
    Barrier barrier (pool.size());
    pool.run ([&] (int t)
    {
        int k0, k1, i0, i1;
        splitRange (air.nz, pool.size(), t, k0, k1);
        splitRange (air.irregularCount(), pool.size(), t, i0, i1);
        for (long long n = 0; n < steps; ++n)
        {
            air.stepRegular (n, k0, k1);
            air.stepIrregular (n, i0, i1);
            barrier.wait();
            // Level n+1 is not written again until step n+3.
            if (t == 0)
                for (std::size_t r = 0; r < receivers.size(); ++r)
                    traces[r].push_back (air.potential (n + 1, receivers[r]));
        }
    });
    if (const char* dump = std::getenv ("TAIKOR_FDTD_DUMP"))
        if (FILE* file = std::fopen (dump, "wb"))
        {
            for (const auto& trace : traces)
                std::fwrite (trace.data(), sizeof (double), trace.size(), file);
            std::fclose (file);
        }
    const auto modes = rigidCylinderModes (drum.radius, drum.depth, drum.soundSpeed, 620.0);
    // Each analytic mode is searched for between the midpoints to its analytic
    // neighbours (capped at +/- 3 %); a maximum on the window edge means the
    // mode was not resolved separately and is skipped.
    double worst = 0.0, worstAxial = 0.0;
    int measured = 0;
    std::string table;
    for (std::size_t m = 0; m < modes.size(); ++m)
    {
        const double f0 = modes[m].frequency;
        if (f0 < 1.0)
            continue;
        double lo = 0.97 * f0, hi = 1.03 * f0;
        for (std::size_t o = 0; o < modes.size(); ++o)
        {
            const double f1 = modes[o].frequency;
            if (o == m || std::abs (f1 - f0) < 1e-9)
                continue;
            if (f1 < f0)
                lo = std::max (lo, 0.5 * (f0 + f1));
            else
                hi = std::min (hi, 0.5 * (f0 + f1));
        }
        const double f = spectralPeakFrequency (traces, k, lo, hi);
        const bool edge = f < lo + 0.3 || f > hi - 0.3;
        const double error = f / f0 - 1.0;
        if (! edge)
        {
            worst = std::max (worst, std::abs (error));
            if (modes[m].m == 0 && modes[m].q == 0)
                worstAxial = std::max (worstAxial, std::abs (error));
            ++measured;
        }
        table += format ("(%d,%d,%d) %.2f/%.2f Hz %s; ", modes[m].m, modes[m].q, modes[m].p, f, f0,
                         edge ? "unresolved" : format ("%+.3f%%", 100.0 * error).c_str());
    }
    report.note (format ("air grid %dx%dx%d, Delta %.3f mm, Courant %.4f", air.nx, air.ny, air.nz,
                         1e3 * geometry.delta, air.courant()));
    report.note (table);
    // Axial modes test the head planes (L / Delta integer, heads midway
    // between air planes); the transverse ones also carry the staircase
    // barrel, whose jagged Neumann wall acts 0.4-1 % wide at a / Delta ~ 31
    // (its area is within 0.13 % of pi a^2; the longer perimeter lowers the
    // modes, Weyl's second term). That error is reported, not corrected.
    report.check ("rigid cylinder axial modes", worstAxial < 0.001,
                  format ("worst |error| %.4f %%, limit 0.1 %%", 100.0 * worstAxial));
    report.check ("rigid cylinder modes below 620 Hz", worst < 0.015 && measured >= 10,
                  format ("worst |error| %.3f %% over %d resolved modes (staircase barrel), limit 1.5 %%", 100.0 * worst, measured));
}
