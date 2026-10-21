// Test 5 (first part): the absorbing layer's reflection on the production
// air grid, for a short (2-cell) and a longer (4-cell) Gaussian pulse at
// normal and oblique incidence, against a reference domain large enough that
// nothing returns within the recording. pmlScan (run only when
// TAIKOR_FDTD_SCAN is set) tabulates thickness against theoretical R.

#include "TestSupport.h"

#include "AirChecks.h"
#include "FdtdDrumConfig.h"
#include "ThreadPool.h"

#include <chrono>
#include <cstdlib>

#include "AirDomain.h"

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::fdtd::tests;
} // namespace

FDTD_TEST (pmlReflection)
{
    const auto& drum = nagado().drum;
    const auto numerics = FdtdNumerics::fromEnvironment();
    const double k = 1.0 / (48000.0 * numerics.oversample);
    const double delta = drum.depth / numerics.airCellsPerDepth;
    for (const double blob : { 2.0, 4.0 })
    {
        const auto r = measurePmlReflection (delta, k, drum.soundSpeed, drum.airDensity, numerics.pmlCells,
                                             numerics.pmlReflection, blob, std::min (8, candidateThreadCount()));
        report.check (format ("PML reflection, %.0f-cell pulse", blob), r.peakReflectionDb <= -40.0,
                      format ("%.1f dB peak, %.1f dB energy (%d cells, theoretical R %.0e; %.1f s), limit -40 dB peak",
                              r.peakReflectionDb, r.energyReflectionDb, numerics.pmlCells, numerics.pmlReflection,
                              r.seconds));
    }
}

FDTD_TEST (pmlScan)
{
    if (! std::getenv ("TAIKOR_FDTD_SCAN"))
        return;
    const auto& drum = nagado().drum;
    const auto numerics = FdtdNumerics::fromEnvironment();
    const double k = 1.0 / (48000.0 * numerics.oversample);
    const double delta = drum.depth / numerics.airCellsPerDepth;
    if (std::getenv ("TAIKOR_FDTD_SCAN_ONE"))
    {
        const double blob = std::atof (std::getenv ("TAIKOR_FDTD_SCAN_ONE"));
        const int cells = std::getenv ("TAIKOR_FDTD_SCAN_CELLS") ? std::atoi (std::getenv ("TAIKOR_FDTD_SCAN_CELLS")) : 10;
        const double reflection = std::getenv ("TAIKOR_FDTD_SCAN_R") ? std::atof (std::getenv ("TAIKOR_FDTD_SCAN_R")) : 1e-3;
        const auto b = measurePmlReflection (delta, k, drum.soundSpeed, drum.airDensity, cells, reflection, blob, 8);
        report.note (format ("one: %.1f / %.1f dB", b.peakReflectionDb, b.energyReflectionDb));
        return;
    }
    for (const int cells : { 4, 6, 8, 10, 14 })
        for (const double reflection : { 1e-2, 1e-3, 1e-4, 1e-6 })
        {
            const auto a = measurePmlReflection (delta, k, drum.soundSpeed, drum.airDensity, cells, reflection, 2.0, 8);
            const auto b = measurePmlReflection (delta, k, drum.soundSpeed, drum.airDensity, cells, reflection, 4.0, 8);
            report.note (format ("cells %2d R %.0e: 2-cell pulse %.1f dB peak %.1f dB energy; 4-cell %.1f / %.1f dB",
                                 cells, reflection, a.peakReflectionDb, a.energyReflectionDb, b.peakReflectionDb,
                                 b.energyReflectionDb));
        }
}

// Development timing (TAIKOR_FDTD_SCAN set): the layer kernel per cell on
// full planes against the strip planes, and the regular kernel.
FDTD_TEST (pmlBench)
{
    if (! std::getenv ("TAIKOR_FDTD_SCAN"))
        return;
    const auto& drum = nagado().drum;
    const auto numerics = FdtdNumerics::fromEnvironment();
    AirGeometry g;
    g.radius = drum.radius;
    g.depth = drum.depth;
    g.wall = drum.shellWallThickness;
    g.delta = drum.depth / numerics.airCellsPerDepth;
    g.margin = numerics.airMargin;
    g.pmlCells = numerics.pmlCells;
    g.pmlReflection = numerics.pmlReflection;
    AirDomain<float> air;
    air.build (g, drum.soundSpeed, drum.airDensity, 1.0 / (48000.0 * numerics.oversample));
    const int border = numerics.pmlCells + 1;
    const auto time = [&] (auto&& f)
    {
        const auto t0 = std::chrono::steady_clock::now();
        for (int n = 0; n < 50; ++n)
            f (n);
        return std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count() / 50.0;
    };
    const double full = time ([&] (int n) { air.stepPml (n, 1, border); air.stepPml (n, air.nz - border, air.nz - 1); });
    const double strips = time ([&] (int n) { air.stepPml (n, border, air.nz - border); });
    const double regular = time ([&] (int n) { air.stepRegular (n, 0, air.nz); });
    const double fullCells = 2.0 * (border - 1) * (air.nx - 2.0) * (air.ny - 2.0);
    const double stripCells = static_cast<double> (air.pmlCellCount()) - fullCells;
    report.note (format ("PML full planes %.2f ns/cell (%.0f cells), strip planes %.2f ns/cell (%.0f cells); regular %.2f ns/cell",
                         1e9 * full / fullCells, fullCells, 1e9 * strips / stripCells, stripCells,
                         1e9 * regular / (air.cells() - air.pmlCellCount())));
}
