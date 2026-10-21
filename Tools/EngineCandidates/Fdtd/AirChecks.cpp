#include "AirChecks.h"

#include "AirDomain.h"
#include "ThreadPool.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace taikor::candidates::fdtd
{
namespace
{
// Runs a drum-free box (half-width `half` cells around the source, `pml`
// layer cells) from a Gaussian blob at rest and records Psi at receivers.
std::vector<std::vector<double>> runBox (double delta, double k, double c, double rho, int half, int pml,
                                         double reflection, double blobCells, int steps,
                                         const std::vector<std::array<double, 3>>& receivers, int threads)
{
    AirGeometry g;
    g.radius = 0.0;
    g.wall = 0.0;
    g.depth = delta;               // one cell "barrel" of zero radius: no solid
    g.delta = delta;
    g.margin = half * delta;
    g.pmlCells = pml;
    g.pmlReflection = reflection;
    AirDomain<double> air;
    air.build (g, c, rho, k);
    const double s = blobCells * delta;
    for (int kz = 0; kz < air.nz; ++kz)
        for (int j = 0; j < air.ny; ++j)
            for (int i = 0; i < air.nx; ++i)
            {
                const int cell = air.index (i, j, kz);
                if (air.type (cell) != CellType::Regular && air.type (cell) != CellType::Irregular)
                    continue;
                const double r2 = air.x (i) * air.x (i) + air.y (j) * air.y (j) + air.z (kz) * air.z (kz);
                air.level (0)[cell] = air.level (-1)[cell] = std::exp (-r2 / (2.0 * s * s));
            }
    std::vector<std::vector<double>> traces (receivers.size());
    ThreadPool pool (threads);
    Barrier barrier (threads);
    pool.run ([&] (int t)
    {
        int k0, k1, i0, i1;
        splitRange (air.nz, threads, t, k0, k1);
        splitRange (air.irregularCount(), threads, t, i0, i1);
        for (long long n = 0; n < steps; ++n)
        {
            air.stepRegular (n, k0, k1);
            air.stepIrregular (n, i0, i1);
            air.stepPml (n, k0, k1);
            barrier.wait();
            if (t == 0)
                for (std::size_t r = 0; r < receivers.size(); ++r)
                    traces[r].push_back (air.potential (n + 1, receivers[r]));
            barrier.wait();
        }
    });
    return traces;
}
} // namespace

ReflectionResult measurePmlReflection (double delta, double k, double c, double rho, int pmlCells,
                                       double pmlReflection, double blobCells, int threads)
{
    const auto t0 = std::chrono::steady_clock::now();
    // Source at the origin, layer from `half` cells out; receivers 6 cells
    // inside the layer's face, on the axis and at 45 and 55 degrees.
    const int half = static_cast<int> (std::ceil (8.0 * blobCells)) + 12;
    const double inner = (half - 6) * delta;
    const std::vector<std::array<double, 3>> receivers {
        { inner, 0.0, 0.0 }, { inner, inner, 0.0 }, { inner, inner, inner }, { 0.0, 0.0, inner }
    };
    // Record until the incident blob and the layer's round trip (6 + pml
    // cells each way, the blob's +/- 4 sigma) have passed the farthest receiver.
    double windowCells = 0.0;
    for (const auto& r : receivers)
        windowCells = std::max (windowCells, std::hypot (r[0], r[1], r[2]) / delta + 2.0 * (6 + pmlCells) + 8.0 * blobCells);
    const int steps = static_cast<int> (std::ceil (windowCells * delta / (c * k)));
    // Reference box: the leading edge (4 sigma ahead of the centre) of its
    // nearest image source (a wall reflection) must reach every receiver
    // after the window.
    int far = half;
    for (bool enough = false; ! enough; ++far)
    {
        enough = true;
        for (const auto& r : receivers)
            for (int axis = 0; axis < 3; ++axis)
            {
                auto image = r;
                image[static_cast<std::size_t> (axis)] = 2.0 * far * delta - r[static_cast<std::size_t> (axis)];
                const double d = std::hypot (image[0], image[1], image[2]) / delta;
                enough = enough && d >= windowCells + 4.0 * blobCells + 4.0;
            }
    }
    const auto layered = runBox (delta, k, c, rho, half, pmlCells, pmlReflection, blobCells, steps, receivers, threads);
    const auto reference = runBox (delta, k, c, rho, far, 0, pmlReflection, blobCells, steps, receivers, threads);
    if (const char* dump = std::getenv ("TAIKOR_FDTD_DUMP"))
        if (FILE* file = std::fopen (dump, "wb"))
        {
            for (const auto* set : { &layered, &reference })
                for (const auto& trace : *set)
                    std::fwrite (trace.data(), sizeof (double), trace.size(), file);
            std::fclose (file);
        }
    ReflectionResult result;
    double peakIncident = 0.0, peakDifference = 0.0, energyIncident = 0.0, energyDifference = 0.0;
    for (std::size_t r = 0; r < receivers.size(); ++r)
        for (int n = 1; n + 1 < steps; ++n)
        {
            // Pressure ~ d Psi / dt: centred differences of the recorded Psi.
            const double incident = reference[r][static_cast<std::size_t> (n + 1)] - reference[r][static_cast<std::size_t> (n - 1)];
            const double difference = (layered[r][static_cast<std::size_t> (n + 1)] - layered[r][static_cast<std::size_t> (n - 1)]) - incident;
            peakIncident = std::max (peakIncident, std::abs (incident));
            peakDifference = std::max (peakDifference, std::abs (difference));
            energyIncident += incident * incident;
            energyDifference += difference * difference;
        }
    result.peakReflectionDb = 20.0 * std::log10 (peakDifference / peakIncident);
    result.energyReflectionDb = 10.0 * std::log10 (energyDifference / energyIncident);
    result.seconds = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
    return result;
}

} // namespace taikor::candidates::fdtd
