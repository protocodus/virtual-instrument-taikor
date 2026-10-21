#pragma once

// The sealed barrel for the wall-loss tests: the production air grid's
// interior, closed by two rigid ends where the heads would be, with the
// Kirchhoff side-wall boundary layer (WallLoss.h), stepped on a thread pool.

#include "TestSupport.h"

#include "AirDomain.h"
#include "FdtdDrumConfig.h"
#include "ThreadPool.h"
#include "WallLoss.h"

#include <cmath>
#include <functional>
#include <memory>
#include <vector>

namespace taikor::candidates::fdtd::tests
{

struct SealedBarrel
{
    AirDomain<double> air;
    WallLoss<double> wall;
    SquareRootFit fit;
    double k = 0.0, radius = 0.0, depth = 0.0;
    bool withWall = true;

    void build (const DrumSpec& drum, double timeStep, bool wallLoss)
    {
        const auto numerics = FdtdNumerics::fromEnvironment();
        k = timeStep;
        radius = drum.radius;
        depth = drum.depth;
        withWall = wallLoss;
        AirGeometry geometry;
        geometry.radius = drum.radius;
        geometry.depth = drum.depth;
        geometry.wall = drum.shellWallThickness;
        geometry.delta = drum.depth / numerics.airCellsPerDepth;
        geometry.margin = geometry.delta;
        geometry.pmlCells = 0;
        air.build (geometry, drum.soundSpeed, drum.airDensity, k);
        const auto columns = air.headColumns (drum.radius);
        air.closeHeadFaces (air.planeAbove (0.0), columns, false);
        air.closeHeadFaces (air.planeAbove (-drum.depth), columns, false);
        air.finalise();
        fit = fitSquareRoot (k);
        if (withWall)
            wall.build (air, drum.radius, drum.depth, fit);
    }

    bool interior (int i, int j, int kz) const
    {
        const int cell = air.index (i, j, kz);
        return air.type (cell) != CellType::Solid && std::hypot (air.x (i), air.y (j)) < radius && air.z (kz) < 0.0
               && air.z (kz) > -depth;
    }

    // Runs `steps` steps from level 0; after each step n (level n + 1 final)
    // every thread calls perThread (thread, threads, n), then thread 0 calls
    // afterStep (n). Returns the wall's dissipated energy over the run (tracking).
    double run (long long steps, const std::function<void (int, int, long long)>& perThread,
                const std::function<void (long long)>& afterStep, int threads = 0)
    {
        ThreadPool pool (threads > 0 ? threads : std::min (8, candidateThreadCount()));
        Barrier barrier (pool.size());
        std::vector<double> lost (static_cast<std::size_t> (pool.size()), 0.0);
        pool.run ([&] (int t)
        {
            int k0, k1, i0, i1, f0 = 0, f1 = 0, c0 = 0, c1 = 0;
            splitRange (air.nz, pool.size(), t, k0, k1);
            splitRange (air.irregularCount(), pool.size(), t, i0, i1);
            if (withWall)
            {
                splitRange (wall.frictionCount(), pool.size(), t, f0, f1);
                splitRange (wall.ordinaryCellCount(), pool.size(), t, c0, c1);
            }
            for (long long n = 0; n < steps; ++n)
            {
                air.stepRegular (n, k0, k1);
                air.stepIrregular (n, i0, i1);
                barrier.wait();
                if (withWall)
                {
                    lost[static_cast<std::size_t> (t)] += wall.stepFriction (n, f0, f1);
                    lost[static_cast<std::size_t> (t)] += wall.finishCells (n, c0, c1);
                    barrier.wait();
                }
                if (perThread)
                {
                    perThread (t, pool.size(), n);
                    barrier.wait();
                }
                if (t == 0 && afterStep)
                    afterStep (n);
                barrier.wait();
            }
        });
        double total = 0.0;
        for (const double l : lost)
            total += l;
        return total;
    }
};

} // namespace taikor::candidates::fdtd::tests
