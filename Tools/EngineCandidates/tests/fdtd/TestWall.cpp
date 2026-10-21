// The barrel wall's Kirchhoff boundary layer (WallLoss.h): the sqrt(s) fit,
// the staircase's areas, and passivity of the sealed barrel with the loss.

#include "WallSupport.h"

#include <algorithm>
#include <cmath>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::fdtd::tests;
constexpr double pi = 3.14159265358979323846;
} // namespace

// The fit's accuracy in the discrete domain, per band (WallLossFit.cpp).
FDTD_TEST (wallSquareRootFit)
{
    const auto numerics = FdtdNumerics::fromEnvironment();
    const double k = 1.0 / (48000.0 * numerics.oversample);
    const auto fit = fitSquareRoot (k);
    const double bands[4][2] = { { 20.0, 1000.0 }, { 1000.0, 9000.0 }, { 9000.0, 12000.0 }, { 12000.0, 16000.0 } };
    double worst[2][4] = {}, reactive[2][4] = {};
    for (int kind = 0; kind < 2; ++kind)
        for (int b = 0; b < 4; ++b)
            for (int s = 0; s <= 400; ++s)
            {
                const double f = bands[b][0] * std::pow (bands[b][1] / bands[b][0], s / 400.0);
                const double w = 2.0 * pi * f;
                auto target = std::sqrt (std::complex<double> (0.0, w));
                if (kind == 1)
                    target *= std::cos (0.5 * w * k) * std::cos (0.5 * w * k);
                const auto got = realisedSquareRoot (fit, kind == 0 ? fit.thermal : fit.viscous, k, w);
                worst[kind][b] = std::max (worst[kind][b], std::abs (got.real() / target.real() - 1.0));
                reactive[kind][b] = std::max (reactive[kind][b], std::abs (got.imag() / target.imag() - 1.0));
            }
    bool positive = true;
    std::string weights;
    for (std::size_t i = 0; i < fit.xi.size(); ++i)
    {
        positive = positive && fit.thermal[i] >= 0.0 && fit.viscous[i] >= 0.0;
        weights += format ("%.3g Hz: %.4g/%.4g; ", fit.xi[i] / (2.0 * pi), fit.thermal[i], fit.viscous[i]);
    }
    report.note (format ("rate %.0f Hz, %zu poles (thermal/viscous weights): %s", 1.0 / k, fit.xi.size(), weights.c_str()));
    for (int kind = 0; kind < 2; ++kind)
        report.note (format ("%s real-part error 20-1k %.3f %%, 1-9k %.3f %%, 9-12k %.2f %%, 12-16k %.1f %% (reactance %.0f/%.0f/%.0f/%.0f %%)",
                             kind == 0 ? "thermal" : "viscous", 100 * worst[kind][0], 100 * worst[kind][1],
                             100 * worst[kind][2], 100 * worst[kind][3], 100 * reactive[kind][0], 100 * reactive[kind][1],
                             100 * reactive[kind][2], 100 * reactive[kind][3]));
    report.check ("wall sqrt(s) weights non-negative (passive branches)", positive, "every mu_i >= 0");
    report.check ("wall sqrt(s) loss within 1.5 % over 20 Hz - 9 kHz",
                  std::max ({ worst[0][0], worst[0][1], worst[1][0], worst[1][1] }) < 0.015,
                  format ("worst thermal %.3f %%, viscous %.3f %%", 100 * std::max (worst[0][0], worst[0][1]),
                          100 * std::max (worst[1][0], worst[1][1])));
}

// The staircase's wall and friction areas against the barrel's side wall.
FDTD_TEST (wallGeometry)
{
    const auto& drum = nagado().drum;
    const auto numerics = FdtdNumerics::fromEnvironment();
    SealedBarrel barrel;
    barrel.build (drum, 1.0 / (48000.0 * numerics.oversample), true);
    const double side = 2.0 * pi * drum.radius * drum.depth;
    const auto& wall = barrel.wall;
    // Each wall face gives half its area to each of its cell's two z faces;
    // at the rigid ends one of them is the end, so the axial area falls
    // short by half a plane's worth at each end: 1 / N_L of the wall.
    const double axialExpected = side * (1.0 - 1.0 / numerics.airCellsPerDepth);
    report.note (format ("%d wall cells, %d friction elements on %d faces, %d cells corrected; wall %.5f m^2, axial %.5f, azimuthal %.5f (side 2 pi a L %.5f)",
                         wall.wallCellCount(), wall.frictionCount(), wall.frictionFaceCount(), wall.ordinaryCellCount(), wall.wallArea(),
                         wall.axialArea(), wall.azimuthalArea(), side));
    report.check ("staircase wall faces carry the barrel's side area", std::abs (wall.rawWallArea() / side - 1.0) < 0.015
                  && std::abs (wall.wallArea() / side - 1.0) < 1e-12,
                  format ("%.4f of 2 pi a L before the normalisation (limit 1.5 %%), %.12f after", wall.rawWallArea() / side,
                          wall.wallArea() / side));
    report.check ("axial friction area", std::abs (wall.axialArea() / axialExpected - 1.0) < 0.01,
                  format ("%.4f of 2 pi a L (1 - 1/N_L), limit 1 %%", wall.axialArea() / axialExpected));
    report.check ("azimuthal friction area", std::abs (wall.azimuthalArea() / side - 1.0) < 1e-12,
                  format ("%.12f of 2 pi a L (one element per wall cell)", wall.azimuthalArea() / side));
}

// Passivity: from a smooth blob (every mode family) the energy of air +
// boundary storage never increases, it plus the dissipated energy stays
// constant, and at zero strength the wall leaves the air bit-identical.
FDTD_TEST (wallPassivity)
{
    const auto& drum = nagado().drum;
    const auto numerics = FdtdNumerics::fromEnvironment();
    const double k = 1.0 / (48000.0 * numerics.oversample);
    const auto seed = [&drum] (SealedBarrel& b)
    {
        const std::array<double, 3> centre { 0.31, 0.12, -0.21 }; // near the wall
        const double s = 1.5 * b.air.delta;
        for (int kz = 0; kz < b.air.nz; ++kz)
            for (int j = 0; j < b.air.ny; ++j)
                for (int i = 0; i < b.air.nx; ++i)
                    if (b.interior (i, j, kz))
                    {
                        const double dx = b.air.x (i) - centre[0], dy = b.air.y (j) - centre[1], dz = b.air.z (kz) - centre[2];
                        const int cell = b.air.index (i, j, kz);
                        b.air.level (0)[cell] = std::exp (-(dx * dx + dy * dy + dz * dz) / (2.0 * s * s));
                        b.air.level (-1)[cell] = 0.0; // moving: the boundary layer works from the first step
                    }
        (void) drum;
    };
    const long long steps = 3000;
    // Lossy run with the energy each step.
    SealedBarrel lossy;
    lossy.build (drum, k, true);
    lossy.wall.setTracking (true);
    seed (lossy);
    std::vector<double> energy;
    std::vector<double> lostSoFar;
    double lostTotal = 0.0;
    {
        ThreadPool pool (std::min (8, candidateThreadCount()));
        Barrier barrier (pool.size());
        std::vector<double> lost (static_cast<std::size_t> (pool.size()), 0.0);
        auto& air = lossy.air;
        auto& wall = lossy.wall;
        pool.run ([&] (int t)
        {
            int k0, k1, i0, i1, f0, f1, c0, c1;
            splitRange (air.nz, pool.size(), t, k0, k1);
            splitRange (air.irregularCount(), pool.size(), t, i0, i1);
            splitRange (wall.frictionCount(), pool.size(), t, f0, f1);
            splitRange (wall.ordinaryCellCount(), pool.size(), t, c0, c1);
            for (long long n = 0; n < steps; ++n)
            {
                air.stepRegular (n, k0, k1);
                air.stepIrregular (n, i0, i1);
                barrier.wait();
                lost[static_cast<std::size_t> (t)] = wall.stepFriction (n, f0, f1) + wall.finishCells (n, c0, c1);
                barrier.wait();
                if (t == 0)
                {
                    for (const double l : lost)
                        lostTotal += l;
                    energy.push_back (air.energy (n) + wall.energy (n));
                    lostSoFar.push_back (lostTotal);
                }
                barrier.wait();
            }
        });
    }
    double worstRise = 0.0, worstBalance = 0.0;
    const double e0 = energy.front();
    for (std::size_t n = 1; n < energy.size(); ++n)
    {
        worstRise = std::max (worstRise, (energy[n] - energy[n - 1]) / e0);
        worstBalance = std::max (worstBalance, std::abs (energy[n] + lostSoFar[n] - energy[0] - lostSoFar[0]) / e0);
    }
    report.note (format ("sealed barrel %lld steps: energy %.6e -> %.6e J, wall dissipated %.6e J (%.3e of E0)",
                         steps, e0, energy.back(), lostTotal, lostTotal / e0));
    report.check ("wall loss: energy never rises", worstRise <= 1.0e-13,
                  format ("largest step rise %.2e of E0, limit 1e-13 (rounding)", worstRise));
    report.check ("wall loss: energy + dissipated constant", worstBalance < 1.0e-11,
                  format ("worst drift %.2e of E0, limit 1e-11", worstBalance));
    // Zero strength: bit-identical to the air without any wall.
    SealedBarrel rigid, zero;
    rigid.build (drum, k, false);
    zero.build (drum, k, true);
    zero.wall.scale (0.0, 0.0);
    seed (rigid);
    seed (zero);
    rigid.run (600, nullptr, nullptr);
    zero.run (600, nullptr, nullptr);
    long long differing = 0;
    for (long long level = 600; level >= 599; --level)
        for (long long cell = 0; cell < rigid.air.cells(); ++cell)
            differing += rigid.air.level (level)[cell] != zero.air.level (level)[cell] ? 1 : 0;
    report.check ("wall loss at zero strength leaves the air bit-identical", differing == 0,
                  format ("%lld cells differ after 600 steps", differing));
}
