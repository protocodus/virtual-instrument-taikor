// Sizing the batter's grid for a stretching nonlinearity whose tension the
// energy does not bound (von Karman, StretchingProbe.h).
//
// The explicit SAV step applies the stretching force like an extra tension
// N* on top of T, and the head's grid is stable only if its spacing is the
// stability bound at T + N* (FdtdDrumConfig.h headSpacing). For Berger N* is
// the uniform tension, bounded a priori by the drum's energy, and the grid is
// built for that bound before the stroke. For von Karman the stress
// concentrates under the stick; the energy bounds only its mean square (a
// pointwise bound from it, sqrt(2 E C) / h, is ~60 T here), so N* can only be
// measured, by running the stroke. With the batter on the Berger-bound grid
// (4.833 mm) the von Karman Don reaches N* ~ 1.4 T at ~1 ms, far past the
// 2558 N/m that grid admits, and ~2 ms in its top modes grow and psi drops
// away from sqrt(2 V) (tests/integration vonKarmanTracking reproduces this).
//
// The sizing, on the drum at rest before the stroke:
//  1. a first estimate of the peak N* over the attack (sizingWindow) from a
//     run of the stroke on the present grids at half the step: halving k
//     lets that grid carry about 50 T of extra tension (HeadGrid.h bound
//     solved for c at k/2, the bending term dominating), so this run
//     resolves the stretching whatever it reaches;
//  2. a fixed point at the production step: build the batter for the design
//     tension (stability bound at T + design), run the stroke with exactly
//     the take's configuration (same Real type, threads, air, rear head),
//     and accept when the run's own peak N* is within the design; otherwise
//     the peak becomes the design and the run repeats.
// The accepted run is the take's first sizingWindow bit for bit (the take
// then reports its own peak against the design over its whole length), so
// the grid is proven to cover the stretching it meets rather than assumed
// to. No constant enters: the design tension is a measurement of the stroke.
//
// How conservative that is. N* is a pointwise bound: it spreads the peak
// stress under the stick over the whole head as a uniform extra tension.
// At the Don's peak state (1.30 ms) a Lanczos run on the operator the SAV
// step actually applies, T K_1 + D B + the full Hessian of V, finds its
// largest eigenvalue at 0.942 of the sized-for (T + design) K_1 + D B
// (VonKarmanPotential.h, peakStretchingTension): the linear part is 0.910,
// and a uniform extra tension raises the eigenvalue by 8.0 per N/m, so the
// stress field needs the equivalent of about 3,300 N/m against the 9,329.5
// N/m design. A Lanczos-based design would give a finer grid, closer to the
// Berger letter's, but would hold only for the states it was measured on;
// the pointwise design bounds the stress stiffness G at every step whose own
// N* stays below it, which the take measures at every step and reports (the
// in-plane material part of the Hessian, which N* does not bound, adds 0.2 %
// of the eigenvalue at the peak state). The price is a batter 1.8 %
// coarser than the Berger letter's (4.9187 against 4.8330 mm on the Don),
// the one difference between the two letters besides the potential;
// tests/integration gridParity measures what it alone does to the take:
// Berger on the two spacings differs by -44.6 to -45.7 dB (difference
// energy below 4.9 kHz, windows over the first 0.1 s; peak force 120.24
// against 118.94 N, impulse within 0.02 %), against -8.4 to -11.1 dB between
// the two potentials on one spacing.

#include "FdtdEngine.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <sstream>

namespace taikor::candidates::fdtd
{
namespace
{
// Numerical cap on the fixed-point runs; each moves the design only by the
// change of the peak with the spacing (the Don accepts the first run:
// tests/integration vonKarmanTracking).
constexpr int maxSizingRuns = 4;
} // namespace

template <typename Real>
bool FdtdEngine<Real>::needsBatterSizing() const
{
    const auto& head = heads[0];
    return options.sizeBatterByStress && options.batterDesignTension <= 0.0 && batterDesign <= 0.0 && head.nonlinear
           && head.potential && head.probe != nullptr && ! head.probe->boundedByEnergy();
}

template <typename Real>
void FdtdEngine<Real>::sizeBatter (const StrokeSpec& stroke)
{
    const auto t0 = std::chrono::steady_clock::now();
    std::ostringstream log;
    char line[512];
    FdtdOptions twin = options;
    twin.selfChecks = false;
    twin.sizeBatterByStress = false;
    twin.strokeEnergy = headEnergy;
    twin.threads = threadCount;
    twin.rearSpacing = heads[1].grid.spacing();
    twin.rearInclusion = heads[1].calibration.inclusionRadius;

    // One run of the stroke over sizingWindow; the twin borrows the batter's
    // own potential (its prepare() re-builds it for the twin's grid).
    struct Result
    {
        bool ok = false;
        double peak = 0.0, spacing = 0.0, when = 0.0;
        InclusionCalibration calibration;
    };
    const auto run = [&] (const FdtdOptions& o, double sampleRate)
    {
        Result result;
        FdtdEngine<Real> engine (std::move (heads[0].potential), o);
        std::string error;
        if (engine.prepare (drumSpec, micSpec, sampleRate, error))
        {
            engine.strike (stroke);
            engine.advance (static_cast<int> (std::ceil (sizingWindow * engine.internalRate())));
            const auto& batter = engine.batterState();
            result.peak = batter.windowSavTension;
            result.when = batter.peakSavTensionTime;
            result.spacing = batter.grid.spacing();
            result.calibration = batter.calibration;
            result.ok = std::isfinite (result.peak) && result.peak > 0.0;
        }
        heads[0].potential = engine.releaseBatterPotential();
        return result;
    };

    const double T = drumSpec.batter.tension;
    FdtdOptions fine = twin;
    fine.batterSpacing = heads[0].grid.spacing();
    fine.batterInclusion = heads[0].calibration.inclusionRadius;
    const Result estimate = run (fine, 2.0 * outputRate);
    bool accepted = false;
    Result last;
    double design = estimate.peak;
    int runs = 0;
    if (estimate.ok)
    {
        std::snprintf (line, sizeof (line), "estimate at half the step on the %.4f mm grid: peak N* %.1f N/m (%.4f T) at %.3f ms",
                       1e3 * estimate.spacing, estimate.peak, estimate.peak / T, 1e3 * estimate.when);
        log << line;
        for (; runs < maxSizingRuns && ! accepted; ++runs)
        {
            FdtdOptions production = twin;
            production.batterDesignTension = design;
            last = run (production, outputRate);
            if (! last.ok)
                break;
            std::snprintf (line, sizeof (line), "; run %d: design %.1f N/m -> h %.4f mm, peak N* %.1f N/m (%.4f of the design) at %.3f ms",
                           runs + 1, design, 1e3 * last.spacing, last.peak, last.peak / design, 1e3 * last.when);
            log << line;
            if (last.peak <= design)
                accepted = true;
            else
                design = last.peak;
        }
    }
    if (! accepted)
    {
        // Keep whatever grid the batter has and say so; the take reports its
        // own peak against the bound it was built for.
        log << (estimate.ok && last.ok ? "; NOT CONVERGED: the design is the last run's peak"
                                       : "; FAILED: a sizing run did not finish, the batter keeps the Berger bound");
        if (! (estimate.ok && last.ok))
        {
            sizingReport = log.str();
            sizingSeconds = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
            std::string reprepare;
            if (! heads[0].potential->prepare (heads[0].grid.grid(), drumSpec.batter, drumSpec.radius, reprepare))
                fail ("batter potential re-prepared after failed sizing", reprepare);
            return;
        }
    }
    batterDesign = design;
    sizingAcceptedPeak = last.peak;
    batterKnownInclusion = last.calibration.inclusionRadius;
    std::string error;
    buildHeads (drumSpec, error, true);
    if (! error.empty())
        fail ("batter rebuild at the sized spacing", error);
    heads[0].calibration = last.calibration;
    buildHeadCoupling (drumSpec);
    heads[0].scratch.assign (static_cast<std::size_t> (pool->size()), std::vector<double> (heads[0].w[0].size(), 0.0));
    sizingReport = log.str();
    sizingSeconds = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
}

template bool FdtdEngine<float>::needsBatterSizing() const;
template bool FdtdEngine<double>::needsBatterSizing() const;
template void FdtdEngine<float>::sizeBatter (const StrokeSpec&);
template void FdtdEngine<double>::sizeBatter (const StrokeSpec&);

} // namespace taikor::candidates::fdtd
