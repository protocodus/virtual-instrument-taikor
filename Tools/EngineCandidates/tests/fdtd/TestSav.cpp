// Test 10: the SAV step keeps simulating the nonlinearity it claims. The
// scheme conserves heads + air + stick + psi^2 / 2 exactly whatever psi does
// (test 3 checks that), so test 3 cannot see psi drifting away from
// sqrt(2 V(w)); this test compares them directly on the production drum
// (lossy heads, float air with the PML, the Don), on both heads:
//  - the tension the SAV step applies, sqrt(2 xi) psi^{n-1/2}, against the
//    Berger tension xi G = sqrt(2 xi (V(w^n) + V(w^{n-1}))), V recomputed by
//    an independent Berger potential from the heads' displacement;
//  - the energy psi^2 / 2 against the mean true potential, relative to the
//    stroke's energy;
//  - the share of int |grad w|^2 held by the time-Nyquist part of w,
//    G((w^{n+1} - 2 w^n + w^{n-1}) / 4) / G(w^n) with G the Berger edge sum:
//    resolved motion at f holds (pi f k)^4 of it, a grid-scale mode ~1.
// With the head spacing at the linear stability bound the Don failed the
// first and last from 2.3 ms on (SAV tension 15-40 N/m against a Berger
// tension of 300-600 N/m, time-Nyquist share 0.39-0.49, FdtdDrumConfig.h
// headSpacing); the stiffened bound keeps them near 1e-3 and 1e-4 and the
// energy near 1e-5 of E0. The limits (1e-2 each) sit between.

#include "TestSupport.h"

#include "FdtdEngine.h"

#include <array>
#include <cmath>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::fdtd::tests;

struct Tracking
{
    double worstTension = 0.0, worstEnergy = 0.0, worstShare = 0.0, peakTension = 0.0, worstRatio = 0.0;
};
} // namespace

FDTD_TEST (savTracking)
{
    const auto& spec = nagado();
    FdtdOptions options;
    options.selfChecks = false;
    FdtdEngine<float> engine (makeBergerPotential(), options);
    std::string error;
    if (! engine.prepare (spec.drum, spec.mics, 48000.0, error))
    {
        report.check ("prepare", false, error);
        return;
    }
    engine.strike (spec.stroke);
    const double e0 = strokeEnergy (spec.stroke);
    const double eps = FdtdNumerics::fromEnvironment().savEpsilon;
    const HeadState* heads[2] = { &engine.batterState(), &engine.rearState() };
    const HeadSpec* specs[2] = { &spec.drum.batter, &spec.drum.rear };
    std::unique_ptr<NonlinearPotential> reference[2] = { makeBergerPotential(), makeBergerPotential() };
    double xi[2];
    for (int h = 0; h < 2; ++h)
    {
        reference[h]->prepare (heads[h]->grid.grid(), *specs[h], spec.drum.radius, error);
        const double nu = specs[h]->poissonRatio;
        xi[h] = specs[h]->youngsModulus * specs[h]->thickness
                / (2.0 * (1.0 - nu * nu) * 3.14159265358979323846 * spec.drum.radius * spec.drum.radius);
    }
    // Numerical: 60 ms (the contact, the attack's largest Berger tension and
    // the rear head's response), every 8th internal step.
    const int every = 8;
    const int steps = static_cast<int> (0.06 * engine.internalRate());
    Tracking t[2];
    std::vector<std::array<double, 3>> samples[2]; // {SAV tension, Berger tension, psi / sqrt(2 V + eps)}
    std::vector<double> gradient, nyquist;
    for (int s = 0; s < steps; s += every)
    {
        engine.advance (every);
        const long long n = engine.stepsTaken();
        for (int h = 0; h < 2; ++h)
        {
            const auto& head = *heads[h];
            const auto& now = head.w[static_cast<std::size_t> (n % 3)];
            const auto& prev = head.w[static_cast<std::size_t> ((n - 1) % 3)];
            const auto& older = head.w[static_cast<std::size_t> ((n + 1) % 3)]; // w^{n-2}
            gradient.assign (now.size(), 0.0);
            nyquist.assign (now.size(), 0.0);
            const double vNow = reference[h]->evaluate (now.data(), gradient.data());
            const double vPrev = reference[h]->evaluate (prev.data(), gradient.data());
            for (std::size_t l = 0; l < now.size(); ++l)
                nyquist[l] = 0.25 * (now[l] - 2.0 * prev[l] + older[l]);
            const double vNyquist = reference[h]->evaluate (nyquist.data(), gradient.data());
            const double psi = head.psi[(n - 1) & 1];
            // The tension the SAV step applies, sqrt(2 xi) psi, and the one an
            // exactly tracked psi = sqrt(2 V) would apply, xi G.
            const double tSav = std::sqrt (2.0 * xi[h]) * psi;
            const double tTrue = std::sqrt (2.0 * xi[h] * (vNow + vPrev));
            samples[h].push_back ({ tSav, tTrue, psi / std::sqrt (vNow + vPrev + eps) });
            t[h].worstEnergy = std::max (t[h].worstEnergy, std::abs (0.5 * (psi * psi - eps) - 0.5 * (vNow + vPrev)) / e0);
            // G ratio = sqrt of the V ratio (V = xi / 4 G^2).
            if (vPrev > 1e3 * eps)
                t[h].worstShare = std::max (t[h].worstShare, std::sqrt (vNyquist / vPrev));
            t[h].peakTension = std::max (t[h].peakTension, tTrue);
        }
    }
    // The tension error relative to the head's peak Berger tension (a ratio
    // psi / sqrt(2 V) is ill-conditioned while V is still tiny: over the first
    // 0.1 ms the explicit SAV lags the onset by about a step, 0.13 against
    // 0.21 N/m at 0.04 ms and 3.75 against 3.90 N/m at 0.08 ms, 4e-4 of the
    // peak); the ratio is reported where the tension exceeds 1 % of its peak,
    // and separately after the first 0.2 ms (past that onset lag).
    double lateRatio[2] = { 0.0, 0.0 };
    for (int h = 0; h < 2; ++h)
        for (std::size_t i = 0; i < samples[h].size(); ++i)
        {
            const auto& sample = samples[h][i];
            t[h].worstTension = std::max (t[h].worstTension, std::abs (sample[0] - sample[1]) / t[h].peakTension);
            if (sample[1] > 0.01 * t[h].peakTension)
            {
                t[h].worstRatio = std::max (t[h].worstRatio, std::abs (sample[2] - 1.0));
                if (static_cast<double> ((i + 1) * every) / engine.internalRate() > 2.0e-4)
                    lateRatio[h] = std::max (lateRatio[h], std::abs (sample[2] - 1.0));
            }
        }
    const char* names[2] = { "batter", "rear" };
    for (int h = 0; h < 2; ++h)
    {
        report.check (format ("%s: the SAV tension is the Berger tension over 60 ms", names[h]), t[h].worstTension < 1e-2,
                      format ("max |sqrt(2 xi) psi - xi G(w)| = %.2e of the peak Berger tension %.1f N/m, limit 1e-2; "
                              "max |psi / sqrt(V(w^n) + V(w^n-1) + eps) - 1| = %.2e where the tension exceeds 1 %% "
                              "of its peak, %.2e after the first 0.2 ms (h %.4f mm)",
                              t[h].worstTension, t[h].peakTension, t[h].worstRatio, lateRatio[h],
                              1e3 * heads[h]->grid.spacing()));
        report.check (format ("%s: SAV energy equals the true potential", names[h]), t[h].worstEnergy < 1e-2,
                      format ("max |psi^2 / 2 - mean V(w)| = %.2e of the stroke's %.4f J, limit 1e-2",
                              t[h].worstEnergy, e0));
        report.check (format ("%s: no grid-scale growth", names[h]), t[h].worstShare < 1e-2,
                      format ("largest time-Nyquist share of int |grad w|^2 %.2e, limit 1e-2", t[h].worstShare));
    }
}
