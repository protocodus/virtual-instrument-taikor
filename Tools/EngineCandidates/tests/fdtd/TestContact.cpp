// Test 4: the collision alone. The stroke's tip mass meets a rigid barrier
// through the same discrete-gradient Hertz step the drum uses; its contact
// time is compared with Hertz's result for alpha = 3/2 (K. L. Johnson,
// Contact Mechanics, CUP 1985, sec. 11.4: t_c = 2.94 (5 M / 4 K)^(2/5)
// v^(-1/5), peak compression (5 M v^2 / 4 K)^(2/5)), and the rebound must
// return the impact speed (the scheme conserves energy exactly).

#include "TestSupport.h"

#include "Collision.h"
#include "FdtdDrumConfig.h"
#include "FdtdEngine.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::fdtd::tests;
} // namespace

FDTD_TEST (rigidBarrierContact)
{
    const auto& stroke = nagado().stroke;
    const auto numerics = FdtdNumerics::fromEnvironment();
    const double M = stroke.tipMass, K = stroke.hertzStiffness, v = stroke.impactSpeed;
    HertzContact hertz;
    hertz.stiffness = K;
    hertz.exponent = stroke.hertzExponent;
    const double analytic = 2.94 * std::pow (5.0 * M / (4.0 * K), 0.4) * std::pow (v, -0.2);
    const double deltaMax = std::pow (5.0 * M * v * v / (4.0 * K), 0.4);
    for (const int oversample : { numerics.oversample, 8 * numerics.oversample })
    {
        const double k = 1.0 / (48000.0 * oversample);
        // eta = -y (the barrier at y = 0), stick moving down at v, touching
        // at t = 0 exactly; the release is the interpolated zero of eta.
        double yNow = 0.0, yPrev = v * k, peakEta = 0.0, duration = 0.0;
        for (int n = 0; n < 1000000; ++n)
        {
            const double etaPrev = -yPrev, etaFree = -(2.0 * yNow - yPrev);
            double f = 0.0;
            (void) hertz.solve (etaPrev, etaFree, k * k / M, f);
            const double yNext = 2.0 * yNow - yPrev + k * k * f / M;
            peakEta = std::max (peakEta, -yNext);
            if (n > 0 && yNow <= 0.0 && yNext > 0.0)
            {
                duration = (n + (-yNow) / (yNext - yNow)) * k;
                yPrev = yNow;
                yNow = yNext;
                break;
            }
            yPrev = yNow;
            yNow = yNext;
        }
        // Two free steps later both levels are out of contact, and the
        // discrete energy is all kinetic again.
        for (int extra = 0; extra < 2; ++extra)
        {
            double f = 0.0;
            (void) hertz.solve (-yPrev, -(2.0 * yNow - yPrev), k * k / M, f);
            const double yNext = 2.0 * yNow - yPrev + k * k * f / M;
            yPrev = yNow;
            yNow = yNext;
        }
        const double exit = (yNow - yPrev) / k;
        report.check (format ("contact time at %.0f kHz", 48.0 * oversample), std::abs (duration / analytic - 1.0) < 0.03,
                      format ("%.4f ms vs Hertz %.4f ms (%+.2f %%; %.0f steps), limit 3 %%",
                              1e3 * duration, 1e3 * analytic, 100.0 * (duration / analytic - 1.0), duration / k));
        report.check (format ("peak compression at %.0f kHz", 48.0 * oversample), std::abs (peakEta / deltaMax - 1.0) < 0.01,
                      format ("%.4f um vs %.4f um, limit 1 %%", 1e6 * peakEta, 1e6 * deltaMax));
        report.check (format ("rebound speed at %.0f kHz", 48.0 * oversample), std::abs (exit / v - 1.0) < 1e-9,
                      format ("%.12f m/s vs impact %.12f m/s (energy-conserving collision)", exit, v));
    }
}

// Development (TAIKOR_FDTD_SCAN set): the production stroke's contact force
// history on the production drum, written to TAIKOR_FDTD_DUMP (float64 pairs
// time, force) and summarised per episode.
FDTD_TEST (strokeForceHistory)
{
    if (! std::getenv ("TAIKOR_FDTD_SCAN"))
        return;
    FdtdOptions options;
    options.selfChecks = false;
    FdtdEngine<float> engine (makeBergerPotential(), options);
    std::string error;
    engine.prepare (nagado().drum, nagado().mics, 48000.0, error);
    engine.strike (nagado().stroke);
    const int steps = static_cast<int> (0.012 * engine.internalRate());
    std::vector<double> dump;
    bool in = false;
    double start = 0.0, impulse = 0.0, peak = 0.0;
    for (int n = 0; n < steps; ++n)
    {
        engine.advance (1);
        const double t = engine.stepsTaken() / engine.internalRate(), f = engine.lastForce();
        dump.push_back (t);
        dump.push_back (f);
        if (f > 0.0 && ! in)
        {
            in = true;
            start = t;
            impulse = peak = 0.0;
        }
        if (in)
        {
            impulse += f / engine.internalRate();
            peak = std::max (peak, f);
        }
        if (f <= 0.0 && in)
        {
            in = false;
            report.note (format ("episode %.4f - %.4f ms: peak %.2f N, impulse %.5f N s", 1e3 * start, 1e3 * t, peak, impulse));
        }
    }
    if (const char* path = std::getenv ("TAIKOR_FDTD_DUMP"))
        if (FILE* file = std::fopen (path, "wb"))
        {
            std::fwrite (dump.data(), sizeof (double), dump.size(), file);
            std::fclose (file);
        }
}
