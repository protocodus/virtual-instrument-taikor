// Test 5: the SAV step keeps simulating the von Karman stretching it claims,
// on the production drum (float air, PML, lossy heads, the Don, the batter
// sized by FdtdEngineSizing.cpp). The scheme conserves heads + air + stick +
// psi^2 / 2 whatever psi does (test 1), so only a direct comparison can show
// psi drifting from sqrt(2 V(w)) (as tests/fdtd savTracking does for Berger):
//  - psi^{n-1/2} against sqrt(V(w^n) + V(w^{n-1})), V from an independent
//    von Karman potential solved to 1e-9 on the batter's displacement;
//  - psi^2 / 2 against the mean true potential, relative to the stroke's
//    energy;
//  - the time-Nyquist share of int |grad w|^2 (the Berger edge sum of
//    (w^{n+1} - 2 w^n + w^{n-1}) / 4 against that of w^n): resolved motion
//    at f holds (pi f k)^4 of it, a growing grid-scale mode ~1;
//  - the sizing: the take's own peak N* over the sizing window equals the
//    accepted sizing run's (the run is the take's first 20 ms), and the peak
//    over the whole test stays within the design tension;
//  - the drum's own SAV consistency monitor (the figure key.md reports for
//    the whole take) against the same measure from the reference.
// The heavy variant vonKarmanTrackingBergerGrid runs the same letter on the
// Berger-bound grid, which the sizing replaces, to show why it must; the
// heavy vonKarmanWholeTake holds the monitor's spurious tension over the
// take's full 4.4 s.

#include "IntegrationSupport.h"

#include "VonKarmanPotential.h"

#include <cmath>
#include <cstdlib>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::integration;

struct Tracking
{
    double worstRatio = 0.0, worstEnergy = 0.0, worstShare = 0.0, peakV = 0.0, peakRatioTime = 0.0;
    double spacing = 0.0, design = 0.0, windowPeak = 0.0, sizedPeak = 0.0, peakSav = 0.0;
    // The drum's own SAV consistency monitor (HeadState::SavConsistency, V
    // from the take's 3e-3 solves) against the same measure from the 1e-9
    // reference: largest difference of the ratio and of the spurious tension.
    double worstMonitorRatio = 0.0, worstMonitorTension = 0.0, worstReferenceTension = 0.0;
    HeadState::SavConsistency consistency;
};

Tracking track (bool sized, double seconds, Report& report)
{
    const auto& spec = nagado();
    FdtdOptions options;
    options.selfChecks = false;
    options.sizeBatterByStress = sized;
    FdtdEngine<float> engine (std::make_unique<VonKarmanPotential>(), options);
    std::string error;
    Tracking t;
    if (! engine.prepare (spec.drum, spec.mics, 48000.0, error))
    {
        report.check ("prepare", false, error);
        return t;
    }
    engine.strike (spec.stroke);
    const auto& head = engine.batterState();
    VonKarmanPotential reference;
    reference.prepare (head.grid.grid(), spec.drum.batter, spec.drum.radius, error);
    reference.solverSettings().relativeTolerance = 1e-9;
    reference.solverSettings().maxIterations = 2000;
    reference.setWarmStartOrder (1);
    auto edges = makeBergerPotential(); // its V is xi / 4 G^2 with G the edge sum
    edges->prepare (head.grid.grid(), spec.drum.batter, spec.drum.radius, error);
    const double e0 = strokeEnergy (spec.stroke);
    const double eps = FdtdNumerics::fromEnvironment().savEpsilon;
    // Numerical: every 8th internal step.
    const int every = 8;
    const int steps = static_cast<int> (seconds * engine.internalRate());
    std::vector<double> gradient (head.w[0].size()), nyquist (head.w[0].size());
    std::vector<std::array<double, 3>> samples; // {t, ratio, V}
    for (int s = 0; s < steps; s += every)
    {
        engine.advance (every);
        const long long n = engine.stepsTaken();
        const auto& now = head.w[static_cast<std::size_t> (n % 3)];
        const auto& prev = head.w[static_cast<std::size_t> ((n - 1) % 3)];
        const auto& older = head.w[static_cast<std::size_t> ((n + 1) % 3)];
        const double vOlder = reference.evaluate (older.data(), gradient.data());
        const double vNow = reference.evaluate (now.data(), gradient.data());
        const double vPrev = reference.evaluate (prev.data(), gradient.data());
        // The drum's monitor ran last at step n - 1: psi^{n-3/2} against
        // V(w^{n-1}) + V(w^{n-2}), with N* of w^{n-1}.
        if (static_cast<double> (n) / engine.internalRate() > 2.0e-4)
        {
            const double ratioRef = head.psi[n & 1] / std::sqrt (vPrev + vOlder + eps);
            const double spuriousRef = reference.peakStretchingTension() * (ratioRef - 1.0);
            t.worstMonitorRatio = std::max (t.worstMonitorRatio, std::abs (head.consistency.ratio - ratioRef));
            t.worstMonitorTension = std::max (t.worstMonitorTension, std::abs (head.consistency.spurious - spuriousRef));
            t.worstReferenceTension = std::max (t.worstReferenceTension, std::abs (spuriousRef));
        }
        for (std::size_t l = 0; l < now.size(); ++l)
            nyquist[l] = 0.25 * (now[l] - 2.0 * prev[l] + older[l]);
        const double gNyquist = edges->evaluate (nyquist.data(), gradient.data());
        const double gPrev = edges->evaluate (prev.data(), gradient.data());
        const double psi = head.psi[(n - 1) & 1];
        const double time = static_cast<double> (n) / engine.internalRate();
        samples.push_back ({ time, psi / std::sqrt (vNow + vPrev + eps), 0.5 * (vNow + vPrev) });
        t.worstEnergy = std::max (t.worstEnergy, std::abs (0.5 * (psi * psi - eps) - 0.5 * (vNow + vPrev)) / e0);
        if (gPrev > 1e3 * eps)
            t.worstShare = std::max (t.worstShare, std::sqrt (gNyquist / gPrev));
        t.peakV = std::max (t.peakV, 0.5 * (vNow + vPrev));
    }
    // The ratio where V exceeds 1e-3 of its peak and after the first 0.2 ms
    // (the explicit SAV lags the onset by about a step, tests/fdtd savTracking).
    for (const auto& sample : samples)
        if (sample[0] > 2.0e-4 && sample[2] > 1e-3 * t.peakV && std::abs (sample[1] - 1.0) > t.worstRatio)
        {
            t.worstRatio = std::abs (sample[1] - 1.0);
            t.peakRatioTime = sample[0];
        }
    t.spacing = head.grid.spacing();
    t.design = engine.designTension (0);
    t.windowPeak = head.windowSavTension;
    t.sizedPeak = engine.sizingPeak();
    t.peakSav = head.peakSavTension;
    t.consistency = head.consistency;
    report.note (format ("%s grid: h %.4f mm for %.0f N/m; worst |psi / sqrt(V^n + V^n-1) - 1| %.2e at %.3f ms; energy "
                         "%.2e E0; time-Nyquist share %.2e; peak V %.4f J; peak SAV tension %.1f N/m",
                         sized ? "sized" : "Berger-bound", 1e3 * t.spacing, t.design, t.worstRatio,
                         1e3 * t.peakRatioTime, t.worstEnergy, t.worstShare, t.peakV, t.peakSav));
    report.note (format ("drum's SAV monitor: spurious tension %.3f N/m at %.3f ms within 20 ms, %.3f N/m at %.3f ms "
                         "after; |psi^2/2 - V| %.2e J; against the 1e-9 reference (every 8th step after 0.2 ms): ratio "
                         "within %.2e, spurious tension within %.3f N/m of the reference's (largest %.3f N/m)",
                         t.consistency.worstAttack, 1e3 * t.consistency.worstAttackTime, t.consistency.worstTail,
                         1e3 * t.consistency.worstTailTime, t.consistency.worstOffset, t.worstMonitorRatio,
                         t.worstMonitorTension, t.worstReferenceTension));
    return t;
}
} // namespace

INTEGRATION_TEST (vonKarmanTracking)
{
    // Numerical: 60 ms (the contact, the attack's largest stretching and the
    // rear head's response), as tests/fdtd savTracking.
    const auto t = track (true, 0.06, report);
    report.check ("psi tracks sqrt(2 V_vK) on the sized grid", t.worstRatio < 1e-2,
                  format ("max |psi / sqrt(V(w^n) + V(w^n-1) + eps) - 1| = %.2e (after 0.2 ms, V > 1e-3 of its peak "
                          "%.4f J), limit 1e-2",
                          t.worstRatio, t.peakV));
    report.check ("SAV energy equals the true von Karman potential", t.worstEnergy < 1e-2,
                  format ("max |psi^2 / 2 - mean V(w)| = %.2e of the stroke's energy, limit 1e-2", t.worstEnergy));
    report.check ("no grid-scale growth", t.worstShare < 1e-2,
                  format ("largest time-Nyquist share of int |grad w|^2 %.2e, limit 1e-2", t.worstShare));
    report.check ("the take's first 20 ms is the accepted sizing run", t.windowPeak == t.sizedPeak && t.sizedPeak > 0.0,
                  format ("peak N* over the window: take %.6f N/m, accepted run %.6f N/m (equal bit for bit)",
                          t.windowPeak, t.sizedPeak));
    report.check ("the grid covers the stretching the take meets", t.peakSav <= t.design,
                  format ("peak SAV tension over 60 ms %.1f N/m = %.4f of the design %.1f N/m (h %.4f mm)", t.peakSav,
                          t.peakSav / t.design, t.design, 1e3 * t.spacing));
    // The monitor the take reports (FdtdEngineReport.cpp "SAV consistency")
    // uses the take's own V; it must measure what the reference measures.
    // Numerical limit: 1e-4 T (0.6 N/m), a tenth of the 1e-3 T the whole
    // take is held to (vonKarmanWholeTake).
    const double T = nagado().drum.batter.tension;
    report.check ("the drum's SAV consistency monitor agrees with the 1e-9 reference",
                  t.worstMonitorTension < 1e-4 * T,
                  format ("spurious tension N* (psi / sqrt(2 V) - 1): monitor within %.3f N/m (%.1e T) of the "
                          "reference's over 60 ms (largest %.3f N/m); ratio within %.2e; limit 1e-4 T",
                          t.worstMonitorTension, t.worstMonitorTension / T, t.worstReferenceTension, t.worstMonitorRatio));
}

INTEGRATION_HEAVY_TEST (vonKarmanWholeTake)
{
    // The SAV variable keeps the offset from sqrt(2 V) the attack leaves
    // (VonKarmanPotential.h "SAV consistency"), so its tracking has to be
    // judged over the whole take, in absolute terms: the spurious tension
    // N* (psi / sqrt(2 V) - 1) the drum's monitor measures at every step.
    // The production letter (float air, 96 kHz, sized batter) plays the Don
    // for the take's 4.4 s after the stroke (the harness renders 4.5 s with
    // the stroke at 0.1 s); TAIKOR_INTEGRATION_WHOLE_TAKE_S shortens it.
    // Numerical limit: after the attack's 20 ms, spurious tension < 1e-3 T,
    // which moves any partial by less than 0.05 % (sqrt(1 + 1e-3)), under a
    // cent; within the attack the 60 ms test holds psi / sqrt(2 V) to 1e-2.
    const auto& spec = nagado();
    const double seconds = std::getenv ("TAIKOR_INTEGRATION_WHOLE_TAKE_S")
                               ? std::atof (std::getenv ("TAIKOR_INTEGRATION_WHOLE_TAKE_S")) : 4.4;
    FdtdOptions options;
    options.selfChecks = false;
    FdtdEngine<float> engine (std::make_unique<VonKarmanPotential>(), options);
    std::string error;
    if (! engine.prepare (spec.drum, spec.mics, 48000.0, error))
        return report.check ("prepare", false, error);
    engine.strike (spec.stroke);
    const int steps = static_cast<int> (std::lround (seconds * engine.internalRate()));
    // Numerical: 0.1 s chunks, each noting the monitor's state.
    const int chunk = static_cast<int> (std::lround (0.1 * engine.internalRate()));
    const double T = spec.drum.batter.tension;
    std::string trace;
    for (int done = 0; done < steps; done += chunk)
    {
        engine.advance (std::min (chunk, steps - done));
        const auto& m = engine.batterState().consistency;
        const double t = static_cast<double> (engine.stepsTaken()) / engine.internalRate();
        if (std::fmod (t + 1e-9, 0.5) < 0.1 - 1e-9 || done + chunk >= steps)
            trace += format ("%.1f s: %+.3f N/m (ratio - 1 %+.2e); ", t, m.spurious, m.ratio - 1.0);
    }
    const auto& m = engine.batterState().consistency;
    const double e0 = strokeEnergy (spec.stroke);
    report.note ("spurious tension at the end of each 0.5 s: " + trace);
    report.check ("after the attack the SAV step misstates the stretching by < 1e-3 T over the whole take",
                  m.worstTail < 1e-3 * T && std::isfinite (m.worstTail),
                  format ("largest |N* (psi / sqrt(2 V) - 1)| after 20 ms %.3f N/m (%.2e T) at %.3f s, over %.2f s; "
                          "within the attack %.2f N/m (%.2e T) at %.3f ms; largest |psi^2 / 2 - V| %.2e J (%.1e E0)",
                          m.worstTail, m.worstTail / T, m.worstTailTime, seconds, m.worstAttack, m.worstAttack / T,
                          1e3 * m.worstAttackTime, m.worstOffset, m.worstOffset / e0));
}

INTEGRATION_HEAVY_TEST (vonKarmanTrackingBergerGrid)
{
    // The letter on the grid Berger's energy bound gives (sizing off): the
    // stress under the stick exceeds what that grid admits and the tracking
    // fails a few ms in. This documents why the sizing exists; the checks
    // assert the failure the sizing removes.
    const auto t = track (false, 0.02, report);
    report.check ("on the Berger-bound grid the SAV step loses the von Karman potential", t.worstRatio > 1e-2 || t.worstShare > 1e-2,
                  format ("max |psi / sqrt(2 V) - 1| = %.2e at %.3f ms, time-Nyquist share %.2e; peak SAV tension %.0f "
                          "N/m against the %.0f N/m the grid admits",
                          t.worstRatio, 1e3 * t.peakRatioTime, t.worstShare, t.peakSav, t.design));
}
