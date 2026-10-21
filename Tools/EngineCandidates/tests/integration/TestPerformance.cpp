// Test 7: cost. The von Karman letter adds one in-plane solve per internal
// step to the drum. It runs at the head of each step's first phase
// (FdtdEngine.h: P1, a dynamic queue), shared by a team of the pool's threads
// (FdtdOptions::potentialTeam, automaticTeamSize), while the other threads
// take the air's leapfrog, the PML and the heads' explicit parts; the letter
// costs more than Berger only by what the solve does not hide behind them.
// Measured here, back to back on the same machine (so all see its load):
//  - wall-clock per internal step of the letter (automatic team) and of
//    Berger over the first 0.1 s after the stroke (from the stroke: the
//    attack's solves are the costliest), with the letter's per-phase profile;
//  - the letter over the attack's first 30 ms with teams of 1 to 4 members
//    (the measurement behind automaticTeamSize);
//  - the von Karman evaluate alone, serial, in thread CPU time (independent of
//    the machine's load), on the take's own batter displacements after 0.1 s;
//  - the projected wall time of the 4.5 s take.
// The check is a loose regression bound (the letter within twice Berger's
// cost per step); the numbers are the report.

#include "IntegrationSupport.h"

#include "VonKarmanPotential.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::integration;

double cpuSeconds()
{
    timespec t {};
    clock_gettime (CLOCK_THREAD_CPUTIME_ID, &t);
    return static_cast<double> (t.tv_sec) + 1e-9 * static_cast<double> (t.tv_nsec);
}

struct Timing
{
    double msPerStep = 0.0, prepare = 0.0, strike = 0.0;
    std::array<double, 10> profile {};
    int threads = 0, team = 1;
    double design = 0.0, batterInclusion = 0.0, rearInclusion = 0.0;
    MembraneGrid batterGrid;
};

template <typename Make>
Timing time (Make make, const FdtdOptions& options, int steps, std::vector<std::vector<double>>* record, int recordSteps)
{
    Timing t;
    setenv ("TAIKOR_FDTD_PROFILE", "1", 1);
    FdtdEngine<float> engine (make(), options);
    std::string error;
    const auto t0 = std::chrono::steady_clock::now();
    engine.prepare (nagado().drum, nagado().mics, 48000.0, error);
    const auto t1 = std::chrono::steady_clock::now();
    engine.strike (nagado().stroke);
    const auto t2 = std::chrono::steady_clock::now();
    engine.advance (steps);
    const auto t3 = std::chrono::steady_clock::now();
    unsetenv ("TAIKOR_FDTD_PROFILE");
    t.prepare = std::chrono::duration<double> (t1 - t0).count();
    t.strike = std::chrono::duration<double> (t2 - t1).count();
    t.msPerStep = 1e3 * std::chrono::duration<double> (t3 - t2).count() / steps;
    for (const auto& p : engine.profileSeconds())
        for (std::size_t k = 0; k < 10; ++k)
            t.profile[k] += p[k] / steps;
    t.threads = static_cast<int> (engine.profileSeconds().size());
    t.team = engine.batterState().teamSize;
    t.design = engine.designTension (0);
    t.batterInclusion = engine.batterState().calibration.inclusionRadius;
    t.rearInclusion = engine.rearState().calibration.inclusionRadius;
    t.batterGrid = engine.batterState().grid.grid();
    if (record != nullptr)
        for (int s = 0; s < recordSteps; ++s)
        {
            engine.advance (1);
            record->push_back (engine.batterState().w[static_cast<std::size_t> (engine.stepsTaken() % 3)]);
        }
    return t;
}
} // namespace

INTEGRATION_TEST (performance)
{
    const int steps = 9600; // 0.1 s at the internal rate
    FdtdOptions options;
    options.selfChecks = false;
    std::vector<std::vector<double>> sequence;
    const auto makeVk = [] { return std::make_unique<VonKarmanPotential>(); };
    const auto vk = time (makeVk, options, steps, &sequence, 400);
    const auto berger = time ([] { return makeBergerPotential(); }, options, steps, nullptr, 0);
    const char* names[10] = { "potentials", "air planes", "head explicit", "air regular", "air irregular", "PML",
                              "P2a", "P2b", "P3", "barrier wait" };
    std::string profile;
    for (int k = 0; k < 10; ++k)
        profile += format ("%s %.3f; ", names[k], 1e3 * vk.profile[static_cast<std::size_t> (k)]);
    double air = 0.0;
    for (int k = 1; k <= 5; ++k)
        air += vk.profile[static_cast<std::size_t> (k)];
    // The attack with teams of 1-4 on the sized grid (the sizing reused).
    std::string teams;
    {
        FdtdOptions fixed = options;
        fixed.batterDesignTension = vk.design;
        fixed.batterInclusion = vk.batterInclusion;
        fixed.rearInclusion = vk.rearInclusion;
        for (int size = 1; size <= vk::maxTeam; ++size)
        {
            fixed.potentialTeam = size;
            const auto t = time (makeVk, fixed, 2880, nullptr, 0);
            teams += format ("%d: %.3f; ", t.team, t.msPerStep);
        }
    }
    // The evaluate alone, thread CPU time, on the take's displacements from
    // 0.1 s on (the ring-down that fills the take).
    VonKarmanPotential potential;
    std::string error;
    double cpuPerCall = 0.0;
    std::string solver;
    {
        potential.prepare (vk.batterGrid, nagado().drum.batter, nagado().drum.radius, error);
        std::vector<double> gradient (sequence.front().size());
        const double c0 = cpuSeconds();
        for (const auto& w : sequence)
            (void) potential.evaluate (w.data(), gradient.data());
        cpuPerCall = (cpuSeconds() - c0) / static_cast<double> (sequence.size());
        solver = potential.diagnostics();
    }
    const double takeSteps = 4.5 * 96000.0;
    const double projected = takeSteps * vk.msPerStep * 1e-3 + vk.prepare + vk.strike;
    report.note (format ("von Karman: %.3f ms per step (%d threads, team of %d), prepare %.1f s, strike (head rebuild "
                         "and grid sizing) %.1f s; per step, thread-ms summed over threads: %s",
                         vk.msPerStep, vk.threads, vk.team, vk.prepare, vk.strike, profile.c_str()));
    report.note (format ("Berger: %.3f ms per step, prepare %.1f s, strike %.1f s", berger.msPerStep, berger.prepare,
                         berger.strike));
    report.note (format ("von Karman over the first 30 ms, ms per step by team size: %s", teams.c_str()));
    report.note (format ("von Karman evaluate alone (serial): %.0f us thread CPU per call on the take's batter from 0.1 s "
                         "(%s); the air's P1 work per step is %.2f thread-ms",
                         1e6 * cpuPerCall, solver.c_str(), 1e3 * air));
    report.note (format ("projected wall time of the 4.5 s take at this load (self-checks excluded): %.1f min", projected / 60.0));
    report.check ("the von Karman letter costs at most twice Berger per step (regression bound)",
                  vk.msPerStep < 2.0 * berger.msPerStep,
                  format ("%.3f against %.3f ms per step (x %.2f), limit x 2", vk.msPerStep, berger.msPerStep,
                          vk.msPerStep / berger.msPerStep));
}
