// Test 6: the batter's evaluate shared by a team of threads (VonKarman/Team.h,
// FdtdOptions::potentialTeam) computes what the serial evaluate computes.
//  1. The single-precision multigrid V-cycle of a team is the serial
//     V-cycle's arithmetic entry for entry: the same z bit for bit, for teams
//     of 2, 3 and 4 (InPlaneSolver.h).
//  2. On the take's own batter displacements (the attack, where the solves
//     are the costliest), a team's V and gradient match the serial ones to
//     rounding: the conjugate-gradient reductions sum the members' partials,
//     so they round differently, and the single-precision V-cycle turns those
//     last-bit differences of its double input into float rounding (~1e-8
//     per entry), far below the solve's own tolerance (3e-3).
//  3. The whole drum played with a team of 2 and serially: the capsules and
//     the contact agree far below the in-plane tolerance's own effect on the
//     take (toleranceInsensitivity).

#include "IntegrationSupport.h"

#include "VonKarmanPotential.h"

#include <cmath>
#include <cstring>
#include <random>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::integration;

// The batter's displacements over `count` steps from `startMs` after the
// stroke, on the von Karman drum (serial evaluate, Berger-bound batter grid:
// the recording needs shapes, not the sized grid).
std::vector<std::vector<double>> recordBatter (double startMs, int count, MembraneGrid& grid, std::string& error)
{
    FdtdOptions options;
    options.selfChecks = false;
    options.sizeBatterByStress = false;
    options.potentialTeam = 1;
    FdtdEngine<float> engine (std::make_unique<VonKarmanPotential>(), options);
    std::vector<std::vector<double>> sequence;
    if (! engine.prepare (nagado().drum, nagado().mics, 48000.0, error))
        return sequence;
    engine.strike (nagado().stroke);
    engine.advance (static_cast<int> (std::lround (startMs * 1e-3 * engine.internalRate())));
    for (int s = 0; s < count; ++s)
    {
        engine.advance (1);
        sequence.push_back (engine.batterState().w[static_cast<std::size_t> (engine.stepsTaken() % 3)]);
    }
    grid = engine.batterState().grid.grid();
    return sequence;
}
} // namespace

INTEGRATION_TEST (teamEvaluate)
{
    std::string error;
    MembraneGrid grid;
    const auto sequence = recordBatter (0.5, 300, grid, error);
    if (! error.empty() || sequence.empty())
    {
        report.check ("record the batter", false, error);
        return;
    }
    const auto& spec = nagado();
    ThreadPool pool (vk::maxTeam);

    // 1. The V-cycle, bit for bit.
    {
        vk::PlaneStress c;
        const auto& head = spec.drum.batter;
        c.c11 = head.youngsModulus * head.thickness / (1.0 - head.poissonRatio * head.poissonRatio);
        c.c12 = head.poissonRatio * c.c11;
        c.c66 = 0.5 * (1.0 - head.poissonRatio) * c.c11;
        vk::InPlaneSolver solver;
        solver.build (grid.nx, grid.interior, c, {}, error);
        const auto& L = solver.fine();
        std::vector<float> r (static_cast<std::size_t> (L.size()), 0.0f), serial (r.size()), team (r.size());
        std::mt19937 random (7);
        std::uniform_real_distribution<float> uniform (-1.0f, 1.0f);
        for (int component = 0; component < 2; ++component)
            for (int k = 0; k < L.N; ++k)
                if ((component == 0 ? L.freeU : L.freeV)[static_cast<std::size_t> (k)] != 0)
                    r[static_cast<std::size_t> (L.base (component) + k)] = uniform (random);
        solver.vcycle (0, r.data(), serial.data());
        std::string detail;
        bool identical = true;
        for (int size = 2; size <= vk::maxTeam; ++size)
        {
            Barrier barrier (size);
            vk::TeamShared shared;
            // z's pads must be zero on entry (the stencils read them), as the
            // solver's own work vectors are.
            std::fill (team.begin(), team.end(), 0.0f);
            pool.run ([&] (int thread)
            {
                if (thread >= size)
                    return;
                vk::Team member;
                member.member = thread;
                member.size = size;
                member.barrier = &barrier;
                member.shared = &shared;
                solver.vcycle (0, r.data(), team.data(), member);
            });
            const bool same = std::memcmp (serial.data(), team.data(), serial.size() * sizeof (float)) == 0;
            identical = identical && same;
            long differing = 0;
            float largest = 0.0f;
            for (std::size_t k = 0; k < serial.size(); ++k)
                if (serial[k] != team[k])
                {
                    ++differing;
                    largest = std::max (largest, std::abs (serial[k] - team[k]));
                }
            detail += format ("team of %d: %s", size, same ? "identical; " : "DIFFERS");
            if (! same)
                detail += format (" in %ld entries, by up to %.2e; ", differing, static_cast<double> (largest));
        }
        report.check ("the team V-cycle is the serial V-cycle bit for bit", identical,
                      detail + format ("%ld level-vector entries", L.size()));
    }

    // 2. The evaluate on the recorded attack, serial against teams.
    {
        const std::size_t n = sequence.front().size();
        std::vector<double> gs (n), gt (n);
        std::string detail;
        bool within = true;
        for (int size = 2; size <= vk::maxTeam; ++size)
        {
            VonKarmanPotential serial, team;
            serial.prepare (grid, spec.drum.batter, spec.drum.radius, error);
            team.prepare (grid, spec.drum.batter, spec.drum.radius, error);
            team.setTeam (size);
            double worstV = 0.0, worstG = 0.0;
            for (const auto& w : sequence)
            {
                const double vs = serial.evaluate (w.data(), gs.data());
                double vt = 0.0;
                pool.run ([&] (int thread)
                {
                    if (thread >= size)
                        return;
                    const double v = team.evaluateMember (w.data(), gt.data(), thread);
                    if (thread == 0)
                        vt = v;
                });
                worstV = std::max (worstV, std::abs (vt / vs - 1.0));
                double num = 0.0, den = 0.0;
                for (std::size_t l = 0; l < n; ++l)
                {
                    num += (gt[l] - gs[l]) * (gt[l] - gs[l]);
                    den += gs[l] * gs[l];
                }
                worstG = std::max (worstG, std::sqrt (num / std::max (1e-300, den)));
            }
            within = within && worstV < 1e-9 && worstG < 1e-6
                     && serial.stats().iterations == team.stats().iterations;
            detail += format ("team of %d: |V/V_serial - 1| %.1e, gradient %.1e, PCG iterations %ld vs %ld; ", size,
                              worstV, worstG, team.stats().iterations, serial.stats().iterations);
        }
        report.check ("a team's evaluate matches the serial one to rounding on the attack", within,
                      detail + format ("%zu steps from 0.5 ms, limits 1e-9 (V) and 1e-6 (relative L2 gradient), same "
                                       "iteration counts", sequence.size()));
    }

    // 3. The drum with a team of 2 against the serial evaluate.
    {
        const int steps = static_cast<int> (std::ceil (0.010 * 96000.0));
        FdtdOptions options;
        options.selfChecks = false;
        options.sizeBatterByStress = false;
        options.potentialTeam = 1;
        auto a = runTake<double> (std::make_unique<VonKarmanPotential>(), options, spec.stroke, steps, error);
        options.potentialTeam = 2;
        options.batterInclusion = a.engine->batterState().calibration.inclusionRadius;
        options.rearInclusion = a.engine->rearState().calibration.inclusionRadius;
        auto b = runTake<double> (std::make_unique<VonKarmanPotential>(), options, spec.stroke, steps, error);
        const double left = differenceDb (a.left, b.left, 0, static_cast<std::size_t> (steps));
        const double right = differenceDb (a.right, b.right, 0, static_cast<std::size_t> (steps));
        const double dPeak = b.engine->contact().peakForce / a.engine->contact().peakForce - 1.0;
        report.check ("the drum with a team of 2 plays the serial take", a.engine->batterState().teamSize == 1
                          && b.engine->batterState().teamSize == 2 && std::max (left, right) < -120.0 && std::abs (dPeak) < 1e-9,
                      format ("difference energy over 10 ms, left / right %.1f / %.1f dB (limit -120 dB, against the "
                              "tolerance's own effect near -60 dB), peak force %+.1e (limit 1e-9); team sizes %d / %d",
                              left, right, dPeak, a.engine->batterState().teamSize, b.engine->batterState().teamSize));
    }
}
