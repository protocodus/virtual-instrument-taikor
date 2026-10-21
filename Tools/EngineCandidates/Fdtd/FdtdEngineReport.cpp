// Energy accounting, self-checks and diagnostics of the finite-difference drum.

#include "FdtdEngine.h"

#include "AirChecks.h"

#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <sstream>

namespace taikor::candidates::fdtd
{
namespace
{
inline int mod3 (long long n) noexcept { return static_cast<int> (((n % 3) + 3) % 3); }

std::string format (const char* pattern, ...)
{
    char buffer[4096];
    va_list args;
    va_start (args, pattern);
    std::vsnprintf (buffer, sizeof (buffer), pattern, args);
    va_end (args);
    return buffer;
}

// The contact's robust figures from the force history f (one value per
// internal step k): the 5-95 % span of its impulse, and the strongest
// spectral peak of the force above 4 kHz within that span (the footprint
// patch's oscillation on the Hertz spring, Decimator.h) with its amplitude
// relative to the mean force over the span. Hann-windowed DFT on a 10 Hz grid.
struct ForceFigures
{
    double t5 = 0.0, t95 = 0.0, oscillationHz = 0.0, oscillationDb = -300.0;
};

ForceFigures forceFigures (const std::vector<double>& f, double k)
{
    ForceFigures r;
    double total = 0.0;
    for (const double v : f)
        total += std::max (0.0, v);
    if (total <= 0.0)
        return r;
    double sum = 0.0;
    std::size_t i5 = 0, i95 = 0;
    for (std::size_t i = 0; i < f.size(); ++i)
    {
        const double before = sum;
        sum += std::max (0.0, f[i]);
        if (before < 0.05 * total && sum >= 0.05 * total)
            i5 = i;
        if (before < 0.95 * total && sum >= 0.95 * total)
            i95 = i;
    }
    r.t5 = (i5 + 1) * k;
    r.t95 = (i95 + 1) * k;
    const std::size_t n = i95 > i5 ? i95 - i5 + 1 : 0;
    if (n < 16)
        return r;
    double mean = 0.0;
    for (std::size_t i = i5; i <= i95; ++i)
        mean += f[i] / static_cast<double> (n);
    const double pi = 3.14159265358979323846;
    double best = 0.0;
    for (double hz = 4000.0; hz < 0.45 / k; hz += 10.0)
    {
        double re = 0.0, im = 0.0, wsum = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            const double w = 0.5 - 0.5 * std::cos (2.0 * pi * (i + 0.5) / static_cast<double> (n));
            const double phase = 2.0 * pi * hz * static_cast<double> (i) * k;
            re += w * f[i5 + i] * std::cos (phase);
            im -= w * f[i5 + i] * std::sin (phase);
            wsum += w;
        }
        const double amplitude = 2.0 * std::sqrt (re * re + im * im) / wsum;
        if (amplitude > best)
        {
            best = amplitude;
            r.oscillationHz = hz;
        }
    }
    r.oscillationDb = 20.0 * std::log10 (std::max (1e-300, best) / std::max (1e-300, mean));
    return r;
}
} // namespace

template <typename Real>
double FdtdEngine<Real>::energyHeads() const
{
    double total = 0.0;
    for (const auto& head : heads)
        total += head.grid.energy (head.w[mod3 (step)].data(), head.w[mod3 (step - 1)].data());
    return total;
}

template <typename Real>
double FdtdEngine<Real>::energyAir() const
{
    if (! options.withAir)
        return 0.0;
    return domain.energy (step - 1) + (wallActive ? wall.energy (step - 1) : 0.0);
}

template <typename Real>
double FdtdEngine<Real>::dissipatedWall() const
{
    double total = 0.0;
    for (const double lost : wallLost)
        total += lost;
    return total;
}

template <typename Real>
double FdtdEngine<Real>::energySav() const
{
    double total = 0.0;
    for (const auto& head : heads)
        if (head.nonlinear)
        {
            const double psi = head.psi[(step - 1) & 1];
            total += 0.5 * (psi * psi - numerics.savEpsilon);
        }
    return total;
}

template <typename Real>
double FdtdEngine<Real>::energyStick() const
{
    if (! struckYet)
        return 0.0;
    const auto& head = heads[0];
    const double y1 = stickY[mod3 (step)], y0 = stickY[mod3 (step - 1)];
    const double eta1 = head.grid.dot (head.footprint.data(), head.w[mod3 (step)].data()) - y1;
    const double eta0 = head.grid.dot (head.footprint.data(), head.w[mod3 (step - 1)].data()) - y0;
    const double v = (y1 - y0) / k;
    return 0.5 * stickMass * v * v + 0.5 * (hertz.potential (eta1) + hertz.potential (eta0));
}

// Self-checks on the production configuration, run once at the stroke (they
// need the stroke's footprint and contact law): the lossless energy drift of
// the same drum in a reflecting box, the absorbing layer's reflection at the
// production air grid, and the staircase heads' doublet split.
template <typename Real>
void FdtdEngine<Real>::runSelfChecks (const DrumSpec& drum, const MicSpec& mics, double sampleRate)
{
    std::ostringstream out;
    const auto t0 = std::chrono::steady_clock::now();
    {
        FdtdOptions twin = options;
        twin.lossless = true;
        twin.rigidBox = true;
        twin.selfChecks = false;
        twin.batterInclusion = heads[0].calibration.inclusionRadius;
        twin.rearInclusion = heads[1].calibration.inclusionRadius;
        twin.strokeEnergy = headEnergy;
        // The batter's grid as the take has it (its sized design, if any).
        twin.sizeBatterByStress = false;
        twin.batterDesignTension = batterDesign;
        const std::string batterName = heads[0].potential ? heads[0].potential->name() : std::string ("no");
        // The check borrows the batter's own potential (the twin's batter grid
        // is identical), so the letter's nonlinearity is the one checked.
        FdtdEngine<double> check (std::move (heads[0].potential), twin);
        std::string error;
        if (check.prepare (drum, mics, sampleRate, error))
        {
            check.strike (strokeSpec);
            const auto total = [&check]
            { return check.energyHeads() + check.energyAir() + check.energyStick() + check.energySav(); };
            check.advance (1);
            const double e0 = total();
            double worst = 0.0;
            // Numerical: 2 ms of internal steps, through the whole first contact
            // pulse's rise; each block is checked.
            const int blocks = 24, per = static_cast<int> (std::ceil (2.0e-3 * check.internalRate() / blocks));
            for (int b = 0; b < blocks; ++b)
            {
                check.advance (per);
                worst = std::max (worst, std::abs (total() / e0 - 1.0));
            }
            out << format ("- Lossless energy check (same grids, reflecting box, %s batter, Berger rear and collision "
                           "on, double air): relative drift %.2e over %lld steps (%.1f ms); stick energy at the end %.3f "
                           "of the start.\n",
                           batterName.c_str(), worst, check.stepsTaken(), 1e3 * check.stepsTaken() / check.internalRate(),
                           check.energyStick() / e0);
        }
        heads[0].potential = check.releaseBatterPotential();
        std::string reprepare;
        if (heads[0].potential
            && ! heads[0].potential->prepare (heads[0].grid.grid(), heads[0].grid.spec(), heads[0].grid.radius(), reprepare))
            fail ("batter potential re-prepared after the self-check", reprepare);
    }
    if (options.withAir && ! options.rigidBox)
    {
        const auto shortPulse = measurePmlReflection (domain.delta, k, drum.soundSpeed, drum.airDensity,
                                                      numerics.pmlCells, numerics.pmlReflection, 2.0, pool->size());
        const auto longPulse = measurePmlReflection (domain.delta, k, drum.soundSpeed, drum.airDensity,
                                                     numerics.pmlCells, numerics.pmlReflection, 4.0, pool->size());
        out << format ("- Absorbing layer (Grote-Sim PML, %d cells, theoretical R %.0e): measured reflection of "
                       "Gaussian pulses of 2 and 4 cells at normal and 45/55-degree incidence %.1f and %.1f dB peak "
                       "(%.1f and %.1f dB energy).\n",
                       numerics.pmlCells, numerics.pmlReflection, shortPulse.peakReflectionDb,
                       longPulse.peakReflectionDb, shortPulse.energyReflectionDb, longPulse.energyReflectionDb);
    }
    {
        // Staircase doublet split of the batter: even-m cos/sin pairs in vacuo.
        const auto& head = heads[0].grid;
        const auto modes = clampedStiffModes (head.spec(), head.radius(), 4, 600.0);
        std::vector<std::vector<double>> shapes;
        std::vector<double> initial (static_cast<std::size_t> (head.size()), 0.0);
        std::vector<StiffMode> used;
        for (const auto& mode : modes)
            if (mode.m == 2 || mode.m == 4)
                for (const bool sine : { false, true })
                {
                    auto shape = modeShapeOnGrid (head, mode, sine);
                    const double norm = std::sqrt (head.dot (shape.data(), shape.data()));
                    for (std::size_t l = 0; l < shape.size(); ++l)
                    {
                        shape[l] /= norm;
                        initial[l] += shape[l];
                    }
                    shapes.push_back (std::move (shape));
                    used.push_back (mode);
                }
        HeadGrid lossless;
        lossless.build (head.spec(), head.radius(), k, head.spacing(), heads[0].calibration.inclusionRadius, true);
        const auto traces = simulateInVacuo (lossless, initial, shapes, static_cast<int> (0.3 / k));
        double worst = 0.0;
        for (std::size_t i = 0; i + 1 < used.size(); i += 2)
        {
            const double f = used[i].frequency;
            const double fc = spectralPeakFrequency (traces[i], k, 0.95 * f, 1.03 * f);
            const double fs = spectralPeakFrequency (traces[i + 1], k, 0.95 * f, 1.03 * f);
            worst = std::max (worst, std::abs (fc - fs) / f);
        }
        out << format ("- Staircase doublet split (batter in vacuo, even m <= 4 below 600 Hz): worst cos/sin "
                       "split %.3f %% of f (odd m exactly degenerate by the grid's symmetry).\n", 100.0 * worst);
    }
    out << format ("- Self-checks took %.1f s.\n",
                   std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count());
    checkReport = out.str();
}

template <typename Real>
std::string FdtdEngine<Real>::diagnostics() const
{
    std::ostringstream out;
    const double fs = 1.0 / k;
    if (! failed.empty())
        out << "- FAILED: " << failed << "; the drum stopped and its output is silent.\n";
    out << format ("- Internal rate %.0f Hz (%d x %.0f Hz); output through a linear-phase Kaiser FIR (group delay "
                   "%d internal samples = %.3f ms, cancelled by running the drum that far ahead of its output) "
                   "limiting it to the air grid's valid band: pass edge %.0f Hz, 65 dB stop edge at the on-axis "
                   "cut-off %.0f Hz; %s air.\n",
                   fs, numerics.oversample, outputRate, decimator[0].delayInternalSamples(),
                   1e3 * decimator[0].delayInternalSamples() * k, decimator[0].passEdgeHz(),
                   decimator[0].stopEdgeHz(), sizeof (Real) == 4 ? "float32" : "double");
    const char* names[2] = { "Batter", "Rear" };
    for (int index = 0; index < 2; ++index)
    {
        const auto& head = heads[index];
        const auto& spec = head.grid.spec();
        const double c2 = spec.tension / spec.arealDensity;
        // A sized batter also names the spacing the Berger letter's batter
        // has (its energy bound), so the two FDTD letters' one difference
        // besides the potential is on the record (FdtdEngineSizing.cpp).
        const double bergerBound = bergerTensionBound (spec, drumSpec.radius, headEnergy);
        const double bergerSpacing = std::max (headSpacing (spec, k, bergerBound), options.minimumHeadSpacing);
        const std::string bound = index == 0 && batterDesign > 0.0
                                      ? format ("the peak stretching tension N* the stroke's sizing runs measured "
                                                "(%.4f T; Berger's energy bound for a %.4f J drum would be %.0f N/m, "
                                                "h %.4f mm, the Berger letter's batter spacing, %.2f %% finer; what "
                                                "that difference alone does: tests/integration gridParity)",
                                                batterDesign / spec.tension, headEnergy, bergerBound,
                                                1e3 * bergerSpacing, 100.0 * (1.0 - bergerSpacing / head.grid.spacing()))
                                      : format ("the largest SAV tension a %.4f J drum allows", headEnergy);
        out << format ("- %s head: h %.4f mm (stability bound at T + %.0f N/m, %s, x 1.000001; the linear bound "
                       "is %.4f mm), a/h %.2f, %d interior nodes, "
                       "membrane Courant c k / h %.4f; staircase correction: inclusion radius a - %.3f h "
                       "moves the in-vacuo (0,1) from %.4f Hz (%+.3f %%) to %.4f Hz against the analytic "
                       "clamped %.4f Hz; nonlinearity %s; coupling mu %.3e, nu %.3e, air-plane cells %d.\n",
                       names[index], 1e3 * head.grid.spacing(), tensionBound (index), bound.c_str(),
                       1e3 * stableHeadSpacing (spec, k), head.grid.radius() / head.grid.spacing(),
                       head.grid.interiorCount(), std::sqrt (c2) * k / head.grid.spacing(),
                       (head.grid.radius() - head.calibration.inclusionRadius) / head.grid.spacing(),
                       head.calibration.uncorrectedHz,
                       100.0 * (head.calibration.uncorrectedHz / std::max (1e-9, head.calibration.analyticHz) - 1.0),
                       head.calibration.correctedHz, head.calibration.analyticHz,
                       head.nonlinear && head.potential ? head.potential->name().c_str() : "off",
                       head.mu, head.nu, head.coupling.planeSize());
        if (head.probe)
            out << format ("  SAV stretching tension N* psi / sqrt(2 V + eps) (StretchingProbe.h) over the take: peak "
                           "%.1f N/m at %.3f ms, %.4f of the %.0f N/m the grid admits; over the first %.0f ms %.1f N/m.\n",
                           head.peakSavTension, 1e3 * head.peakSavTensionTime,
                           head.peakSavTension / std::max (1e-30, tensionBound (index)), tensionBound (index),
                           1e3 * sizingWindow, head.windowSavTension);
        if (head.monitorSav)
        {
            const auto& m = head.consistency;
            const double e0 = strokeEnergy (strokeSpec);
            out << format ("  SAV consistency (psi^{n-1/2} against sqrt(V^n + V^{n-1}) at every step, "
                           "VonKarmanPotential.h): the force psi misstates, as a spurious tension N* (psi / sqrt(2 V) - "
                           "1), peaks at %.3f N/m (%.1e T) at %.3f ms within the first %.0f ms (the contact, whose onset "
                           "the explicit step lags by about a step) and at %.3f N/m (%.1e T) at %.3f s after them; at the "
                           "end of the take %.3f N/m with psi / sqrt(2 V) - 1 = %.3g; largest |psi^2 / 2 - V| %.2e J "
                           "(%.1e of the stroke's energy).\n",
                           m.worstAttack, m.worstAttack / spec.tension, 1e3 * m.worstAttackTime, 1e3 * sizingWindow,
                           m.worstTail, m.worstTail / spec.tension, m.worstTailTime, m.spurious, m.ratio - 1.0,
                           m.worstOffset, m.worstOffset / std::max (1e-30, e0));
        }
        if (head.probe && ! head.probe->report().empty())
            out << "  " << head.probe->report()
                << (head.teamSize > 1 ? format ("; each evaluate shared by a team of %d threads", head.teamSize) : std::string())
                << ".\n";
        if (index == 0 && ! sizingReport.empty())
            out << format ("  Grid sizing (FdtdEngineSizing.cpp, %.1f s): ", sizingSeconds) << sizingReport << ".\n";
    }
    if (options.withAir)
        out << format ("- Air: %d x %d x %d cells (%.2f M), Delta %.3f mm (L / %d), Courant %.4f (limit 0.5774); "
                       "margin %.3f m (%.3f m above the batter) to a %d-cell Grote-Sim PML (theoretical R %.0e); %d irregular cells; "
                       "Neumann order %d.\n",
                       domain.nx, domain.ny, domain.nz, domain.cells() / 1e6, 1e3 * domain.delta,
                       numerics.airCellsPerDepth, domain.courant(), numerics.airMargin, numerics.airMarginAbove, numerics.pmlCells,
                       numerics.pmlReflection, domain.irregularCount(), numerics.neumannOrder);
    if (wallActive)
    {
        // Kirchhoff's plane-wave decay in this barrel, (1/a) sqrt(w nu / 2)(1 + (gamma - 1) / sqrt(Pr)), in dB/s.
        const Viscothermal gas;
        const auto planeWave = [&] (double hz)
        {
            const double nu = gas.viscosity / drumSpec.airDensity;
            return 20.0 / std::log (10.0) * std::sqrt (2.0 * 3.14159265358979323846 * hz * nu / 2.0)
                   * (1.0 + (gas.gamma - 1.0) / std::sqrt (gas.prandtl)) / drumSpec.radius;
        };
        out << format ("- Barrel wall: Kirchhoff viscothermal boundary layer on the interior side wall (WallLoss.h; "
                       "viscosity %.3g Pa s, Prandtl %.2f, gamma %.2f), %d wall cells whose staircase faces carry "
                       "2 pi a L = %.4f m^2 (%.4f before normalising), shear on %d elements (axial per z face, "
                       "azimuthal per wall cell); sqrt(j w) by %d diffusive branches (poles %.0f Hz - %.0f kHz) fitted "
                       "in the discrete domain, loss within %.2f %% (thermal) and %.2f %% (viscous) over %.0f Hz - %.0f kHz; "
                       "plane-wave decay %.1f / %.1f / %.1f / %.1f dB/s at 3 / 4 / 5 / 6 kHz (tests/fdtd wallDuctModes: "
                       "axial modes within 0.4 %% to first order, single transverse modes within 5 %%). The heads' ring "
                       "of plane cells carries it inside the coupling solve. Bulk atmospheric absorption is not "
                       "modelled (it needs a humidity the drum spec does not carry).\n",
                       gas.viscosity, gas.prandtl, gas.gamma, wall.wallCellCount(), wall.wallArea(), wall.rawWallArea(),
                       wall.frictionCount(), wall.poles(), wallFit.xi.front() / (2.0 * 3.14159265358979323846),
                       1e-3 * wallFit.xi.back() / (2.0 * 3.14159265358979323846), 100.0 * wallFit.thermalError,
                       100.0 * wallFit.viscousError, wallFit.bandLow, 1e-3 * wallFit.errorHigh, planeWave (3000.0),
                       planeWave (4000.0), planeWave (5000.0), planeWave (6000.0));
    }
    else if (options.withAir)
        out << "- Barrel wall: rigid and lossless (no boundary layer).\n";
    out << format ("- Stick: tip mass %.5f kg, Hertz K %.4e N/m^%.2f, footprint radius %.2f mm over %.2f nodes' "
                   "cells, impact speed %.4f m/s at (%.4f, %.4f) m.\n",
                   strokeSpec.tipMass, strokeSpec.hertzStiffness, strokeSpec.hertzExponent,
                   1e3 * strokeSpec.footprintRadius, heads[0].footprintCells, strokeSpec.impactSpeed,
                   strokeSpec.strikeX, strokeSpec.strikeY);
    const double mv = strokeSpec.tipMass * strokeSpec.impactSpeed;
    const auto figures = forceFigures (stats.force, k);
    // The resolution study quoted after the take's own figures belongs to
    // the batter's potential: the Berger letter's, or the von Karman
    // letter's, whose stretching stiffens the patch under the stick. The von
    // Karman figures are the render tool's dev take of the Don at the
    // default rate and with TAIKOR_FDTD_OVERSAMPLE=4 (192 kHz; the batter
    // sized at each rate, h 4.9187 / 3.4074 mm): chatter 9.44 / 10.10 kHz,
    // impulse 0.27115 / 0.27128 N s, peak force 130.41 / 131.91 N, 5-95 %
    // span 0.302-3.688 / 0.292-3.688 ms (first contact one 96 kHz step
    // later), peak N* 9008.5 / 8938.7 N/m.
    const auto* batterProbe = dynamic_cast<const StretchingProbe*> (heads[0].potential.get());
    const bool stretchingSized = batterProbe != nullptr && ! batterProbe->boundedByEnergy();
    const char* resolution = stretchingSized
                                 ? "that frequency depends on the resolution (the von Karman Don: 9.44 kHz at 96 kHz "
                                   "internal, 10.10 kHz at 192 kHz), while its impulse agrees to 0.05 %, its 5-95 % "
                                   "span to one 96 kHz step and its peak force to 1.2 % between the two rates"
                                 : "that frequency depends on the resolution (the Don: 8.83 kHz at 96 kHz internal, "
                                   "9.45 kHz at 192 kHz, 9.55 kHz in a continuum membrane/disk/Hertz reference), "
                                   "while its impulse and 5-95 % span agree to 0.02 % and its peak force to 1.2 % "
                                   "between the two rates";
    out << format ("- Contact: %d episode(s); first contact %.4f ms after the strike, last %.3f ms; "
                   "peak force %.2f N; impulse %.5f N s = %.3f m v0; 5-95 %% of the impulse between %.3f and "
                   "%.3f ms; stick exit speed %.3f m/s. The force carries the footprint patch's oscillation on "
                   "the Hertz spring at %.2f kHz, %.1f dB re the mean contact force (Decimator.h); %s.\n",
                   stats.episodes, 1e3 * stats.firstContact, 1e3 * stats.lastContact, stats.peakForce,
                   stats.impulse, stats.impulse / std::max (1e-30, mv), 1e3 * figures.t5, 1e3 * figures.t95,
                   stats.exitSpeed, 1e-3 * figures.oscillationHz, figures.oscillationDb, resolution);
    if (strokeBeyondBound)
        out << format ("- WARNING: a stroke on the moving drum brought more energy than the %.4f J the head grids "
                       "admit; the SAV step may lose track of the nonlinear potential.\n", headEnergy);
    out << checkReport;
    if (profiling)
    {
        const char* names[10] = { "potential", "air planes", "head explicit", "air regular", "air irregular", "PML",
                                  "P2a", "P2b", "P3", "barrier wait" };
        out << "- Profile (thread-seconds):";
        for (int kind = 0; kind < 10; ++kind)
        {
            double total = 0.0;
            for (const auto& p : profile)
                total += p[static_cast<std::size_t> (kind)];
            out << format (" %s %.2f;", names[kind], total);
        }
        out << "\n";
    }
    out << format ("- Cost: prepare %.1f s (grids, staircase calibration, coupling), head rebuild for the stroke "
                   "%.1f s, %lld internal steps in %.1f s (%.3f ms per step, %d threads).",
                   prepareSeconds, rebuildSeconds, processedSteps, processSeconds,
                   processedSteps > 0 ? 1e3 * processSeconds / processedSteps : 0.0, pool ? pool->size() : 0);
    return out.str();
}

#define TAIKOR_REPORT_INSTANTIATE(Real)                                                              \
    template double FdtdEngine<Real>::energyHeads() const;                                           \
    template double FdtdEngine<Real>::energyAir() const;                                             \
    template double FdtdEngine<Real>::dissipatedWall() const;                                        \
    template double FdtdEngine<Real>::energySav() const;                                             \
    template double FdtdEngine<Real>::energyStick() const;                                           \
    template void FdtdEngine<Real>::runSelfChecks (const DrumSpec&, const MicSpec&, double);         \
    template std::string FdtdEngine<Real>::diagnostics() const;
TAIKOR_REPORT_INSTANTIATE (float)
TAIKOR_REPORT_INSTANTIATE (double)
#undef TAIKOR_REPORT_INSTANTIATE

} // namespace taikor::candidates::fdtd
