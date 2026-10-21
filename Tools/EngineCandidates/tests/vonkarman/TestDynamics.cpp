// (6) SAV time-stepping with the potential conserves energy, (6b) the output
// does not move when the in-plane tolerance is tightened, and (7) the cost of
// one evaluate() on the drum's grids.

#include "InVacuoHead.h"

#include <algorithm>

namespace vktest
{

namespace
{
// Forwards to a potential and times every call.
struct TimedPotential final : NonlinearPotential
{
    NonlinearPotential* inner = nullptr;
    double seconds = 0.0;
    double worst = 0.0;
    long calls = 0;
    std::string name() const override { return inner->name(); }
    bool prepare (const MembraneGrid&, const HeadSpec&, double, std::string&) override { return true; }
    double evaluate (const double* w, double* gradient) override
    {
        const double t0 = vktest::cpuSeconds();
        const double v = inner->evaluate (w, gradient);
        const double dt = vktest::cpuSeconds() - t0;
        seconds += dt;
        worst = std::max (worst, dt);
        ++calls;
        return v;
    }
};

std::vector<double> dentAndBowl (const MembraneGrid& g, double a)
{
    const auto& s = shipping().stroke;
    return sampleField (g, [&] (double x, double y)
    {
        const double dx = x - s.strikeX, dy = y - s.strikeY;
        const double bowl = 1.0 - (x * x + y * y) / (a * a);
        return 0.008 * std::exp (-(dx * dx + dy * dy) / (0.04 * 0.04)) * bowl + 0.003 * bowl;
    });
}

struct StrokeRun
{
    std::vector<double> probe, force, finalField, volumeVelocity;
    double evalSeconds = 0.0, worstCall = 0.0;
    long calls = 0;
    VonKarmanPotential::Stats stats;
    double peakDisplacement = 0.0;
    double peakV = 0.0, psiGap = 0.0; // max |psi - sqrt(2V + eps)|
    double peakStress = 0.0;          // max mean stretching stress, N/m
    std::string diagnostics;          // VonKarmanPotential::diagnostics() at the end
};

// The test stroke (tip mass, speed, footprint and Hertz law of the resolved
// StrokeSpec) on an in-vacuo lossless batter at spacing h and sample rate fs.
StrokeRun runStroke (double h, double fs, int steps, double tolerance, bool nonlinear, bool tensionRelative, bool loss = false)
{
    const double a = headRadius();
    StrokeRun out;
    Fixture f = makeFixture (a, h, tolerance, tensionRelative);
    if (! f.vk)
    {
        check (false, "prepare: " + f.error);
        return out;
    }
    TimedPotential timed;
    timed.inner = f.vk.get();
    InVacuoHead head;
    head.init (f.grid, batter(), fs, nonlinear ? &timed : nullptr);
    head.loss = loss;
    head.setDisplacement (std::vector<double> (static_cast<std::size_t> (f.grid.size()), 0.0));
    head.armStick (shipping().stroke);
    // Probe node: a point away from the stick, (-0.15, 0.10) m.
    const int pi = static_cast<int> (std::lround ((-0.15 - f.grid.x0) / f.grid.h));
    const int pj = static_cast<int> (std::lround ((0.10 - f.grid.y0) / f.grid.h));
    const auto probe = static_cast<std::size_t> (f.grid.index (pi, pj));
    for (int n = 0; n < steps; ++n)
    {
        head.step();
        out.probe.push_back (head.w[probe]);
        out.force.push_back (head.contactForce);
        double q = 0.0; // volume velocity sum (w^{n+1} - w^n) h^2 / k, the head's monopole source
        for (std::size_t k = 0; k < head.w.size(); ++k)
            q += head.w[k] - head.wm[k];
        out.volumeVelocity.push_back (q * f.grid.h * f.grid.h * fs);
        for (double v : head.w)
            out.peakDisplacement = std::max (out.peakDisplacement, std::abs (v));
        if (nonlinear)
        {
            out.peakV = std::max (out.peakV, head.lastV);
            out.psiGap = std::max (out.psiGap, std::abs (head.psi - std::sqrt (2.0 * head.lastV + head.epsR)));
            for (int k = 0; k < f.grid.size(); ++k)
                out.peakStress = std::max (out.peakStress, 0.5 * (f.vk->stressXX()[k] + f.vk->stressYY()[k]));
        }
    }
    out.finalField = head.w;
    out.evalSeconds = timed.seconds;
    out.worstCall = timed.worst;
    out.calls = timed.calls;
    out.stats = f.vk->stats();
    out.diagnostics = f.vk->diagnostics();
    return out;
}

// Max |x - ref| over each window of `window` samples, relative to max |ref|
// over the whole run.
std::string windowedDifference (const std::vector<double>& x, const std::vector<double>& ref, std::size_t window)
{
    double den = 0.0;
    for (double v : ref)
        den = std::max (den, std::abs (v));
    std::string out;
    for (std::size_t w0 = 0; w0 < ref.size(); w0 += window)
    {
        double num = 0.0;
        for (std::size_t k = w0; k < std::min (ref.size(), w0 + window); ++k)
            num = std::max (num, std::abs (x[k] - ref[k]));
        out += format (" %.1e", den > 0.0 ? num / den : num);
    }
    return out;
}

double relativeMaxDifference (const std::vector<double>& x, const std::vector<double>& ref)
{
    double num = 0.0, den = 0.0;
    for (std::size_t k = 0; k < ref.size(); ++k)
    {
        num = std::max (num, std::abs (x[k] - ref[k]));
        den = std::max (den, std::abs (ref[k]));
    }
    return den > 0.0 ? num / den : num;
}
} // namespace

void testSav()
{
    const double a = headRadius();
    for (double tolerance : { 1.0e-10, 1.0e-2 })
    {
        Fixture f = makeFixture (a, a / 20.0, tolerance);
        if (! f.vk)
            return check (false, "prepare: " + f.error);
        InVacuoHead head;
        head.init (f.grid, batter(), 48000.0, f.vk.get());
        head.setDisplacement (dentAndBowl (f.grid, a));
        head.step();
        const double e0 = head.energy();
        double drift = 0.0, psiGap = 0.0, vPeak = 0.0;
        for (int n = 0; n < 4800; ++n)
        {
            head.step();
            drift = std::max (drift, std::abs (head.energy() / e0 - 1.0));
            vPeak = std::max (vPeak, head.lastV);
            psiGap = std::max (psiGap, std::abs (head.psi - std::sqrt (2.0 * head.lastV + head.epsR)));
        }
        check (drift < 1e-11, format ("tolerance %.0e: lossless energy conserved over 0.1 s (4800 steps): max |E/E0 - 1| = %.2e (< 1e-11)",
                                      tolerance, drift));
        note (format ("E0 = %.4e J, peak V = %.4e J, max |psi - sqrt(2V + eps)| = %.2e sqrt(J) (%.1e of sqrt(2 V_peak)), mean PCG its %.2f",
                      e0, vPeak, psiGap, psiGap / std::sqrt (2 * vPeak),
                      static_cast<double> (f.vk->stats().iterations) / static_cast<double> (f.vk->stats().calls)));
    }
}

void testTolerance()
{
    // The regime in which the explicit SAV step resolves the stroke: the
    // 48 kHz stability-limit grid stepped at 96 kHz (at 48 kHz the local von
    // Karman tension under the stick exceeds the explicit step's margin and
    // psi decouples from sqrt(2V); see testSavResolution). 50 ms covers the
    // whole contact (about 4.5 ms) and the start of the ring-down.
    const double fs = 96000.0;
    const double h = 1.001 * InVacuoHead::minimumSpacing (batter(), 48000.0);
    const int steps = static_cast<int> (0.05 * fs);
    note (format ("grid h = %.3f mm stepped at %.0f kHz, %d steps", h * 1e3, fs / 1000, steps));
    const StrokeRun ref = runStroke (h, fs, steps, 1.0e-10, true, false);
    const StrokeRun lin = runStroke (h, fs, steps, 1.0e-10, false, false);
    note (format ("reference (tol 1e-10): peak |w| = %.2f mm, peak contact force %.0f N, mean PCG its %.2f, max %d",
                  ref.peakDisplacement * 1e3, *std::max_element (ref.force.begin(), ref.force.end()),
                  static_cast<double> (ref.stats.iterations) / static_cast<double> (ref.stats.calls), ref.stats.maxIterations));
    const double vkEffect = relativeMaxDifference (lin.volumeVelocity, ref.volumeVelocity);
    note (format ("size of the von Karman effect (linear head against von Karman): volume velocity %.2e, probe %.2e, force %.2e",
                  vkEffect, relativeMaxDifference (lin.probe, ref.probe), relativeMaxDifference (lin.force, ref.force)));
    const double shipping = VonKarmanPotential().solverSettings().relativeTolerance;
    const auto window = static_cast<std::size_t> (0.01 * fs);
    {
        // How the dynamics amplify a perturbation: a tighter solve against the
        // reference, per 10 ms window.
        const StrokeRun tight = runStroke (h, fs, steps, 1.0e-12, true, false);
        note ("tol 1e-12 against 1e-10, volume velocity per 10 ms:" + windowedDifference (tight.volumeVelocity, ref.volumeVelocity, window));
    }
    double chosen = -1.0;
    for (bool rule : { false, true })
        for (double tolerance : { 1.0e-2, 3.0e-3, 1.0e-3, 1.0e-4 })
        {
            const StrokeRun run = runStroke (h, fs, steps, tolerance, true, rule);
            const double dQ = relativeMaxDifference (run.volumeVelocity, ref.volumeVelocity);
            const double dProbe = relativeMaxDifference (run.probe, ref.probe);
            const double dForce = relativeMaxDifference (run.force, ref.force);
            note (format ("%s tol %.0e: volume velocity %.2e, probe %.2e, contact force %.2e (max diff / max); PCG its mean %.2f max %d; %.0f us/call",
                          rule ? "tension-relative" : "plain relative  ", tolerance, dQ, dProbe, dForce,
                          static_cast<double> (run.stats.iterations) / static_cast<double> (run.stats.calls), run.stats.maxIterations,
                          1e6 * run.evalSeconds / static_cast<double> (run.calls)));
            note ("    volume velocity per 10 ms:" + windowedDifference (run.volumeVelocity, ref.volumeVelocity, window));
            if (rule == VonKarmanPotential().tensionRelativeStopping() && std::abs (tolerance / shipping - 1.0) < 1e-9)
                chosen = std::max ({ dQ, dProbe, dForce });
        }
    check (chosen >= 0.0 && chosen < 1e-3,
           format ("shipping setting (tol %.0e): output within %.1e of the 1e-10 reference (< 1e-3; the von Karman effect itself is %.1e)",
                   shipping, chosen, vkEffect));
}

void testSavResolution()
{
    // Whether the explicit SAV step resolves the test stroke's local von
    // Karman stiffening: psi must track sqrt(2 V). On a grid at the linear
    // stability limit of its own step the stretching stress under the stick
    // (2-3 T) exceeds the step's margin and psi decouples; with the step
    // halved on the same grid it tracks. This bounds how the drum may run the
    // potential; the potential itself is exact either way (tests 1-6).
    struct Case { double gridRate, stepRate; };
    double tracked = 1.0;
    for (bool loss : { false, true })
        for (const Case c : { Case { 48000.0, 96000.0 }, Case { 48000.0, 48000.0 }, Case { 96000.0, 96000.0 }, Case { 96000.0, 192000.0 } })
        {
            const double h = 1.001 * InVacuoHead::minimumSpacing (batter(), c.gridRate, loss);
            const StrokeRun run = runStroke (h, c.stepRate, static_cast<int> (0.006 * c.stepRate), 1.0e-8, true, false, loss);
            const double gap = run.psiGap / std::sqrt (2.0 * run.peakV);
            note (format ("%s grid h = %.3f mm stepped at %3.0f kHz: peak V %.4f J, peak stretching stress %5.0f N/m (%.1f T), max |psi - sqrt(2V)| / sqrt(2 V_peak) = %.3f",
                          loss ? "hide loss:" : "lossless: ", h * 1e3, c.stepRate / 1000, run.peakV, run.peakStress, run.peakStress / batter().tension, gap));
            if (! loss && c.gridRate == 48000.0 && c.stepRate == 96000.0)
                tracked = gap;
        }
    check (tracked < 0.02, format ("48 kHz grid stepped at 96 kHz: psi tracks sqrt(2V) to %.3f (< 0.02)", tracked));
}

void testBenchmark()
{
    // Cost of evaluate() in the shipping setup on the drum's two grids: the
    // 48 kHz stability-limit grid (h ~ 6.9 mm) and the 96 kHz one (h ~ 4.8 mm),
    // both stepped at 96 kHz over the 50 ms that hold the whole contact, and
    // the 48 kHz grid at 48 kHz for comparison.
    struct Case { double gridRate, stepRate, target; };
    for (const Case c : { Case { 48000.0, 96000.0, 0.3e-3 }, Case { 96000.0, 96000.0, 1.0e-3 }, Case { 48000.0, 48000.0, 0.3e-3 } })
    {
        const double h = 1.001 * InVacuoHead::minimumSpacing (batter(), c.gridRate);
        const int steps = static_cast<int> (0.05 * c.stepRate);
        const VonKarmanPotential shippingSetup;
        const StrokeRun run = runStroke (h, c.stepRate, steps, shippingSetup.solver().settings().relativeTolerance, true,
                                         shippingSetup.tensionRelativeStopping());
        const MembraneGrid g = makeDiskGrid (headRadius(), h);
        int interior = 0;
        for (auto flag : g.interior)
            interior += flag;
        const double perCall = run.evalSeconds / static_cast<double> (run.calls);
        note (format ("grid h = %.3f mm (%d x %d, %d interior nodes), stepped at %.0f kHz: %.0f us per evaluate (worst %.0f us), PCG its mean %.2f max %d",
                      h * 1e3, g.nx, g.ny, interior, c.stepRate / 1000, perCall * 1e6, run.worstCall * 1e6,
                      static_cast<double> (run.stats.iterations) / static_cast<double> (run.stats.calls), run.stats.maxIterations));
        note (format ("   target %.1f ms: %s (thread CPU time on a shared machine)", c.target * 1e3, perCall <= c.target ? "met" : "NOT met"));
        note ("   " + run.diagnostics);
        // Loose regression bound only (the machine is shared): 10x the target.
        check (perCall < 10.0 * c.target, format ("%.0f us per evaluate < 10 x target", perCall * 1e6));
    }
}

} // namespace vktest
