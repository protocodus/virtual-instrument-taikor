// (5) lossless energy conservation of stick + bank + Berger through a
// collision; (6) mass on a rigid barrier against Hertz's contact time;
// (7) convergence of the contact impulse with oversampling and mode count;
// and the numerical tolerances of the bank (Berger sparsity, pruning).

#include "ModalBuild.h"
#include "ModalFixtures.h"
#include "ModalRuntime.h"
#include "ModalTestUtil.h"

#include <algorithm>
#include <map>
#include <memory>

using namespace taikor::candidates;
using namespace taikor::candidates::modal;

namespace
{
// Models of the nagado at a given batter band, built once per band.
const ModalModel& nagadoModel (double headMaxHz)
{
    static std::map<double, std::unique_ptr<ModalModel>> cache;
    auto& slot = cache[headMaxHz];
    if (! slot)
    {
        BuildOptions o;
        o.headMaxHz = headMaxHz;
        o.coupledMaxHz = std::min (o.coupledMaxHz, headMaxHz);
        std::string error;
        slot = buildModel (modaltest::nagadoDrum(), modaltest::nagadoMics(), o, error);
    }
    return *slot;
}

struct Run
{
    ContactReport report;
    std::vector<double> left;
    double exitVelocity = 0.0;
};

Run strikeAndRun (const ModalModel& model, const RuntimeOptions& ro, double seconds, bool lossless = false)
{
    Bank bank = makeBank (model, modaltest::nagadoDrum(), modaltest::nagadoMics(), modaltest::nagadoStroke(), ro);
    if (lossless)
    {
        for (auto& m : bank.modes)
            m.alpha = 0.0;
        bank.clusters.clear();
    }
    Stepper stepper (std::move (bank), 48000.0, ro);
    stepper.strike();
    Run run;
    run.left.resize (static_cast<std::size_t> (seconds * 48000.0));
    stepper.processDouble (run.left.data(), nullptr, static_cast<int> (run.left.size()));
    run.report = stepper.report();
    run.exitVelocity = stepper.stickVelocity();
    return run;
}
} // namespace

MODAL_TEST (contact_rigid_barrier)
{
    // No modes: the tip mass meets a rigid half-space through the Hertz law.
    // Hertz (Johnson, Contact Mechanics, 1985, sec. 11.4): maximum compression
    // d = (5 m v^2 / (4 K))^(2/5) and contact time t = 2.9432 d / v, where
    // 2.9432 = 2 int_0^1 dx / sqrt(1 - x^(5/2)); the rebound is elastic.
    const StrokeSpec st = modaltest::nagadoStroke();
    Bank bank;
    bank.tipMass = st.tipMass;
    bank.hertzStiffness = st.hertzStiffness;
    bank.hertzExponent = st.hertzExponent;
    bank.impactSpeed = st.impactSpeed;
    RuntimeOptions ro;
    Stepper stepper (bank, 48000.0, ro);
    stepper.strike();
    std::vector<double> out (480);
    stepper.processDouble (out.data(), nullptr, 480);
    const auto& r = stepper.report();
    const double v = st.impactSpeed, m = st.tipMass, k = st.hertzStiffness;
    const double depth = std::pow (5.0 * m * v * v / (4.0 * k), 0.4);
    const double tc = 2.9432 * depth / v;
    const double dt = 1.0 / (48000.0 * ro.oversample);
    const double measured = r.lastContact - r.firstContact + dt;
    modaltest::nearRel ("contact time vs Hertz 2.9432 d/v", measured, tc, 3.0 * dt / tc);
    modaltest::nearRel ("maximum compression vs (5 m v^2 / 4K)^(2/5)", r.maxPenetration, depth, 2e-3);
    modaltest::nearRel ("rebound speed (elastic)", stepper.stickVelocity(), v, 1e-10);
    modaltest::nearRel ("impulse / (2 m v)", r.impulse / (2.0 * m * v), 1.0, 1e-10);
}

MODAL_TEST (contact_energy_lossless)
{
    const ModalModel& model = nagadoModel (8000.0);
    // Base rate only (one numerical energy throughout), damping removed, Berger on.
    RuntimeOptions ro;
    ro.forceRate = 1;
    ro.oversample = 1;
    Bank bank = makeBank (model, modaltest::nagadoDrum(), modaltest::nagadoMics(), modaltest::nagadoStroke(), ro);
    for (auto& m : bank.modes)
        m.alpha = 0.0;
    bank.clusters.clear();
    Stepper stepper (std::move (bank), 48000.0, ro);
    stepper.strike();
    const double e0 = stepper.totalEnergy();
    double worst = 0.0;
    for (int i = 0; i < 960; ++i)
    {
        double l = 0.0;
        stepper.processDouble (&l, nullptr, 1);
        worst = std::max (worst, std::abs (stepper.totalEnergy() - e0) / e0);
    }
    std::printf ("  %d modes, %d contact episode(s), peak Berger tension rise %.1f N/m\n", stepper.modeCount(),
                 stepper.report().episodes, stepper.report().peakTensionRise[0]);
    modaltest::below ("max relative energy drift over 20 ms through the collision (lossless)", worst, 1e-12);
    // Oversampled throughout (rate conversions every sample): the window balance.
    RuntimeOptions over;
    over.forceRate = 2;
    over.oversample = 8;
    const Run run = strikeAndRun (model, over, 0.02, true);
    modaltest::below ("oversampled x8: |E_end + dissipated - E_start| / E_start (lossless)",
                      std::abs (run.report.energyEnd + run.report.dissipated - run.report.energyStart) / run.report.energyStart, 1e-12);
}

MODAL_TEST (contact_impulse_convergence)
{
    const double mv = modaltest::nagadoStroke().tipMass * modaltest::nagadoStroke().impactSpeed;
    // Oversampling on the full-band model.
    const ModalModel& full = nagadoModel (21600.0);
    std::vector<double> impulse;
    for (int os : { 4, 8, 16, 32 })
    {
        RuntimeOptions ro;
        ro.oversample = os;
        const Run run = strikeAndRun (full, ro, 0.012);
        impulse.push_back (run.report.impulse / mv);
        std::printf ("  x%-2d: impulse %.6f m v0, peak force %.2f N, contact %.3f-%.3f ms, %d episode(s)\n", os,
                     impulse.back(), run.report.peakForce, run.report.firstContact * 1e3, run.report.lastContact * 1e3,
                     run.report.episodes);
    }
    modaltest::below ("impulse change x16 -> x32 (relative)", std::abs (impulse[2] - impulse[3]) / impulse[3], 2e-3);
    modaltest::check (std::abs (impulse[2] - impulse[3]) <= std::abs (impulse[1] - impulse[2]) + 1e-6,
                      "impulse differences shrink with oversampling", std::abs (impulse[2] - impulse[3]),
                      std::abs (impulse[1] - impulse[2]), 0.0);
    // Mode count: batter band 4 -> 8 -> 12 -> 21.6 kHz at x16.
    std::vector<double> byBand;
    for (double band : { 4000.0, 8000.0, 12000.0, 21600.0 })
    {
        RuntimeOptions ro;
        const Run run = strikeAndRun (nagadoModel (band), ro, 0.012);
        byBand.push_back (run.report.impulse / mv);
        std::printf ("  band %5.0f Hz: impulse %.6f m v0, peak force %.2f N\n", band, byBand.back(), run.report.peakForce);
    }
    modaltest::below ("impulse change 12 -> 21.6 kHz band (relative)", std::abs (byBand[2] - byBand[3]) / byBand[3], 5e-3);
}

MODAL_TEST (modal_berger_sparsity)
{
    const ModalModel& model = nagadoModel (8000.0);
    RuntimeOptions dense;
    dense.bergerTolerance = 0.0;
    dense.bergerRankTolerance = 0.0;
    RuntimeOptions sparse;
    const Run a = strikeAndRun (model, dense, 0.1), b = strikeAndRun (model, sparse, 0.1);
    modaltest::below ("peak batter tension rise: sparse vs dense Berger form (relative)",
                      std::abs (a.report.peakTensionRise[0] - b.report.peakTensionRise[0]) / a.report.peakTensionRise[0], 2e-3);
    double diff = 0.0, ref = 0.0;
    for (std::size_t i = 0; i < a.left.size(); ++i)
    {
        diff += (a.left[i] - b.left[i]) * (a.left[i] - b.left[i]);
        ref += a.left[i] * a.left[i];
    }
    modaltest::below ("capsule signal, first 100 ms: rms difference / rms", std::sqrt (diff / ref), 1e-3);
}

MODAL_TEST (modal_pruning)
{
    const ModalModel& model = nagadoModel (8000.0);
    RuntimeOptions all;
    all.pruneResidue = 0.0;
    all.pruneContact = 0.0;
    RuntimeOptions pruned;
    const Run a = strikeAndRun (model, all, 0.2), b = strikeAndRun (model, pruned, 0.2);
    double diff = 0.0, ref = 0.0;
    for (std::size_t i = 0; i < a.left.size(); ++i)
    {
        diff += (a.left[i] - b.left[i]) * (a.left[i] - b.left[i]);
        ref += a.left[i] * a.left[i];
    }
    modaltest::below ("capsule signal, first 200 ms: rms difference pruned vs all / rms", std::sqrt (diff / ref), 1e-5);
    modaltest::nearRel ("contact impulse pruned vs all", b.report.impulse, a.report.impulse, 1e-5);
}

MODAL_TEST (build_damping_coupling)
{
    // Damping couplings are kept between modes closer than B (alpha_k + alpha_l).
    // Band envelopes (20 ms frames, 40 Hz-8 kHz in five bands, within 50 dB of
    // each band's peak) for B = 8 against B = 16, and for first-order projection
    // alone (no clusters) against B = 16. Waveforms are not compared: a coupling
    // across a wide gap pulls a frequency by ~Ct^2 / (4 dw), a phase drift that
    // changes no envelope.
    // A fourth build lowers dampingFloor tenfold (couplings kept down to 1e-4
    // sqrt(alpha_k alpha_l) instead of 1e-3).
    std::vector<double> signal[4];
    const double bandwidth[4] = { 8.0, 16.0, 8.0, 8.0 };
    const double floor[4] = { 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-4 };
    for (int i = 0; i < 4; ++i)
    {
        BuildOptions o;
        o.headMaxHz = 8000.0;
        o.dampingBandwidth = bandwidth[i];
        o.dampingFloor = floor[i];
        std::string error;
        const auto model = buildModel (modaltest::nagadoDrum(), modaltest::nagadoMics(), o, error);
        Bank bank = makeBank (*model, modaltest::nagadoDrum(), modaltest::nagadoMics(), modaltest::nagadoStroke(), RuntimeOptions {});
        if (i == 2)
            bank.clusters.clear();
        Stepper stepper (std::move (bank), 48000.0, RuntimeOptions {});
        stepper.strike();
        signal[i].resize (24000);
        stepper.processDouble (signal[i].data(), nullptr, 24000);
    }
    // Butterworth-free band split: FFT masks per 20 ms frame (Hann, 50 % overlap).
    const auto envelopeError = [&] (const std::vector<double>& x, const std::vector<double>& ref)
    {
        const int frame = 960, hop = 480;
        const double bands[6] = { 40.0, 250.0, 1000.0, 2000.0, 4000.0, 8000.0 };
        std::vector<std::vector<double>> ex (5), er (5);
        for (int start = 0; start + frame <= static_cast<int> (x.size()); start += hop)
        {
            double sx[5] = { 0, 0, 0, 0, 0 }, sr[5] = { 0, 0, 0, 0, 0 };
            for (int kbin = 1; kbin < frame / 2; ++kbin)
            {
                const double f = kbin * 48000.0 / frame;
                int band = -1;
                for (int b = 0; b < 5; ++b)
                    if (f >= bands[b] && f < bands[b + 1])
                        band = b;
                if (band < 0)
                    continue;
                std::complex<double> ax = 0.0, ar = 0.0;
                for (int nidx = 0; nidx < frame; ++nidx)
                {
                    const double w = 0.5 - 0.5 * std::cos (2.0 * pi * nidx / frame);
                    const std::complex<double> e = std::polar (1.0, -2.0 * pi * kbin * nidx / frame);
                    ax += w * x[static_cast<std::size_t> (start + nidx)] * e;
                    ar += w * ref[static_cast<std::size_t> (start + nidx)] * e;
                }
                sx[band] += std::norm (ax);
                sr[band] += std::norm (ar);
            }
            for (int b = 0; b < 5; ++b)
            {
                ex[static_cast<std::size_t> (b)].push_back (10.0 * std::log10 (sx[b] + 1e-300));
                er[static_cast<std::size_t> (b)].push_back (10.0 * std::log10 (sr[b] + 1e-300));
            }
        }
        double worst = 0.0;
        for (int b = 0; b < 5; ++b)
        {
            const auto& r = er[static_cast<std::size_t> (b)];
            const double top = *std::max_element (r.begin(), r.end());
            for (std::size_t i = 0; i < r.size(); ++i)
                if (r[i] > top - 50.0)
                    worst = std::max (worst, std::abs (ex[static_cast<std::size_t> (b)][i] - r[i]));
        }
        return worst;
    };
    const double clustered = envelopeError (signal[0], signal[1]);
    const double projection = envelopeError (signal[2], signal[1]);
    const double floorError = envelopeError (signal[0], signal[3]);
    std::printf ("  first-order projection alone vs B = 16: max band-envelope difference %.2f dB\n", projection);
    modaltest::below ("B = 8 vs B = 16: max band-envelope difference (dB)", clustered, 0.5);
    modaltest::below ("damping floor 1e-3 vs 1e-4: max band-envelope difference (dB)", floorError, 0.5);
}
