// The capsule pressure while the stick drives the head (ModalFeedthrough.cpp):
// its time integral over a stroke vanishes (the head starts and ends at rest),
// nothing arrives before sound can, and the output equals a brute-force
// retarded Rayleigh sum of the modal accelerations, forcing included.

#include "ModalBuild.h"
#include "ModalFixtures.h"
#include "ModalRuntime.h"
#include "ModalTestUtil.h"

#include <algorithm>
#include <cmath>
#include <memory>

using namespace taikor::candidates;
using namespace taikor::candidates::modal;

namespace
{
const ModalModel& fullModel()
{
    static std::unique_ptr<ModalModel> model;
    if (! model)
    {
        std::string error;
        model = buildModel (modaltest::nagadoDrum(), modaltest::nagadoMics(), BuildOptions {}, error);
    }
    return *model;
}

struct Take
{
    std::vector<double> left, right;
    std::array<double, 2> tail {};
    double lastContact = 0.0;
    int baseRateContactSteps = 0;
};

Take render (const ModalModel& model, const RuntimeOptions& ro, double seconds)
{
    Bank bank = makeBank (model, modaltest::nagadoDrum(), modaltest::nagadoMics(), modaltest::nagadoStroke(), ro, 48000.0);
    Stepper stepper (std::move (bank), 48000.0, ro);
    stepper.strike();
    Take take;
    const auto count = static_cast<std::size_t> (seconds * 48000.0);
    take.left.resize (count);
    take.right.resize (count);
    stepper.processDouble (take.left.data(), take.right.data(), static_cast<int> (count));
    take.tail = stepper.freeDecayIntegral();
    take.lastContact = stepper.report().lastContact;
    take.baseRateContactSteps = stepper.report().baseRateContactSteps;
    return take;
}

// int_0^T p dt by the trapezoid (p = 0 at the strike) plus the free decay's rest.
double pressureIntegral (const std::vector<double>& p, double tail)
{
    double sum = 0.0;
    for (double v : p)
        sum += v;
    return (sum - 0.5 * p.back()) / 48000.0 + tail;
}
} // namespace

MODAL_TEST (output_pressure_integral)
{
    // int p dt = rho int [w_t] / (2 pi R) dS = 0 for a head at rest before and
    // after (the rear monopole likewise). Measured over 0.25 s plus the analytic
    // integral of the remaining free decay, on the 8 kHz model. With Berger off
    // the contact force is the only external force; with it on, the Berger
    // force's own direct term (not in the output) leaves a remainder.
    std::string error;
    BuildOptions o;
    o.headMaxHz = 8000.0;
    const auto model = buildModel (modaltest::nagadoDrum(), modaltest::nagadoMics(), o, error);
    for (const bool berger : { false, true })
    {
        RuntimeOptions with, without;
        with.berger = without.berger = berger;
        without.feedthrough = false;
        const Take a = render (*model, with, 0.25), b = render (*model, without, 0.25);
        const double ia[2] = { pressureIntegral (a.left, a.tail[0]), pressureIntegral (a.right, a.tail[1]) };
        const double ib[2] = { pressureIntegral (b.left, b.tail[0]), pressureIntegral (b.right, b.tail[1]) };
        std::printf ("  Berger %s: int p dt without the contact-force term %.5f / %.5f Pa s, with it %.3g / %.3g Pa s\n",
                     berger ? "on" : "off", ib[0], ib[1], ia[0], ia[1]);
        // The contact lookahead (ModalRuntime.cpp, 3 samples of the gap's rate) must
        // put every contact step on the substep grid.
        modaltest::below ("contact steps taken at the base rate", a.baseRateContactSteps, 0.0);
        if (! berger)
        {
            modaltest::below ("Berger off: |int p dt| with the kernel / without, left", std::abs (ia[0]) / std::abs (ib[0]), 2e-3);
            modaltest::below ("Berger off: |int p dt| with the kernel / without, right", std::abs (ia[1]) / std::abs (ib[1]), 2e-3);
        }
        else
            modaltest::note ("Berger on: |int p dt| with the kernel / without (Berger direct term), left: %.4f",
                             std::abs (ia[0]) / std::abs (ib[0]));
    }
}

namespace
{
// The capsule pressure as the retarded Rayleigh sum of the modal accelerations,
// p(t) = rho sum_k sum_nodes w y_k''(t - R/c) (+ the rear monopole), from a run
// at fs * oversample with no oversampling, whose base step is the 48 kHz run's
// contact substep, against that 48 kHz run. Returns the rms and peak of the
// reference over the first `samples` output samples and the rms difference.
struct Comparison
{
    double rmsReference = 0.0, peakReference = 0.0, rmsError = 0.0, rmsErrorWithout = 0.0;
};

Comparison retardedReference (bool berger, int samples)
{
    BuildOptions o;
    o.headMaxHz = 4000.0;
    o.coupledMaxHz = 1500.0;
    o.airMaxHz = 2000.0;
    std::string error;
    const auto model = buildModel (modaltest::nagadoDrum(), modaltest::nagadoMics(), o, error);
    const auto drum = modaltest::nagadoDrum();
    const auto mics = modaltest::nagadoMics();
    const auto stroke = modaltest::nagadoStroke();
    RuntimeOptions ro;
    ro.berger = berger;
    ro.pruneResidue = ro.pruneContact = 0.0;
    RuntimeOptions without = ro;
    without.feedthrough = false;
    const int os = ro.oversample;
    const double fine = 48000.0 * os, delta = 1.0 / fine;
    // The 48 kHz engine, with and without the contact-force term.
    std::vector<double> engine[2], plain[2];
    for (int pass = 0; pass < 2; ++pass)
    {
        Bank bank = makeBank (*model, drum, mics, stroke, pass == 0 ? ro : without, 48000.0);
        Stepper stepper (std::move (bank), 48000.0, pass == 0 ? ro : without);
        stepper.strike();
        auto& out = pass == 0 ? engine : plain;
        out[0].resize (static_cast<std::size_t> (samples));
        out[1].resize (static_cast<std::size_t> (samples));
        stepper.processDouble (out[0].data(), out[1].data(), samples);
    }
    // The fine run: modal displacements at every substep.
    RuntimeOptions fineOptions = ro;
    fineOptions.oversample = 1;
    Bank fineBank = makeBank (*model, drum, mics, stroke, fineOptions);
    std::vector<BankMode> modes = fineBank.modes;
    // Bank index -> model index, as makeBank keeps modes in model order.
    std::vector<int> original;
    {
        std::size_t next = 0;
        for (std::size_t k = 0; k < model->modes.size() && next < modes.size(); ++k)
            if (model->modes[k].omega == modes[next].omega && model->modes[k].alpha == modes[next].alpha)
            {
                original.push_back (static_cast<int> (k));
                ++next;
            }
    }
    Stepper fineStepper (std::move (fineBank), fine, fineOptions);
    fineStepper.strike();
    const int steps = samples * os + 2;
    const std::size_t n = modes.size();
    std::vector<double> history (static_cast<std::size_t> (steps + 1) * n, 0.0); // index 0: before the strike
    for (int i = 1; i <= steps; ++i)
    {
        fineStepper.processDouble (nullptr, nullptr, 1);
        const auto& y = fineStepper.displacements();
        std::copy (y.begin(), y.end(), history.begin() + static_cast<std::ptrdiff_t> (static_cast<std::size_t> (i) * n));
    }
    const auto acceleration = [&] (std::size_t k, int i) // at time i * delta
    {
        if (i <= 0 || i >= steps)
            return 0.0;
        const auto at = [&] (int j) { return history[static_cast<std::size_t> (j) * n + k]; };
        return (at (i + 1) - 2.0 * at (i) + at (i - 1)) / (delta * delta);
    };
    Comparison result;
    double sumRef = 0.0, sumErr = 0.0, sumErrWithout = 0.0;
    const double thetaS = std::atan2 (stroke.strikeY, stroke.strikeX);
    for (int j = 0; j < 2; ++j)
    {
        const auto& p = mics.capsule[static_cast<std::size_t> (j)];
        std::vector<double> reference (static_cast<std::size_t> (samples), 0.0);
        for (std::size_t b = 0; b < n; ++b)
        {
            const CoupledMode& cm = model->modes[static_cast<std::size_t> (original[b])];
            std::vector<double> hist;
            const auto add = [&] (double weight, double dist)
            {
                const double x = dist / drum.soundSpeed / delta;
                const auto l = static_cast<std::size_t> (std::floor (x));
                if (hist.size() < l + 2)
                    hist.resize (l + 2, 0.0);
                hist[l] += weight * (1.0 - (x - static_cast<double> (l)));
                hist[l + 1] += weight * (x - static_cast<double> (l));
            };
            double b1Max = 0.0;
            const auto profile = batterProfile (cm, model->batter, drum.radius, b1Max);
            if (! cm.batterIndex.empty())
            {
                const double angular = std::cos (cm.m * (std::atan2 (p[1], p[0]) - thetaS));
                rayleighNodes (cm.m, drum.radius, profile, b1Max, std::hypot (cm.omega, cm.alpha) / drum.soundSpeed, std::hypot (p[0], p[1]), p[2], 1,
                               [&] (double w, double d) { add (angular * w, d); });
            }
            const double volume = rearVolume (cm, model->rear);
            if (volume != 0.0)
            {
                const double path = rearDistance (drum.depth, p);
                add (volume / (4.0 * pi * path), path);
            }
            for (int s = 0; s < samples; ++s)
            {
                const int now = (s + 1) * os;   // output sample s is at (s + 1) / 48000
                double sum = 0.0;
                for (std::size_t l = 0; l < hist.size() && static_cast<int> (l) < now; ++l)
                    sum += hist[l] * acceleration (b, now - static_cast<int> (l));
                reference[static_cast<std::size_t> (s)] += drum.airDensity * sum;
            }
        }
        double window[3] = { 0.0, 0.0, 0.0 }, windowRef[3] = { 0.0, 0.0, 0.0 };
        for (int s = 0; s < samples; ++s)
        {
            const double r = reference[static_cast<std::size_t> (s)];
            sumRef += r * r;
            result.peakReference = std::max (result.peakReference, std::abs (r));
            const double e = engine[j][static_cast<std::size_t> (s)] - r;
            sumErr += e * e;
            sumErrWithout += (plain[j][static_cast<std::size_t> (s)] - r) * (plain[j][static_cast<std::size_t> (s)] - r);
            const int w = s < 96 ? 0 : s < 288 ? 1 : 2;
            window[w] += e * e;
            windowRef[w] += r * r;
        }
        std::printf ("    capsule %d: rms error / rms reference in 0-2 / 2-6 / 6-10 ms: %.2e / %.2e / %.2e\n", j,
                     std::sqrt (window[0] / windowRef[0]), std::sqrt (window[1] / windowRef[1]), std::sqrt (window[2] / windowRef[2]));
    }
    result.rmsReference = std::sqrt (sumRef / (2.0 * samples));
    result.rmsError = std::sqrt (sumErr / (2.0 * samples));
    result.rmsErrorWithout = std::sqrt (sumErrWithout / (2.0 * samples));
    return result;
}
} // namespace

MODAL_TEST (output_retarded_reference)
{
    // First 10 ms, 4 kHz model, Berger off (the contact force is then the only
    // force besides each mode's own), and on (its own direct term is not in the
    // output; this measures what that leaves out).
    const Comparison linear = retardedReference (false, 480);
    std::printf ("  Berger off: reference rms %.3f Pa, peak %.3f Pa; error rms with the kernel %.4f Pa, without %.3f Pa\n",
                 linear.rmsReference, linear.peakReference, linear.rmsError, linear.rmsErrorWithout);
    modaltest::below ("Berger off: rms error / rms reference, first 10 ms", linear.rmsError / linear.rmsReference, 1e-2);
    const Comparison berger = retardedReference (true, 480);
    std::printf ("  Berger on: reference rms %.3f Pa; error rms with the kernel %.4f Pa, without %.3f Pa\n",
                 berger.rmsReference, berger.rmsError, berger.rmsErrorWithout);
    modaltest::note ("Berger on: rms error / rms reference (the Berger force's own direct term), first 10 ms: %.4f",
                     berger.rmsError / berger.rmsReference);
}

MODAL_TEST (output_causality)
{
    // Nothing can arrive before the nearest head point's delay z_mic / c. With
    // Berger off and the damping clusters dissolved, the contact force is the only
    // force besides each mode's own and the cancellation is exact; the clusters'
    // coupling forces and the Berger force carry no direct term of their own and
    // leave a small remainder.
    const auto mics = modaltest::nagadoMics();
    for (int variant = 0; variant < 3; ++variant)
    {
        RuntimeOptions with, without;
        with.berger = without.berger = variant == 2;
        without.feedthrough = false;
        const auto run = [&] (const RuntimeOptions& ro)
        {
            Bank bank = makeBank (fullModel(), modaltest::nagadoDrum(), mics, modaltest::nagadoStroke(), ro, 48000.0);
            if (variant == 0)
                bank.clusters.clear();
            Stepper stepper (std::move (bank), 48000.0, ro);
            stepper.strike();
            std::array<std::vector<double>, 2> out { std::vector<double> (480), std::vector<double> (480) };
            stepper.processDouble (out[0].data(), out[1].data(), 480);
            return out;
        };
        const auto a = run (with), b = run (without);
        double peak = 0.0;
        for (std::size_t i = 0; i < a[0].size(); ++i)
            peak = std::max ({ peak, std::abs (a[0][i]), std::abs (a[1][i]) });
        double early = 0.0, earlyWithout = 0.0;
        for (int j = 0; j < 2; ++j)
        {
            // Output sample i is the pressure at (i + 1) / fs after the strike.
            const double arrival = mics.capsule[static_cast<std::size_t> (j)][2] / modaltest::nagadoDrum().soundSpeed;
            for (std::size_t i = 0; (static_cast<double> (i) + 1.0) / 48000.0 < arrival; ++i)
            {
                early = std::max (early, std::abs (a[static_cast<std::size_t> (j)][i]));
                earlyWithout = std::max (earlyWithout, std::abs (b[static_cast<std::size_t> (j)][i]));
            }
        }
        const char* names[3] = { "no clusters, Berger off", "Berger off", "Berger on" };
        const double limits[3] = { 1e-6, 1e-3, 2e-3 };
        std::printf ("  %s: peak %.2f Pa; max |p| before the first arrival %.2e Pa (without the kernel %.3f Pa)\n",
                     names[variant], peak, early, earlyWithout);
        const std::string what = std::string (names[variant]) + ": max |p| before the first arrival / peak";
        modaltest::below (what.c_str(), early / peak, limits[variant]);
    }
}
