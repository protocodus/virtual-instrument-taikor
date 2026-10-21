// Non-proportional damping of near-degenerate modes. A hide-damped head mode
// (alpha_h) and a lightly damped barrel air mode (alpha_a) with a small
// conservative coupling kappa form, in the conservative basis, the pair
// w0 +/- kappa with decay matrix [[abar, d], [d, abar]], abar = (alpha_h +
// alpha_a)/2, d = (alpha_h - alpha_a)/2. When d > kappa the true modes do not
// hybridise (damping-induced decoupling); first-order projection keeps only
// abar and gets both decays wrong. The dense cluster stepping must follow the
// exact damped two-degree-of-freedom solution.

#include "ModalRuntime.h"
#include "ModalTestUtil.h"

#include <cmath>

using namespace taikor::candidates;
using namespace taikor::candidates::modal;

namespace
{
// RK4 reference of y'' + 2 A y' + W^2 y = 0 for two coordinates.
std::vector<double> reference (double w1, double w2, const double a[4], double v1, double v2, double seconds)
{
    const double fs = 48000.0;
    const int sub = 64;
    const double h = 1.0 / (fs * sub);
    double s[4] = { 0.0, 0.0, v1, v2 };
    const auto f = [&] (const double* x, double* d)
    {
        d[0] = x[2];
        d[1] = x[3];
        d[2] = -w1 * w1 * x[0] - 2.0 * (a[0] * x[2] + a[1] * x[3]);
        d[3] = -w2 * w2 * x[1] - 2.0 * (a[2] * x[2] + a[3] * x[3]);
    };
    std::vector<double> out;
    const int samples = static_cast<int> (seconds * fs);
    for (int i = 0; i < samples; ++i)
    {
        for (int k = 0; k < sub; ++k)
        {
            double k1[4], k2[4], k3[4], k4[4], t[4];
            f (s, k1);
            for (int j = 0; j < 4; ++j) t[j] = s[j] + 0.5 * h * k1[j];
            f (t, k2);
            for (int j = 0; j < 4; ++j) t[j] = s[j] + 0.5 * h * k2[j];
            f (t, k3);
            for (int j = 0; j < 4; ++j) t[j] = s[j] + h * k3[j];
            f (t, k4);
            for (int j = 0; j < 4; ++j) s[j] += h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
        }
        out.push_back (s[0] + s[1]); // the head coordinate (x+ + x-)/sqrt2, times sqrt2
    }
    return out;
}

double relativeError (const std::vector<double>& x, const std::vector<double>& ref)
{
    double e = 0.0, r = 0.0;
    for (std::size_t i = 0; i < ref.size(); ++i)
    {
        e += (x[i] - ref[i]) * (x[i] - ref[i]);
        r += ref[i] * ref[i];
    }
    return std::sqrt (e / r);
}

std::vector<double> stepped (double w1, double w2, double abar, double d, bool cluster, double seconds)
{
    Bank bank;
    for (double w : { w1, w2 })
    {
        BankMode m;
        m.omega = w;
        m.alpha = abar;
        m.mic[0] = 1.0; // p = y of each mode; the sum is the head coordinate
        bank.modes.push_back (m);
    }
    if (cluster)
        bank.clusters.push_back ({ { 0, 1 }, { abar, d, d, abar } });
    RuntimeOptions ro;
    ro.berger = false;
    Stepper stepper (std::move (bank), 48000.0, ro);
    stepper.setModeState (0, 0.0, 1.0);
    stepper.setModeState (1, 0.0, 1.0);
    std::vector<double> out (static_cast<std::size_t> (seconds * 48000.0));
    stepper.processDouble (out.data(), nullptr, static_cast<int> (out.size()));
    return out;
}
} // namespace

MODAL_TEST (damping_cluster_pair)
{
    const double w0 = 2.0 * 3.14159265358979323846 * 3000.0, seconds = 0.15;
    for (const auto& [kappa, label] : { std::pair { 2.0, "damping-dominated (kappa 2 1/s)" },
                                        std::pair { 150.0, "coupling-dominated (kappa 150 1/s)" } })
    {
        const double ah = 75.0, aa = 1.5, abar = 0.5 * (ah + aa), d = 0.5 * (ah - aa);
        const double a[4] = { abar, d, d, abar };
        const auto ref = reference (w0 + kappa, w0 - kappa, a, 1.0, 1.0, seconds);
        const double eCluster = relativeError (stepped (w0 + kappa, w0 - kappa, abar, d, true, seconds), ref);
        const double eDiagonal = relativeError (stepped (w0 + kappa, w0 - kappa, abar, d, false, seconds), ref);
        std::printf ("  %s: rms error dense cluster %.2e, proportional projection %.2e\n", label, eCluster, eDiagonal);
        modaltest::below ("dense cluster vs exact damped 2-DOF solution (rms / rms), 150 ms", eCluster, 2e-3);
    }
}
