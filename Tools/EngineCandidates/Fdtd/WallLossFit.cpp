// Weights of the diffusive representation of sqrt(s) (WallLoss.h), fitted in
// the discrete domain the scheme actually realises.
//
// A branch (d_t. + xi mu_t.) theta = p responds at angular frequency w as
//   B(w) = j sin(wk)/k / (j sin(wk)/k + xi cos(wk)),
// the trapezoidal rule over 2k of s / (s + xi): the analogue section at the
// warped frequency tan(wk)/k (10 % high at 8 kHz at 96 kHz), so the weights are
// fitted to the discrete B_i directly rather than taken from a quadrature of
// the analogue integral. The targets follow from the decay rate of a discrete
// air mode, alpha = P / 2H with the leapfrog energy H: for the thermal part
// the power ~ |rho d_t.Psi|^2 and H carry the same factor sin^2(wk)/k^2, so the
// target is sqrt(j w) itself; the viscous power ~ |G|^2 while H ~
// cos^2(wk/2) sum |G|^2, so its target is cos^2(wk/2) sqrt(j w). A mode at w
// then decays at the Kirchhoff rate of w.
//
// Weights are non-negative (Lawson-Hanson NNLS), which keeps every branch
// positive real whatever the fit's error. The fit minimises the relative error
// of the real part (the loss) and, weighted by 0.1, of the imaginary part (the
// boundary layer's reactance: it shifts the barrel's modes by the loss rate
// over w, < 1e-4 relative above 100 Hz, so its accuracy is immaterial).
//
// Numerical parameters, all justified by tests/fdtd wallSquareRootFit: 12
// poles log-spaced over 1 Hz - 300 kHz (more than a decade beyond the band at
// each end), a 20 Hz - 12 kHz band (past the air grid's on-axis cut-off,
// 8.87 kHz; above it the branches stay passive, only less exact) on 300
// log-spaced points. At 96 kHz the realised loss is then within 0.06 % /
// 0.34 % (thermal / viscous) of the target over 20 Hz - 1 kHz and 0.08 % /
// 1.08 % over 1 - 9 kHz. The target itself cannot be met near the top of the
// grid's band: the discrete branches see the frequency warped to tan(wk)/k,
// and a positive sum of first-order sections cannot keep a 45-degree phase
// while its magnitude flattens, so a 16 kHz band leaves 15-40 % errors.
// With the same band 10 poles give 0.30 / 0.36 % and 0.45 / 0.75 %, 8 poles
// 1.2 / 1.2 % and 1.5 / 3.1 % (scratch scan of the pole count and range).

#include "WallLoss.h"

#include <algorithm>
#include <cmath>

namespace taikor::candidates::fdtd
{
namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr int poleCount = 12;
constexpr double poleLow = 1.0, poleHigh = 3.0e5, bandLow = 20.0, bandHigh = 12000.0;
constexpr int bandPoints = 300;
constexpr double errorHigh = 9000.0; // reported accuracy: up to the air grid's on-axis cut-off
constexpr double reactiveWeight = 0.1;

// Least squares min |A x - b| over the columns in `use` (Householder QR); the
// other entries of x are zero. A is m x n, column-major.
std::vector<double> leastSquares (const std::vector<double>& a, const std::vector<double>& b, int m, int n,
                                  const std::vector<int>& use)
{
    const int p = static_cast<int> (use.size());
    std::vector<double> q (static_cast<std::size_t> (m * p)), r = b;
    for (int j = 0; j < p; ++j)
        for (int i = 0; i < m; ++i)
            q[static_cast<std::size_t> (j * m + i)] = a[static_cast<std::size_t> (use[static_cast<std::size_t> (j)] * m + i)];
    for (int j = 0; j < p; ++j)
    {
        double norm = 0.0;
        for (int i = j; i < m; ++i)
            norm += q[static_cast<std::size_t> (j * m + i)] * q[static_cast<std::size_t> (j * m + i)];
        norm = std::sqrt (norm);
        double& top = q[static_cast<std::size_t> (j * m + j)];
        const double alpha = top > 0.0 ? -norm : norm;
        std::vector<double> v (static_cast<std::size_t> (m - j));
        for (int i = j; i < m; ++i)
            v[static_cast<std::size_t> (i - j)] = q[static_cast<std::size_t> (j * m + i)];
        v[0] -= alpha;
        double vv = 0.0;
        for (const double e : v)
            vv += e * e;
        if (vv <= 0.0)
            continue;
        const auto reflect = [&] (double* column)
        {
            double dot = 0.0;
            for (int i = j; i < m; ++i)
                dot += v[static_cast<std::size_t> (i - j)] * column[i];
            const double f = 2.0 * dot / vv;
            for (int i = j; i < m; ++i)
                column[i] -= f * v[static_cast<std::size_t> (i - j)];
        };
        for (int l = j; l < p; ++l)
            reflect (q.data() + l * m);
        reflect (r.data());
    }
    std::vector<double> x (static_cast<std::size_t> (n), 0.0), z (static_cast<std::size_t> (p), 0.0);
    for (int j = p - 1; j >= 0; --j)
    {
        double s = r[static_cast<std::size_t> (j)];
        for (int l = j + 1; l < p; ++l)
            s -= q[static_cast<std::size_t> (l * m + j)] * z[static_cast<std::size_t> (l)];
        z[static_cast<std::size_t> (j)] = s / q[static_cast<std::size_t> (j * m + j)];
    }
    for (int j = 0; j < p; ++j)
        x[static_cast<std::size_t> (use[static_cast<std::size_t> (j)])] = z[static_cast<std::size_t> (j)];
    return x;
}

// Lawson & Hanson, Solving Least Squares Problems (SIAM 1995), ch. 23, NNLS.
std::vector<double> nonNegativeLeastSquares (const std::vector<double>& a, const std::vector<double>& b, int m, int n)
{
    std::vector<double> x (static_cast<std::size_t> (n), 0.0);
    std::vector<bool> passive (static_cast<std::size_t> (n), false);
    const auto gradient = [&] (std::vector<double>& w)
    {
        std::vector<double> residual = b;
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i)
                residual[static_cast<std::size_t> (i)] -= a[static_cast<std::size_t> (j * m + i)] * x[static_cast<std::size_t> (j)];
        w.assign (static_cast<std::size_t> (n), 0.0);
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i)
                w[static_cast<std::size_t> (j)] += a[static_cast<std::size_t> (j * m + i)] * residual[static_cast<std::size_t> (i)];
    };
    std::vector<double> w;
    for (int outer = 0; outer < 3 * n; ++outer)
    {
        gradient (w);
        int best = -1;
        for (int j = 0; j < n; ++j)
            if (! passive[static_cast<std::size_t> (j)] && w[static_cast<std::size_t> (j)] > 1e-12
                && (best < 0 || w[static_cast<std::size_t> (j)] > w[static_cast<std::size_t> (best)]))
                best = j;
        if (best < 0)
            break;
        passive[static_cast<std::size_t> (best)] = true;
        for (int inner = 0; inner < 3 * n; ++inner)
        {
            std::vector<int> use;
            for (int j = 0; j < n; ++j)
                if (passive[static_cast<std::size_t> (j)])
                    use.push_back (j);
            const auto z = leastSquares (a, b, m, n, use);
            bool feasible = true;
            for (const int j : use)
                feasible = feasible && z[static_cast<std::size_t> (j)] > 0.0;
            if (feasible)
            {
                x = z;
                break;
            }
            double step = 1.0;
            for (const int j : use)
                if (z[static_cast<std::size_t> (j)] <= 0.0)
                    step = std::min (step, x[static_cast<std::size_t> (j)] / (x[static_cast<std::size_t> (j)] - z[static_cast<std::size_t> (j)]));
            for (int j = 0; j < n; ++j)
            {
                x[static_cast<std::size_t> (j)] += step * (z[static_cast<std::size_t> (j)] - x[static_cast<std::size_t> (j)]);
                if (passive[static_cast<std::size_t> (j)] && x[static_cast<std::size_t> (j)] <= 1e-15)
                {
                    passive[static_cast<std::size_t> (j)] = false;
                    x[static_cast<std::size_t> (j)] = 0.0;
                }
            }
        }
    }
    return x;
}

std::complex<double> branch (double xi, double k, double w)
{
    const std::complex<double> s (0.0, std::sin (w * k) / k);
    return s / (s + xi * std::cos (w * k));
}

std::complex<double> target (double k, double w, bool viscous)
{
    const std::complex<double> root = std::sqrt (std::complex<double> (0.0, w));
    const double half = std::cos (0.5 * w * k);
    return viscous ? root * half * half : root;
}
} // namespace

std::complex<double> realisedSquareRoot (const SquareRootFit& fit, const std::vector<double>& weights, double k, double w)
{
    std::complex<double> sum = 0.0;
    for (std::size_t i = 0; i < fit.xi.size(); ++i)
        sum += weights[i] * branch (fit.xi[i], k, w);
    return sum;
}

SquareRootFit fitSquareRoot (double k)
{
    SquareRootFit fit;
    fit.bandLow = bandLow;
    fit.bandHigh = bandHigh;
    fit.errorHigh = errorHigh;
    for (int i = 0; i < poleCount; ++i)
        fit.xi.push_back (2.0 * pi * poleLow * std::pow (poleHigh / poleLow, i / (poleCount - 1.0)));
    const int m = 2 * bandPoints, n = poleCount;
    std::vector<double> frequencies;
    for (int f = 0; f < bandPoints; ++f)
        frequencies.push_back (bandLow * std::pow (bandHigh / bandLow, f / (bandPoints - 1.0)));
    for (const bool viscous : { false, true })
    {
        std::vector<double> a (static_cast<std::size_t> (m * n)), b (static_cast<std::size_t> (m));
        for (int f = 0; f < bandPoints; ++f)
        {
            const double w = 2.0 * pi * frequencies[static_cast<std::size_t> (f)];
            const auto t = target (k, w, viscous);
            const double scale = 1.0 / std::abs (t);
            b[static_cast<std::size_t> (2 * f)] = t.real() * scale;
            b[static_cast<std::size_t> (2 * f + 1)] = reactiveWeight * t.imag() * scale;
            for (int j = 0; j < n; ++j)
            {
                const auto e = branch (fit.xi[static_cast<std::size_t> (j)], k, w);
                a[static_cast<std::size_t> (j * m + 2 * f)] = e.real() * scale;
                a[static_cast<std::size_t> (j * m + 2 * f + 1)] = reactiveWeight * e.imag() * scale;
            }
        }
        auto& weights = viscous ? fit.viscous : fit.thermal;
        weights = nonNegativeLeastSquares (a, b, m, n);
        double worst = 0.0;
        for (const double f : frequencies)
        {
            if (f > errorHigh)
                break;
            const double w = 2.0 * pi * f;
            worst = std::max (worst, std::abs (realisedSquareRoot (fit, weights, k, w).real() / target (k, w, viscous).real() - 1.0));
        }
        (viscous ? fit.viscousError : fit.thermalError) = worst;
    }
    return fit;
}

} // namespace taikor::candidates::fdtd
