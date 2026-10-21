#include "ModalMath.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace taikor::candidates::modal
{
namespace
{
// J_m and J_m' at x from one Miller sequence.
void besselJAndPrime (int m, double x, double& j, double& jPrime)
{
    thread_local std::vector<double> seq;
    special::besselJSequence (std::abs (x), m + 1, seq);
    j = seq[static_cast<std::size_t> (m)];
    const double below = m == 0 ? -seq[1] : seq[static_cast<std::size_t> (m - 1)];
    jPrime = 0.5 * (below - seq[static_cast<std::size_t> (m + 1)]);
}

// e^{-x} I_m(x) and I_m'(x)/I_m(x).
void scaledIAndRatio (int m, double x, double& iScaled, double& ratio)
{
    thread_local std::vector<double> seq;
    special::scaledBesselISequence (x, m + 1, seq);
    iScaled = seq[static_cast<std::size_t> (m)];
    const double below = m == 0 ? seq[1] : seq[static_cast<std::size_t> (m - 1)];
    ratio = 0.5 * (below + seq[static_cast<std::size_t> (m + 1)]) / iScaled;
}

double b2Of (double b1, const HeadMaterial& h) { return std::sqrt (b1 * b1 + h.tension / h.bending); }

double omegaOf (double b1, const HeadMaterial& h)
{
    return std::sqrt ((h.tension * b1 * b1 + h.bending * b1 * b1 * b1 * b1) / h.arealDensity);
}

// The clamped frequency function divided by I_m(b2 a).
double clampedFunction (int m, double x, const HeadMaterial& h)
{
    const double b1 = x / h.radius;
    const double b2 = b2Of (b1, h);
    double j = 0.0, jp = 0.0, is = 0.0, ratio = 0.0;
    besselJAndPrime (m, x, j, jp);
    scaledIAndRatio (m, b2 * h.radius, is, ratio);
    return b1 * jp - b2 * j * ratio;
}
} // namespace

ShapeValue evaluateShape (const HeadMode& mode, double radius, double r)
{
    const int m = mode.m;
    double j = 0.0, jp = 0.0;
    besselJAndPrime (m, mode.b1 * r, j, jp);
    thread_local std::vector<double> seq;
    special::scaledBesselISequence (mode.b2 * r, m + 1, seq);
    const double below = m == 0 ? seq[1] : seq[static_cast<std::size_t> (m - 1)];
    const double iPrimeScaled = 0.5 * (below + seq[static_cast<std::size_t> (m + 1)]);
    const double scale = mode.jEdge / mode.iEdgeScaled * std::exp (mode.b2 * (r - radius));
    ShapeValue out;
    out.value = j - scale * seq[static_cast<std::size_t> (m)];
    out.slope = mode.b1 * jp - scale * mode.b2 * iPrimeScaled;
    return out;
}

std::vector<std::vector<HeadMode>> solveHeadModes (const HeadMaterial& h, double omegaMax)
{
    if (h.tension <= 0.0 || h.bending <= 0.0 || h.arealDensity <= 0.0 || h.radius <= 0.0)
        throw std::invalid_argument ("head material must be positive");
    const double a = h.radius;
    const double w2 = omegaMax * omegaMax;
    const double b1Max = std::sqrt ((std::sqrt (h.tension * h.tension + 4.0 * h.bending * h.arealDensity * w2)
                                     - h.tension) / (2.0 * h.bending));
    const double xMax = b1Max * a;
    const int maxOrder = static_cast<int> (std::ceil (xMax)) + 2;

    std::vector<std::vector<HeadMode>> result;
    // Scan all orders at once: one J sequence at x and one scaled I sequence at b2 a per step.
    constexpr double step = 0.05;
    std::vector<double> jPrev, jCur, iSeq;
    std::vector<double> fPrev (static_cast<std::size_t> (maxOrder + 1)), fCur (fPrev.size());
    const auto evaluateAll = [&] (double x, std::vector<double>& jSeq, std::vector<double>& f)
    {
        const double b1 = x / a;
        const double b2 = b2Of (b1, h);
        special::besselJSequence (x, maxOrder + 1, jSeq);
        special::scaledBesselISequence (b2 * a, maxOrder + 1, iSeq);
        for (int m = 0; m <= maxOrder; ++m)
        {
            const auto k = static_cast<std::size_t> (m);
            const double jpm = 0.5 * ((m == 0 ? -jSeq[1] : jSeq[k - 1]) - jSeq[k + 1]);
            const double ratio = 0.5 * ((m == 0 ? iSeq[1] : iSeq[k - 1]) + iSeq[k + 1]) / iSeq[k];
            f[k] = b1 * jpm - b2 * jSeq[k] * ratio;
        }
    };
    double xPrev = step;
    evaluateAll (xPrev, jPrev, fPrev);
    std::vector<std::vector<double>> roots (static_cast<std::size_t> (maxOrder + 1));
    for (double x = xPrev + step; x <= xMax + step; x += step)
    {
        evaluateAll (x, jCur, fCur);
        for (int m = 0; m <= maxOrder; ++m)
        {
            const auto k = static_cast<std::size_t> (m);
            // No clamped root lies below j'_m1 > m: there J_m, J_m' > 0 and
            // b1 J_m'/J_m < m/a < b2 I_m'/I_m, so the function is negative. Below
            // x = m, J_m(x) also underflows towards denormals whose sign is noise, so
            // the scan for order m starts there (tests/modal heads_root_scan).
            if (xPrev < m)
                continue;
            if ((fPrev[k] < 0.0) != (fCur[k] < 0.0) && fPrev[k] != 0.0)
            {
                const auto f = [m, &h] (double t) { return clampedFunction (m, t, h); };
                roots[k].push_back (special::detail::refineRoot (f, xPrev, x, fPrev[k], fCur[k]));
            }
        }
        std::swap (fPrev, fCur);
        std::swap (jPrev, jCur);
        xPrev = x;
    }

    for (int m = 0; m <= maxOrder; ++m)
    {
        std::vector<HeadMode> modes;
        int n = 0;
        for (double x : roots[static_cast<std::size_t> (m)])
        {
            HeadMode mode;
            mode.m = m;
            mode.n = ++n;
            mode.b1 = x / a;
            mode.b2 = b2Of (mode.b1, h);
            mode.omega = omegaOf (mode.b1, h);
            if (mode.omega > omegaMax)
                break;
            double jp = 0.0;
            besselJAndPrime (m, x, mode.jEdge, jp);
            mode.q = mode.b1 * jp;
            scaledIAndRatio (m, mode.b2 * a, mode.iEdgeScaled, mode.iRatio);
            // int (J - cI)^2 r dr = int J^2 + c^2 int I^2: the cross term
            // a [b2 J_m(b1 a) I_m'(b2 a) - b1 J_m'(b1 a) I_m(b2 a)] / (b1^2 + b2^2)
            // vanishes by the clamped edge condition. Lommel's squared integrals:
            // int J_m(b1 r)^2 r dr = a^2/2 [J_m'^2 + (1 - m^2/x^2) J_m^2],
            // int I_m(b2 r)^2 r dr = a^2/2 [(1 + m^2/y^2) I_m^2 - I_m'^2].
            const double y = mode.b2 * a;
            const double jTerm = jp * jp + (1.0 - double (m * m) / (x * x)) * mode.jEdge * mode.jEdge;
            const double iTerm = mode.jEdge * mode.jEdge * ((1.0 + double (m * m) / (y * y)) - mode.iRatio * mode.iRatio);
            mode.normSq = angularNorm (m) * 0.5 * a * a * (jTerm + iTerm);
            mode.gradDiag = mode.b1 * mode.b1 * mode.normSq;
            if (m == 0)
            {
                // int J_0(b1 r) r dr = a J_1(b1 a)/b1; int I_0(b2 r) r dr = a I_1(b2 a)/b2.
                std::vector<double> iseq;
                special::scaledBesselISequence (y, 1, iseq);
                const double is1 = iseq[1] / iseq[0];
                mode.volume = 2.0 * pi * a * (special::besselJ (1, x) / mode.b1 - mode.jEdge * is1 / mode.b2);
            }
            if (! (mode.normSq > 0.0))
            {
                --n; // not a mode (a sign change of noise); never expected after the scan start above
                continue;
            }
            modes.push_back (mode);
        }
        if (modes.empty())
            break; // the lowest root rises with m, so no higher order has one either
        result.push_back (std::move (modes));
    }
    return result;
}

double hankelTransform (const HeadMode& mode, double radius, double kappa, double jm, double jmPrime)
{
    const double b1s = mode.b1 * mode.b1;
    const double k2 = kappa * kappa;
    const double delta = 1.0e-6 * mode.b1;
    if (std::abs (kappa - mode.b1) < delta)
    {
        // The removable singularity: the exact limit at kappa = b1 (L'Hopital
        // with Bessel's equation) is W(b1) = a^2/2 [J_m'^2 + (1 - m^2/x^2) J_m^2],
        // x = b1 a; within delta of it the transform is the quadratic through
        // W(b1 - delta), W(b1), W(b1 + delta), whose own error is O((delta a)^3).
        const double x = mode.b1 * radius;
        const double jp = mode.q / mode.b1;
        const double w0 = 0.5 * radius * radius
                        * (jp * jp + (1.0 - double (mode.m * mode.m) / (x * x)) * mode.jEdge * mode.jEdge);
        double side[2];
        for (int s = 0; s < 2; ++s)
        {
            const double k = mode.b1 + (s == 0 ? -2.0 : 2.0) * delta;
            double j = 0.0, jpk = 0.0;
            besselJAndPrime (mode.m, k * radius, j, jpk);
            side[s] = hankelTransform (mode, radius, k, j, jpk);
        }
        const double t = (kappa - mode.b1) / (2.0 * delta);
        return w0 + 0.5 * (side[1] - side[0]) * t + 0.5 * (side[1] - 2.0 * w0 + side[0]) * t * t;
    }
    const double b2s = mode.b2 * mode.b2;
    const double numerator = radius * (b1s + b2s) * (kappa * jmPrime * mode.jEdge - mode.q * jm);
    return numerator / ((b1s - k2) * (b2s + k2));
}

double diskAverage (const HeadMode& mode, double radius, double rs, double b)
{
    const double x1 = mode.b1 * b, x2 = mode.b2 * b;
    const double jMean = x1 > 0.0 ? 2.0 * special::besselJ (1, x1) / x1 : 1.0;
    thread_local std::vector<double> seq;
    special::scaledBesselISequence (x2, 1, seq);
    const double iMean = x2 > 0.0 ? 2.0 * seq[1] * std::exp (x2) / x2 : 1.0;
    special::scaledBesselISequence (mode.b2 * rs, mode.m, seq);
    const double iPart = mode.jEdge * seq[static_cast<std::size_t> (mode.m)] / mode.iEdgeScaled
                       * std::exp (mode.b2 * (rs - radius));
    return special::besselJ (mode.m, mode.b1 * rs) * jMean - iPart * iMean;
}

GradientForm gradientForm (const std::vector<HeadMode>& modes, const HeadMaterial& h, double relativeTolerance)
{
    GradientForm form;
    const std::size_t n = modes.size();
    form.diagonal.resize (n);
    if (n == 0)
        return form;
    const int m = modes.front().m;
    const double a = h.radius;
    const double lambda = angularNorm (m);
    const double tOverD = h.tension / h.bending;
    // E_ij = (c_i A_i)(c_j A_j) P_ij / (b1_i^2 + b1_j^2 + T/D), with c_i c_j P_ij in
    // closed form from the modified-Bessel Lommel integral divided by I_m(b2_i a) I_m(b2_j a).
    const auto entry = [&] (std::size_t i, std::size_t j)
    {
        const HeadMode& p = modes[i];
        const HeadMode& r = modes[j];
        double pij = 0.0;
        if (i == j)
        {
            const double y = p.b2 * a;
            pij = lambda * 0.5 * a * a * ((1.0 + double (m * m) / (y * y)) - p.iRatio * p.iRatio);
        }
        else
        {
            pij = lambda * a * (p.b2 * p.iRatio - r.b2 * r.iRatio) / (p.b2 * p.b2 - r.b2 * r.b2);
        }
        const double ai = p.b1 * p.b1 + p.b2 * p.b2, aj = r.b1 * r.b1 + r.b2 * r.b2;
        return p.jEdge * r.jEdge * ai * aj * pij / (p.b1 * p.b1 + r.b1 * r.b1 + tOverD);
    };
    std::vector<double> residual (n);
    for (std::size_t i = 0; i < n; ++i)
    {
        form.diagonal[i] = modes[i].gradDiag;
        residual[i] = entry (i, i);
    }
    std::vector<char> used (n, 0);
    for (std::size_t rank = 0; rank < n; ++rank)
    {
        std::size_t pivot = 0;
        double best = -1.0;
        for (std::size_t i = 0; i < n; ++i)
            if (! used[i] && residual[i] / form.diagonal[i] > best)
            {
                best = residual[i] / form.diagonal[i];
                pivot = i;
            }
        if (best <= relativeTolerance)
        {
            form.truncation = std::max (best, 0.0);
            break;
        }
        std::vector<double> column (n);
        const double root = std::sqrt (residual[pivot]);
        for (std::size_t i = 0; i < n; ++i)
        {
            double value = entry (i, pivot);
            for (const auto& l : form.lowRank)
                value -= l[i] * l[pivot];
            column[i] = used[i] ? 0.0 : value / root;
        }
        column[pivot] = root;
        used[pivot] = 1;
        form.pivots.push_back (best);
        for (std::size_t i = 0; i < n; ++i)
            residual[i] = used[i] ? 0.0 : residual[i] - column[i] * column[i];
        form.lowRank.push_back (std::move (column));
    }
    return form;
}

} // namespace taikor::candidates::modal
