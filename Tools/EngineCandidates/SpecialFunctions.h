#pragma once

// Bessel functions of integer order and their zeros, for orders and arguments
// far past the shipping engine's ascending series (TaikoEngine::besselJ loses
// accuracy above x ~ 25; a 78 cm head needs j_mn ~ 90 at 5 kHz and ~ 200 at
// 20 kHz).
//
// J_m: Miller's backward recurrence J_{k-1} = (2k/x) J_k - J_{k+1}, normalised
// by J_0 + 2 sum_k J_2k = 1 (Abramowitz & Stegun 9.1.46 / 9.1.87; Numerical
// Recipes 3rd ed. sec. 6.5). It is stable for every order at once, so one call
// returns J_0..J_M.
// I_m: the same recurrence with I_{k-1} = (2k/x) I_k + I_{k+1}, normalised by
// e^x = I_0 + 2 sum_k I_k (A&S 9.6.37), returned scaled by e^-x so a stiff
// head's b2 a ~ 60 cannot overflow.
//
// The start index N = max(M, x) + 40 + 10 cbrt(max(M, x)) puts the discarded
// J_N below 1e-50 of the retained values for every x, M the candidates use.

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace taikor::candidates::special
{

inline constexpr double pi = 3.14159265358979323846;

[[nodiscard]] inline int millerStart (double x, int maxOrder) noexcept
{
    const double top = std::max (static_cast<double> (maxOrder), x);
    int start = static_cast<int> (std::ceil (top + 40.0 + 10.0 * std::cbrt (top)));
    return start + (start & 1);
}

// out[k] = J_k(x) for k = 0..maxOrder, x >= 0.
inline void besselJSequence (double x, int maxOrder, std::vector<double>& out)
{
    out.assign (static_cast<std::size_t> (maxOrder + 1), 0.0);
    if (x == 0.0)
    {
        out[0] = 1.0;
        return;
    }
    const int start = millerStart (x, maxOrder);
    double above = 0.0, current = 1.0e-300, sum = 0.0;
    for (int k = start; k >= 1; --k)
    {
        const double below = (2.0 * k / x) * current - above;
        above = current;
        current = below;              // current is now J_{k-1} (unnormalised)
        if (k - 1 <= maxOrder)
            out[static_cast<std::size_t> (k - 1)] = current;
        if (((k - 1) & 1) == 0 && k - 1 > 0)
            sum += 2.0 * current;
        if (std::abs (current) > 1.0e250)
        {
            current *= 1.0e-250;
            above *= 1.0e-250;
            sum *= 1.0e-250;
            for (auto& value : out)
                value *= 1.0e-250;
        }
    }
    sum += current;                   // J_0
    const double scale = 1.0 / sum;
    for (auto& value : out)
        value *= scale;
}

[[nodiscard]] inline double besselJ (int order, double x)
{
    if (order < 0)
        return (order & 1) ? -besselJ (-order, x) : besselJ (-order, x);
    if (x < 0.0)
        return (order & 1) ? -besselJ (order, -x) : besselJ (order, -x);
    std::vector<double> values;
    besselJSequence (x, order + 1, values);
    return values[static_cast<std::size_t> (order)];
}

// dJ_m/dx = (J_{m-1} - J_{m+1}) / 2, with J_{-1} = -J_1.
[[nodiscard]] inline double besselJPrime (int order, double x)
{
    std::vector<double> values;
    besselJSequence (std::abs (x), order + 2, values);
    const double below = order == 0 ? -values[1] : values[static_cast<std::size_t> (order - 1)];
    double result = 0.5 * (below - values[static_cast<std::size_t> (order + 1)]);
    if (x < 0.0 && (order & 1) == 0)
        result = -result;
    return result;
}

// out[k] = e^-x I_k(x) for k = 0..maxOrder, x >= 0.
inline void scaledBesselISequence (double x, int maxOrder, std::vector<double>& out)
{
    out.assign (static_cast<std::size_t> (maxOrder + 1), 0.0);
    if (x == 0.0)
    {
        out[0] = 1.0;
        return;
    }
    const int start = millerStart (x, maxOrder);
    double above = 0.0, current = 1.0e-300, sum = 0.0;
    for (int k = start; k >= 1; --k)
    {
        const double below = (2.0 * k / x) * current + above;
        above = current;
        current = below;
        if (k - 1 <= maxOrder)
            out[static_cast<std::size_t> (k - 1)] = current;
        if (k - 1 > 0)
            sum += 2.0 * current;
        if (current > 1.0e250)
        {
            current *= 1.0e-250;
            above *= 1.0e-250;
            sum *= 1.0e-250;
            for (auto& value : out)
                value *= 1.0e-250;
        }
    }
    sum += current;                   // I_0; sum is now e^x in the same scale
    const double scale = 1.0 / sum;
    for (auto& value : out)
        value *= scale;
}

[[nodiscard]] inline double scaledBesselI (int order, double x)
{
    std::vector<double> values;
    scaledBesselISequence (std::abs (x), std::abs (order) + 1, values);
    return values[static_cast<std::size_t> (std::abs (order))];
}

namespace detail
{
// Refines a bracketed root of f on [lo, hi] (f(lo), f(hi) of opposite sign)
// by bisection safeguarded secant steps to near machine precision.
template <typename F>
[[nodiscard]] double refineRoot (F&& f, double lo, double hi, double flo, double fhi)
{
    for (int iteration = 0; iteration < 200; ++iteration)
    {
        double mid = hi - fhi * (hi - lo) / (fhi - flo);
        if (! (mid > lo && mid < hi) || iteration % 3 == 2)
            mid = 0.5 * (lo + hi);
        const double fm = f (mid);
        if (fm == 0.0)
            return mid;
        if ((fm < 0.0) == (flo < 0.0))
        {
            lo = mid;
            flo = fm;
        }
        else
        {
            hi = mid;
            fhi = fm;
        }
        if (hi - lo < 1.0e-14 * std::max (1.0, std::abs (mid)))
            break;
    }
    return 0.5 * (lo + hi);
}
} // namespace detail

// Zeros of J_m for m = 0..maxOrder that lie below maxValue, ascending per order:
// result[m] = { j_m1, j_m2, ... }. A bracketing scan at a step of 0.2 cannot
// skip a zero: consecutive zeros of J_m are never closer than about 2.4.
[[nodiscard]] inline std::vector<std::vector<double>> besselJZeros (int maxOrder, double maxValue)
{
    std::vector<std::vector<double>> zeros (static_cast<std::size_t> (maxOrder + 1));
    constexpr double step = 0.2;
    std::vector<double> previous, current;
    double xPrevious = step * 0.5;
    besselJSequence (xPrevious, maxOrder, previous);
    for (double x = xPrevious + step; x <= maxValue + step; x += step)
    {
        besselJSequence (x, maxOrder, current);
        for (int m = 0; m <= maxOrder; ++m)
        {
            const double a = previous[static_cast<std::size_t> (m)];
            const double b = current[static_cast<std::size_t> (m)];
            if ((a < 0.0) != (b < 0.0) && a != 0.0)
            {
                const auto f = [m] (double t) { return besselJ (m, t); };
                const double root = detail::refineRoot (f, xPrevious, x, a, b);
                if (root < maxValue)
                    zeros[static_cast<std::size_t> (m)].push_back (root);
            }
        }
        std::swap (previous, current);
        xPrevious = x;
    }
    return zeros;
}

// Zeros of J_m' (the rigid-duct cross-mode eigenvalues alpha'_mq) below
// maxValue. For m = 0 the first entry is the plane-wave zero at 0.
[[nodiscard]] inline std::vector<std::vector<double>> besselJPrimeZeros (int maxOrder, double maxValue)
{
    std::vector<std::vector<double>> zeros (static_cast<std::size_t> (maxOrder + 1));
    zeros[0].push_back (0.0);
    constexpr double step = 0.2;
    std::vector<double> previous, current;
    const auto derivative = [] (const std::vector<double>& j, int m)
    {
        const double below = m == 0 ? -j[1] : j[static_cast<std::size_t> (m - 1)];
        return 0.5 * (below - j[static_cast<std::size_t> (m + 1)]);
    };
    double xPrevious = step * 0.5;
    besselJSequence (xPrevious, maxOrder + 1, previous);
    for (double x = xPrevious + step; x <= maxValue + step; x += step)
    {
        besselJSequence (x, maxOrder + 1, current);
        for (int m = 0; m <= maxOrder; ++m)
        {
            const double a = derivative (previous, m);
            const double b = derivative (current, m);
            if ((a < 0.0) != (b < 0.0) && a != 0.0)
            {
                const auto f = [m] (double t) { return besselJPrime (m, t); };
                const double root = detail::refineRoot (f, xPrevious, x, a, b);
                if (root < maxValue)
                    zeros[static_cast<std::size_t> (m)].push_back (root);
            }
        }
        std::swap (previous, current);
        xPrevious = x;
    }
    return zeros;
}

// Gauss-Legendre nodes and weights on [-1, 1] (Newton on P_n from Chebyshev
// guesses; Numerical Recipes 3rd ed. sec. 4.6).
inline void gaussLegendre (int count, std::vector<double>& nodes, std::vector<double>& weights)
{
    nodes.assign (static_cast<std::size_t> (count), 0.0);
    weights.assign (static_cast<std::size_t> (count), 0.0);
    const int half = (count + 1) / 2;
    for (int i = 0; i < half; ++i)
    {
        double z = std::cos (pi * (i + 0.75) / (count + 0.5));
        double derivative = 1.0;
        for (int iteration = 0; iteration < 100; ++iteration)
        {
            double p1 = 1.0, p2 = 0.0;
            for (int j = 0; j < count; ++j)
            {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
            }
            derivative = count * (z * p1 - p2) / (z * z - 1.0);
            const double previous = z;
            z = previous - p1 / derivative;
            if (std::abs (z - previous) < 1.0e-15)
                break;
        }
        nodes[static_cast<std::size_t> (i)] = -z;
        nodes[static_cast<std::size_t> (count - 1 - i)] = z;
        const double w = 2.0 / ((1.0 - z * z) * derivative * derivative);
        weights[static_cast<std::size_t> (i)] = w;
        weights[static_cast<std::size_t> (count - 1 - i)] = w;
    }
}

} // namespace taikor::candidates::special
