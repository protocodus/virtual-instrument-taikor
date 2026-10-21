// (1) The returned gradient is the derivative of the returned energy, and
// (2) the energy is exactly quartic in the amplitude.

#include "TestUtil.h"

#include <algorithm>

namespace vktest
{

namespace
{
// A field exercising every stencil entry: a (0,1)-like bowl, an asymmetric
// (1,1)-like term, a dent at the strike point and small node-to-node noise.
std::vector<double> testField (const MembraneGrid& g, double a, double amplitude, unsigned seed)
{
    std::uint64_t state = 0x853c49e6748fea9bULL + seed;
    auto rnd = [&state]
    {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double> (state >> 11) * 1.1102230246251565e-16 - 0.5;
    };
    const auto& s = shipping().stroke;
    return sampleField (g, [&] (double x, double y)
    {
        const double r2 = (x * x + y * y) / (a * a);
        const double bowl = 1.0 - r2;
        const double dx = x - s.strikeX, dy = y - s.strikeY;
        const double dent = std::exp (-(dx * dx + dy * dy) / (0.03 * 0.03));
        return amplitude * (0.6 * bowl + 0.3 * bowl * bowl * x / a + 0.4 * dent * bowl + 0.1 * rnd());
    });
}

// Five-point central difference: exact for a quartic in t, and V(w + t d) is
// a quartic in t (V = q^T S q / 2 with S fixed and q quadratic in w).
double directional (VonKarmanPotential& vk, const std::vector<double>& w,
                    const std::vector<double>& d, double eps)
{
    std::vector<double> p (w.size());
    auto at = [&] (double t)
    {
        for (std::size_t k = 0; k < w.size(); ++k)
            p[k] = w[k] + t * d[k];
        return vk.evaluate (p.data(), nullptr);
    };
    return (-at (2 * eps) + 8 * at (eps) - 8 * at (-eps) + at (-2 * eps)) / (12 * eps);
}
} // namespace

void testGradient()
{
    const double a = headRadius();
    Fixture f = makeFixture (a, a / 20.0, 1.0e-13);
    if (! f.vk)
        return check (false, "prepare: " + f.error);
    auto& vk = *f.vk;
    vk.setWarmStartOrder (0);
    const double amplitude = 0.005;
    const auto w = testField (f.grid, a, amplitude, 1);
    std::vector<double> grad (w.size());
    const double V = vk.evaluate (w.data(), grad.data());
    note (format ("grid n = %d, h = %.2f mm, V = %.6e J, PCG %d its, residual %.1e", f.grid.nx,
                  f.grid.h * 1e3, V, vk.stats().lastIterations, vk.stats().lastResidual));

    // Random directions.
    double worst = 0.0;
    std::uint64_t state = 12345;
    for (int trial = 0; trial < 6; ++trial)
    {
        std::vector<double> d (w.size(), 0.0);
        for (int k = 0; k < f.grid.size(); ++k)
        {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            if (f.grid.interior[static_cast<std::size_t> (k)] != 0)
                d[static_cast<std::size_t> (k)] = static_cast<double> (state >> 11) * 1.1102230246251565e-16 - 0.5;
        }
        const double fd = directional (vk, w, d, 0.2 * amplitude);
        double an = 0.0;
        for (std::size_t k = 0; k < w.size(); ++k)
            an += grad[k] * d[k];
        worst = std::max (worst, std::abs (fd - an) / std::abs (an));
    }
    check (worst < 1.0e-9, format ("random directions: max |FD - grad.d| / |grad.d| = %.2e (< 1e-9)", worst));

    // Single nodes: centre, strike point, and nodes on the rim staircase.
    double worstNode = 0.0, scale = 0.0;
    for (double g : grad)
        scale = std::max (scale, std::abs (g));
    int tested = 0;
    for (int j = 0; j < f.grid.ny; ++j)
        for (int i = 0; i < f.grid.nx; ++i)
        {
            const int k = f.grid.index (i, j);
            if (f.grid.interior[static_cast<std::size_t> (k)] == 0)
                continue;
            const bool rim = f.grid.interior[static_cast<std::size_t> (k + 1)] == 0
                          || f.grid.interior[static_cast<std::size_t> (k - 1)] == 0
                          || f.grid.interior[static_cast<std::size_t> (k + f.grid.nx)] == 0
                          || f.grid.interior[static_cast<std::size_t> (k - f.grid.nx)] == 0;
            if (! rim && (i + 3 * j) % 17 != 0)
                continue;
            std::vector<double> d (w.size(), 0.0);
            d[static_cast<std::size_t> (k)] = 1.0;
            const double fd = directional (vk, w, d, 0.2 * amplitude);
            worstNode = std::max (worstNode, std::abs (fd - grad[static_cast<std::size_t> (k)]) / scale);
            ++tested;
        }
    check (worstNode < 1.0e-9, format ("%d single nodes (all rim nodes + a lattice): max |FD - dV/dw_l| / max|dV/dw| = %.2e (< 1e-9)",
                                       tested, worstNode));
    double clampedMax = 0.0;
    for (int k = 0; k < f.grid.size(); ++k)
        if (f.grid.interior[static_cast<std::size_t> (k)] == 0)
            clampedMax = std::max (clampedMax, std::abs (grad[static_cast<std::size_t> (k)]));
    check (clampedMax == 0.0, "gradient is exactly zero at clamped nodes");
}

void testScaling()
{
    const double a = headRadius();
    Fixture f = makeFixture (a, a / 28.0, 1.0e-10);
    if (! f.vk)
        return check (false, "prepare: " + f.error);
    auto& vk = *f.vk;
    vk.setWarmStartOrder (0); // cold starts: PCG iterates scale exactly with b
    const auto w = testField (f.grid, a, 0.004, 2);
    std::vector<double> g1 (w.size()), gs (w.size()), ws (w.size());
    const double V1 = vk.evaluate (w.data(), g1.data());
    double worstV = 0.0, worstG = 0.0;
    for (double s : { 1.0e-3, 0.5, 2.0, 10.0 })
    {
        for (std::size_t k = 0; k < w.size(); ++k)
            ws[k] = s * w[k];
        const double Vs = vk.evaluate (ws.data(), gs.data());
        worstV = std::max (worstV, std::abs (Vs / (s * s * s * s * V1) - 1.0));
        double num = 0.0, den = 0.0;
        for (std::size_t k = 0; k < w.size(); ++k)
        {
            num = std::max (num, std::abs (gs[k] - s * s * s * g1[k]));
            den = std::max (den, std::abs (s * s * s * g1[k]));
        }
        worstG = std::max (worstG, num / den);
    }
    check (worstV < 1.0e-12, format ("V(s w) / (s^4 V(w)) - 1 over s = 1e-3, 0.5, 2, 10: max %.2e (< 1e-12)", worstV));
    check (worstG < 1.0e-12, format ("grad V(s w) against s^3 grad V(w): max rel %.2e (< 1e-12)", worstG));
    std::vector<double> zero (w.size(), 0.0), gz (w.size(), 1.0);
    const double V0 = vk.evaluate (zero.data(), gz.data());
    const double gzMax = *std::max_element (gz.begin(), gz.end(), [] (double x, double y) { return std::abs (x) < std::abs (y); });
    check (V0 == 0.0 && gzMax == 0.0, "V(0) = 0 and grad V(0) = 0 exactly");
}

} // namespace vktest
