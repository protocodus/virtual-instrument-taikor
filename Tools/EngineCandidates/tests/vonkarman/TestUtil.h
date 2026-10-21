#pragma once

// Shared helpers of the von Karman verification suite.

#include "ShippingSpec.h"
#include "VonKarmanPotential.h"

#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <functional>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

namespace vktest
{
using namespace taikor::candidates;

inline int& failures()
{
    static int count = 0;
    return count;
}

inline std::string format (const char* f, ...)
{
    char buffer[1024];
    va_list args;
    va_start (args, f);
    std::vsnprintf (buffer, sizeof (buffer), f, args);
    va_end (args);
    return buffer;
}

inline void check (bool ok, const std::string& what)
{
    std::printf ("  [%s] %s\n", ok ? "pass" : "FAIL", what.c_str());
    std::fflush (stdout);
    if (! ok)
        ++failures();
}

inline void note (const std::string& what)
{
    std::printf ("         %s\n", what.c_str());
    std::fflush (stdout);
}

// The drum and stroke the harness resolves for the listening-test stroke
// (nagado-daiko Don, MIDI 60, velocity 127): the suite tests the potential
// on the batter it will actually run with.
inline const ResolvedShipping& shipping()
{
    static const ResolvedShipping resolved = []
    {
        ShippingStroke stroke;
        stroke.parameters = listeningTestParameters();
        stroke.midiNote = 60;
        stroke.velocity = 1.0f;
        return resolveShipping (stroke);
    }();
    return resolved;
}

inline const HeadSpec& batter() { return shipping().drum.batter; }
inline double headRadius() { return shipping().drum.radius; }
inline double c11()
{
    const auto& b = batter();
    return b.youngsModulus * b.thickness / (1.0 - b.poissonRatio * b.poissonRatio);
}

struct Fixture
{
    MembraneGrid grid;
    std::unique_ptr<VonKarmanPotential> vk;
    std::string error;
};

// A potential on the staircase disk of radius a with node spacing h. The
// tension-relative stopping rule is off unless asked for: the static tests
// check the exact discrete potential, the dynamic ones the shipping setup.
inline Fixture makeFixture (double a, double h, double tolerance, bool tensionRelative = false)
{
    Fixture f;
    f.grid = makeDiskGrid (a, h);
    f.vk = std::make_unique<VonKarmanPotential>();
    f.vk->solverSettings().relativeTolerance = tolerance;
    f.vk->setTensionRelativeStopping (tensionRelative);
    if (! f.vk->prepare (f.grid, batter(), a, f.error))
        f.vk.reset();
    return f;
}

// Discrete int |grad w|^2 dA over every edge (the Berger letter's measure).
inline double gradientSquared (const MembraneGrid& g, const double* w)
{
    double s = 0.0;
    for (int j = 0; j + 1 < g.ny; ++j)
        for (int i = 0; i + 1 < g.nx; ++i)
        {
            const double dx = w[g.index (i + 1, j)] - w[g.index (i, j)];
            const double dy = w[g.index (i, j + 1)] - w[g.index (i, j)];
            s += dx * dx + dy * dy;
        }
    return s; // (slope h)^2 summed = h^2 sum slope^2
}

// Exact Berger energy with an immovable edge, V_B = C11 / (8 pi a^2) I^2
// (Berger, J. Appl. Mech. 22, 1955; FdtdContracts.h).
inline double bergerEnergy (const MembraneGrid& g, const double* w, double a)
{
    const double I = gradientSquared (g, w);
    return c11() / (8.0 * std::numbers::pi * a * a) * I * I;
}

// Fills w with f(x, y) at interior nodes and zero elsewhere.
template <typename F>
std::vector<double> sampleField (const MembraneGrid& g, F&& f)
{
    std::vector<double> w (static_cast<std::size_t> (g.size()), 0.0);
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i)
            if (g.interior[static_cast<std::size_t> (g.index (i, j))] != 0)
                w[static_cast<std::size_t> (g.index (i, j))] = f (g.x (i), g.y (j));
    return w;
}

// Solves -Lap_h w = f (5-point, w = 0 off the interior) by conjugate
// gradients to a relative residual of 1e-13; f and w are grid-sized.
inline std::vector<double> solvePoisson (const MembraneGrid& g, const std::vector<double>& f)
{
    const auto size = static_cast<std::size_t> (g.size());
    auto A = [&g] (const std::vector<double>& x, std::vector<double>& y)
    {
        const double ih2 = 1.0 / (g.h * g.h);
        for (int j = 1; j + 1 < g.ny; ++j)
            for (int i = 1; i + 1 < g.nx; ++i)
            {
                const int k = g.index (i, j);
                y[static_cast<std::size_t> (k)] = g.interior[static_cast<std::size_t> (k)] == 0 ? 0.0
                    : ih2 * (4.0 * x[static_cast<std::size_t> (k)] - x[static_cast<std::size_t> (k - 1)] - x[static_cast<std::size_t> (k + 1)]
                             - x[static_cast<std::size_t> (k - g.nx)] - x[static_cast<std::size_t> (k + g.nx)]);
            }
    };
    std::vector<double> x (size, 0.0), r (size, 0.0), p (size, 0.0), q (size, 0.0);
    for (std::size_t k = 0; k < size; ++k)
        r[k] = g.interior[k] != 0 ? f[k] : 0.0;
    p = r;
    double rr = 0.0, bb = 0.0;
    for (double v : r)
        rr += v * v;
    bb = rr;
    for (int it = 0; it < 20000 && rr > 1e-26 * bb; ++it)
    {
        A (p, q);
        double pq = 0.0;
        for (std::size_t k = 0; k < size; ++k)
            pq += p[k] * q[k];
        const double alpha = rr / pq;
        double rrNew = 0.0;
        for (std::size_t k = 0; k < size; ++k)
        {
            x[k] += alpha * p[k];
            r[k] -= alpha * q[k];
            rrNew += r[k] * r[k];
        }
        for (std::size_t k = 0; k < size; ++k)
            p[k] = r[k] + rrNew / rr * p[k];
        rr = rrNew;
    }
    return x;
}

// Continuum stretching energy and exact Berger energy of an axisymmetric
// field from its slope w'(r) (TestAnalytic.cpp).
std::pair<double, double> axisymmetricEnergy (const std::function<double (double)>& slope, double a, double nu, double C);

// The same for w = f(r) cos(m theta) by Fourier decomposition and linear
// finite elements over `elements` intervals (TestFourier.cpp).
std::pair<double, double> fourierEnergy (const std::function<double (double)>& f, const std::function<double (double)>& slope,
                                         int m, double a, double nu, double C, int elements);

inline double seconds()
{
    using clock = std::chrono::steady_clock;
    return std::chrono::duration<double> (clock::now().time_since_epoch()).count();
}

// CPU time of the calling thread: the benchmark's measure, because the
// machine is shared and wall time includes time the thread was descheduled.
inline double cpuSeconds()
{
    timespec t {};
    clock_gettime (CLOCK_THREAD_CPUTIME_ID, &t);
    return static_cast<double> (t.tv_sec) + 1e-9 * static_cast<double> (t.tv_nsec);
}

// Test entry points (one file each).
void testGradient();
void testScaling();
void testAnalytic();
void testHencky();
void testRatios();
void testDent();
void testSav();
void testTolerance();
void testBenchmark();
void testProfile();
void testSolver();
void testSavResolution();

} // namespace vktest
