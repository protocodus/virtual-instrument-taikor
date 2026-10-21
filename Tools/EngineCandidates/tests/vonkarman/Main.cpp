// Verification suite of the von Karman potential (engine id "vonkarman").
// Usage: TaikorCandidateTests_vonkarman [test ...]; with no arguments every
// test runs. Exit status is the number of failed checks (0 = all passed).

#include "TestUtil.h"

#include <cstring>
#include <functional>
#include <utility>

int main (int argc, char** argv)
{
    using namespace vktest;
    const std::vector<std::pair<const char*, std::function<void()>>> tests {
        { "gradient", testGradient },   // (1) finite-difference check of dV/dw
        { "scaling", testScaling },     // (2) V(s w) = s^4 V(w)
        { "hencky", testHencky },       // (3) Hencky's pressurised membrane
        { "analytic", testAnalytic },   // (4) exact continuum energies, grid convergence
        { "ratios", testRatios },       // (5) von Karman against Berger
        { "dent", testDent },           // (5b) local dent against its continuum value
        { "solver", testSolver },       // in-plane solver algebra and convergence
        { "sav", testSav },             // (6) SAV in-vacuo membrane: energy conservation
        { "resolution", testSavResolution }, // (6b) SAV tracking of the stroke's local stiffening
        { "tolerance", testTolerance }, // (6c) insensitivity to the solver tolerance
        { "benchmark", testBenchmark }, // (7) cost per evaluate at h = 6.9 and 4.8 mm
        { "profile", testProfile },     // cost breakdown of the solver kernels (report only)
    };
    for (const auto& [name, run] : tests)
    {
        bool selected = argc < 2 && std::strcmp (name, "profile") != 0;
        for (int i = 1; i < argc; ++i)
            selected = selected || std::strcmp (argv[i], name) == 0;
        if (! selected)
            continue;
        std::printf ("== %s\n", name);
        const double t0 = seconds();
        run();
        std::printf ("   (%.1f s)\n", seconds() - t0);
    }
    std::printf ("%s: %d failed check(s)\n", failures() == 0 ? "OK" : "FAILED", failures());
    return failures() == 0 ? 0 : 1;
}
