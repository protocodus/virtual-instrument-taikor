#pragma once

// Verification suite of the von Karman letter (engine id "vonkarman"): the
// finite-difference drum of Fdtd/ with the von Karman potential of VonKarman/
// on its batter. The two parts have their own suites (tests/fdtd,
// tests/vonkarman); this one checks what only the combination can show.
// Usage: TaikorCandidateTests_integration [substring ...] runs the tests whose
// names contain any of the substrings (every default test when none is
// given; the heavy ones run only when named in full). Exit status non-zero if
// any check fails.

#include "ShippingSpec.h"

#include "FdtdEngine.h"

#include <array>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace taikor::candidates::integration
{

struct Check
{
    std::string name;
    bool passed = false;
    std::string detail;
};

class Report
{
public:
    void check (const std::string& name, bool passed, const std::string& detail)
    {
        checks.push_back ({ name, passed, detail });
        std::printf ("  [%s] %s: %s\n", passed ? "PASS" : "FAIL", name.c_str(), detail.c_str());
        std::fflush (stdout);
    }
    void note (const std::string& text)
    {
        std::printf ("  %s\n", text.c_str());
        std::fflush (stdout);
    }
    std::vector<Check> checks;
};

using TestFunction = void (*) (Report&);

struct TestEntry
{
    const char* name;
    TestFunction function;
    bool heavy;
};

std::vector<TestEntry>& registry();

struct Registrar
{
    Registrar (const char* name, TestFunction function, bool heavy) { registry().push_back ({ name, function, heavy }); }
};

// The nagado-daiko drum, Don stroke and capsules the render tool gives the
// engines (MIDI 60, velocity 127, 48 kHz), resolved once.
const ResolvedShipping& nagado();

std::string format (const char* pattern, ...);

// A run of the drum from its stroke: the raw capsule pressures at the internal
// rate (Pa, left then right per step) and the engine for inspection.
template <typename Real>
struct Take
{
    std::unique_ptr<fdtd::FdtdEngine<Real>> engine;
    std::vector<double> left, right;
};

// Builds, prepares and strikes an engine with the given batter potential and
// options, then advances `steps` internal steps recording the capsules.
template <typename Real>
Take<Real> runTake (std::unique_ptr<NonlinearPotential> batter, const fdtd::FdtdOptions& options,
                    const StrokeSpec& stroke, int steps, std::string& error);

// Energy of a signal over [begin, end) and of the difference of two.
double energy (const std::vector<double>& a, std::size_t begin, std::size_t end);
double differenceEnergy (const std::vector<double>& a, const std::vector<double>& b, std::size_t begin, std::size_t end);

// 10 log10 of the difference energy over the reference energy.
double differenceDb (const std::vector<double>& reference, const std::vector<double>& other, std::size_t begin,
                     std::size_t end);

} // namespace taikor::candidates::integration

#define INTEGRATION_TEST_IMPL(function, heavy)                                                      \
    static void function (taikor::candidates::integration::Report&);                               \
    static const taikor::candidates::integration::Registrar function##Registrar (#function,          \
                                                                                  function, heavy); \
    static void function (taikor::candidates::integration::Report& report)
#define INTEGRATION_TEST(function) INTEGRATION_TEST_IMPL (function, false)
#define INTEGRATION_HEAVY_TEST(function) INTEGRATION_TEST_IMPL (function, true)
