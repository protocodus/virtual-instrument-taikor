#pragma once

// A minimal registry for the finite-difference drum's verification suite.
// Every test records named checks with their measured and expected numbers;
// the runner (Main.cpp) prints them and exits non-zero if any check fails.
// Usage: TaikorCandidateTests_fdtd [substring ...] runs the tests whose names
// contain any of the substrings (all tests when none is given).

#include "ShippingSpec.h"

#include <cstdio>
#include <string>
#include <vector>

namespace taikor::candidates::fdtd::tests
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
};

std::vector<TestEntry>& registry();

struct Registrar
{
    Registrar (const char* name, TestFunction function) { registry().push_back ({ name, function }); }
};

// The nagado-daiko drum, Don stroke and capsules the render tool gives the
// engines (MIDI 60, velocity 127, 48 kHz), resolved once.
const ResolvedShipping& nagado();

std::string format (const char* pattern, ...);

} // namespace taikor::candidates::fdtd::tests

#define FDTD_TEST(function)                                                                  \
    static void function (taikor::candidates::fdtd::tests::Report&);                         \
    static const taikor::candidates::fdtd::tests::Registrar function##Registrar (#function,   \
                                                                                function);  \
    static void function (taikor::candidates::fdtd::tests::Report& report)
