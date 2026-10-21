#pragma once

// A minimal registry for the coupled-mode drum's verification suite. Every
// test prints what it measured against what it expected; the suite exits
// non-zero if any check fails.

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace modaltest
{

struct Registry
{
    std::vector<std::pair<std::string, std::function<void()>>> tests;
    int failures = 0, checks = 0;
    static Registry& get()
    {
        static Registry registry;
        return registry;
    }
};

struct Register
{
    Register (const char* name, std::function<void()> body)
    {
        Registry::get().tests.emplace_back (name, std::move (body));
    }
};

inline void check (bool ok, const char* what, double measured, double expected, double tolerance)
{
    auto& r = Registry::get();
    ++r.checks;
    if (! ok)
        ++r.failures;
    std::printf ("  [%s] %s: measured %.10g, expected %.10g (tol %.3g)\n", ok ? " ok " : "FAIL", what,
                 measured, expected, tolerance);
}

inline void near (const char* what, double measured, double expected, double absTolerance)
{
    check (std::abs (measured - expected) <= absTolerance && std::isfinite (measured), what, measured, expected,
           absTolerance);
}

inline void nearRel (const char* what, double measured, double expected, double relTolerance)
{
    check (std::abs (measured - expected) <= relTolerance * std::abs (expected) && std::isfinite (measured), what,
           measured, expected, relTolerance);
}

inline void below (const char* what, double measured, double limit)
{
    check (measured <= limit && std::isfinite (measured), what, measured, limit, 0.0);
}

inline void note (const char* format, double value)
{
    std::printf ("  (info) ");
    std::printf (format, value);
    std::printf ("\n");
}

} // namespace modaltest

#define MODAL_TEST(name)                                                      \
    static void name();                                                       \
    static modaltest::Register name##_registration (#name, name);             \
    static void name()
