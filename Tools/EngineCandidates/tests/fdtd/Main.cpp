// Verification suite of the finite-difference drum (engine id "berger").
// See TestSupport.h for the command line.

#include "TestSupport.h"

#include <chrono>
#include <cstdarg>
#include <cstring>

namespace taikor::candidates::fdtd::tests
{

std::vector<TestEntry>& registry()
{
    static std::vector<TestEntry> entries;
    return entries;
}

const ResolvedShipping& nagado()
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

std::string format (const char* pattern, ...)
{
    char buffer[2048];
    va_list args;
    va_start (args, pattern);
    std::vsnprintf (buffer, sizeof (buffer), pattern, args);
    va_end (args);
    return buffer;
}

} // namespace taikor::candidates::fdtd::tests

int main (int argc, char** argv)
{
    using namespace taikor::candidates::fdtd::tests;
    int failed = 0, passed = 0;
    for (const auto& entry : registry())
    {
        bool selected = argc < 2;
        for (int a = 1; a < argc; ++a)
            selected = selected || std::strstr (entry.name, argv[a]) != nullptr;
        if (! selected)
            continue;
        std::printf ("%s\n", entry.name);
        std::fflush (stdout);
        Report report;
        const auto t0 = std::chrono::steady_clock::now();
        entry.function (report);
        const double seconds = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
        for (const auto& check : report.checks)
            (check.passed ? passed : failed) += 1;
        std::printf ("  (%.1f s)\n", seconds);
    }
    std::printf ("%d checks passed, %d failed\n", passed, failed);
    return failed == 0 && passed > 0 ? 0 : 1;
}
