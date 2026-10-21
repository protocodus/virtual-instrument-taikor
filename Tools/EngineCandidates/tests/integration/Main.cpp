// Verification suite of the von Karman letter; see IntegrationSupport.h.

#include "IntegrationSupport.h"

#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstring>

namespace taikor::candidates::integration
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
    char buffer[4096];
    va_list args;
    va_start (args, pattern);
    std::vsnprintf (buffer, sizeof (buffer), pattern, args);
    va_end (args);
    return buffer;
}

template <typename Real>
Take<Real> runTake (std::unique_ptr<NonlinearPotential> batter, const fdtd::FdtdOptions& options,
                    const StrokeSpec& stroke, int steps, std::string& error)
{
    Take<Real> take;
    take.engine = std::make_unique<fdtd::FdtdEngine<Real>> (std::move (batter), options);
    if (! take.engine->prepare (nagado().drum, nagado().mics, 48000.0, error))
        return take;
    take.engine->strike (stroke);
    take.left.reserve (static_cast<std::size_t> (steps));
    take.right.reserve (static_cast<std::size_t> (steps));
    for (int s = 0; s < steps; ++s)
    {
        take.engine->advance (1);
        take.left.push_back (take.engine->micPressure (0));
        take.right.push_back (take.engine->micPressure (1));
    }
    return take;
}

template Take<float> runTake<float> (std::unique_ptr<NonlinearPotential>, const fdtd::FdtdOptions&, const StrokeSpec&,
                                     int, std::string&);
template Take<double> runTake<double> (std::unique_ptr<NonlinearPotential>, const fdtd::FdtdOptions&,
                                       const StrokeSpec&, int, std::string&);

double energy (const std::vector<double>& a, std::size_t begin, std::size_t end)
{
    double sum = 0.0;
    for (std::size_t i = begin; i < end && i < a.size(); ++i)
        sum += a[i] * a[i];
    return sum;
}

double differenceEnergy (const std::vector<double>& a, const std::vector<double>& b, std::size_t begin, std::size_t end)
{
    double sum = 0.0;
    for (std::size_t i = begin; i < end && i < a.size() && i < b.size(); ++i)
        sum += (a[i] - b[i]) * (a[i] - b[i]);
    return sum;
}

double differenceDb (const std::vector<double>& reference, const std::vector<double>& other, std::size_t begin,
                     std::size_t end)
{
    return 10.0 * std::log10 (std::max (1e-300, differenceEnergy (reference, other, begin, end))
                              / std::max (1e-300, energy (reference, begin, end)));
}

} // namespace taikor::candidates::integration

int main (int argc, char** argv)
{
    using namespace taikor::candidates::integration;
    int failed = 0, passed = 0;
    for (const auto& entry : registry())
    {
        bool selected = argc < 2 && ! entry.heavy;
        for (int a = 1; a < argc; ++a)
            selected = selected
                       || (entry.heavy ? std::strcmp (entry.name, argv[a]) == 0 : std::strstr (entry.name, argv[a]) != nullptr);
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
