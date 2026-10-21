// Verification suite of the complete coupled-mode drum (Modal/). Runs every
// registered test, or those whose names contain argv[1]; exits non-zero if any
// check fails.

#include "ModalTestUtil.h"

#include <chrono>
#include <cstring>

int main (int argc, char** argv)
{
    auto& registry = modaltest::Registry::get();
    int run = 0;
    for (auto& [name, body] : registry.tests)
    {
        if (argc > 1 && std::strstr (name.c_str(), argv[1]) == nullptr)
            continue;
        std::printf ("%s\n", name.c_str());
        const auto t0 = std::chrono::steady_clock::now();
        body();
        std::printf ("  (%.2f s)\n", std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count());
        ++run;
    }
    std::printf ("\n%d tests, %d checks, %d failures\n", run, registry.checks, registry.failures);
    return registry.failures == 0 ? 0 : 1;
}
