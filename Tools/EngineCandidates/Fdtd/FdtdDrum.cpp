// The finite-difference drum's factory (engine ids "berger" and, with the
// von Karman potential injected, "vonkarman"). See FdtdEngine.h.

#include "CandidateEngine.h"
#include "FdtdContracts.h"
#include "FdtdEngine.h"

#include <cstdlib>

namespace taikor::candidates
{
std::unique_ptr<Engine> makeFdtdEngine (std::unique_ptr<NonlinearPotential> batterNonlinearity)
{
    if (! batterNonlinearity)
        return nullptr;
    fdtd::FdtdOptions options;
    // Development only: TAIKOR_FDTD_SELF_CHECKS=0 skips the diagnostics' checks.
    if (const char* text = std::getenv ("TAIKOR_FDTD_SELF_CHECKS"))
        options.selfChecks = std::atoi (text) != 0;
    // Development only: TAIKOR_FDTD_WALL_LOSS=0 leaves the barrel's wall
    // lossless (the take before the wall loss, bit for bit).
    if (const char* text = std::getenv ("TAIKOR_FDTD_WALL_LOSS"))
        options.wallLoss = std::atoi (text) != 0;
    if (fdtd::FdtdNumerics::fromEnvironment().doubleAir)
        return std::make_unique<fdtd::FdtdEngine<double>> (std::move (batterNonlinearity), options);
    return std::make_unique<fdtd::FdtdEngine<float>> (std::move (batterNonlinearity), options);
}
} // namespace taikor::candidates
