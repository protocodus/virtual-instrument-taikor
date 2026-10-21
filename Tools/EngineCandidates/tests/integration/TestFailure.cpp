// Test 9: a rebuild at the stroke that fails is reported, not played. The
// drum rebuilds its heads at the stroke (for the stroke's energy, and the
// von Karman batter again at its sized spacing), and each rebuild calls the
// potential's prepare() for the new grid. If that fails, the head is left
// half-built, so the drum must stop: failure() says why, diagnostics() leads
// with it, and process() outputs silence, which the render tool rejects.
// The batter here is a Berger potential that refuses every prepare() after
// the first (the drum's own prepare), i.e. the rebuild for the Don's energy.

#include "IntegrationSupport.h"

#include <cmath>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::integration;

class RefusesRebuild final : public NonlinearPotential
{
public:
    std::string name() const override { return "Berger refusing a rebuild (test)"; }
    bool prepare (const MembraneGrid& grid, const HeadSpec& head, double radius, std::string& error) override
    {
        if (++calls > 1)
        {
            error = "prepare() refused by the test";
            return false;
        }
        return inner->prepare (grid, head, radius, error);
    }
    double evaluate (const double* w, double* gradient) override { return inner->evaluate (w, gradient); }

private:
    std::unique_ptr<NonlinearPotential> inner = makeBergerPotential();
    int calls = 0;
};
} // namespace

INTEGRATION_TEST (rebuildFailure)
{
    const auto& spec = nagado();
    FdtdOptions options;
    options.selfChecks = false;
    FdtdEngine<float> engine (std::make_unique<RefusesRebuild>(), options);
    std::string error;
    if (! engine.prepare (spec.drum, spec.mics, 48000.0, error))
        return report.check ("prepare", false, error);
    engine.strike (spec.stroke);
    std::vector<float> left (256, 1.0f), right (256, 1.0f);
    engine.process (left.data(), right.data(), 256);
    double peak = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i)
        peak = std::max ({ peak, std::abs (static_cast<double> (left[i])), std::abs (static_cast<double> (right[i])) });
    const std::string text = engine.diagnostics();
    report.check ("a failed rebuild at the stroke stops the drum and says why",
                  ! engine.failure().empty() && peak == 0.0 && engine.stepsTaken() == 0
                      && text.rfind ("- FAILED: ", 0) == 0,
                  format ("failure() \"%s\"; output peak %.1f after a 256-sample block; %lld steps taken; diagnostics "
                          "lead with the failure: %s",
                          engine.failure().c_str(), peak, engine.stepsTaken(),
                          text.rfind ("- FAILED: ", 0) == 0 ? "yes" : "NO"));
}
