#include "DSP/TaikoEngine.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

namespace
{
int failures = 0;
void expect (bool condition, const char* message)
{
    if (! condition)
    {
        if (failures < 20)
            std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}
bool near (double a, double b, double tolerance = 3.0e-6)
{
    return std::abs (a - b) <= tolerance * std::max ({ 1.0, std::abs (a), std::abs (b) });
}
}

namespace taikor
{
struct BachiModalIntegrationTestAccess
{
    static void run()
    {
        TaikoEngine engine;
        engine.prepare (48000.0, 256);
        auto voice = std::make_unique<TaikoEngine::Voice>();
        auto otherVoice = std::make_unique<TaikoEngine::Voice>();
        const auto& stroke = TaikoEngine::strikeProfile (Articulation::Don);
        constexpr float strikeRadius = 0.31f;

        const auto findMode = [] (const TaikoEngine::Voice& bank, int id)
            -> const TaikoEngine::Mode*
        {
            for (int i = 0; i < bank.modeCount; ++i)
                if (bank.modes[static_cast<std::size_t> (i)].physicalIndex == id)
                    return &bank.modes[static_cast<std::size_t> (i)];
            return nullptr;
        };

        for (int family = 0; family < drumCount; ++family)
        {
            const auto drum = TaikoEngine::resolveDrumFor (EngineParameters {}, 0.0f, family);
            voice->strikeRadius = strikeRadius;
            voice->strikeAngle = 0.0f;
            engine.buildVoiceModes (*voice, drum, stroke, 0.0f, false);

            // Includes both original and appended radial orders. In particular
            // entries 20 and 25 are (1,4)/(2,4), not the fourth global row.
            for (const auto [entryIndex, radialOrder] :
                    { std::pair { 4, 1 }, { 5, 2 }, { 7, 1 },
                      { 20, 4 }, { 25, 4 }, { 32, 5 } })
            {
                const auto& entry = TaikoEngine::membraneModes()[static_cast<std::size_t> (entryIndex)];
                const auto basis = TaikoEngine::basisForEntry (drum, entryIndex);
                const auto expectedBasis = physical::basisForMode (
                    drum.physicalFamily, drum.physicalFamilyMix,
                    entry.circumferentialOrder, radialOrder);
                for (int branch = 0; branch < 2; ++branch)
                {
                    expect (basis.rotation[static_cast<std::size_t> (branch)]
                                == expectedBasis.rotation[static_cast<std::size_t> (branch)],
                            "appended membrane entries must retain their physical radial order");
                    const auto* mode = findMode (*voice, 2 * entryIndex + branch);
                    expect (mode != nullptr, "all selected principal modes must be in the live builder");
                    if (mode == nullptr)
                        continue;
                    const auto observed = TaikoEngine::observeMode (
                        drum, entryIndex, branch, strikeRadius, 0.0f);
                    expect (near (mode->omega / 6.2831853071795864769, observed.frequencyHz),
                            "readout and renderer must use identical principal-mode frequencies");
                    expect (mode->angularBasis == basis.rotation[static_cast<std::size_t> (branch)],
                            "live modes must retain the physical basis needed by state remapping");
                    if (entryIndex < TaikoEngine::legacyModeEntryCount)
                    {
                        const double radial = TaikoEngine::besselJ (
                            entry.circumferentialOrder, entry.besselZero * strikeRadius);
                        const double patch = bachi::diskProjection (entry.besselZero
                            * bachi::contactRadius (drum.bachi, drum.radius, strikeRadius)
                            / drum.radius);
                        const double projection = radial * patch * basis.project (0.0f, branch)
                                                * stroke.membraneGain * stroke.levelScale;
                        expect (near (mode->contactShape, projection, 2.0e-7),
                                "live contact must use the finite footprint in its principal basis");
                    }
                }
            }

            constexpr int entryIndex = 4; // first m=1 pair, mechanically coupled
            const auto basis = TaikoEngine::basisForEntry (drum, entryIndex);
            const auto& row = basis.rotation[0];
            const float nodalAngle = std::atan2 (row[0], -row[1]);
            auto nodalDrum = drum;
            nodalDrum.micAngleLeft = nodalAngle;
            otherVoice->strikeRadius = strikeRadius;
            otherVoice->strikeAngle = nodalAngle;
            engine.buildVoiceModes (*otherVoice, nodalDrum, stroke, 0.0f, true);
            const auto* nodal = findMode (*otherVoice, 2 * entryIndex);
            const auto* orthogonal = findMode (*otherVoice, 2 * entryIndex + 1);
            expect (nodal != nullptr && orthogonal != nullptr,
                    "both principal modes must survive a nodal gesture");
            if (nodal != nullptr && orthogonal != nullptr)
            {
                expect (std::abs (nodal->contactShape) < 2.0e-7f,
                        "striking a principal nodal diameter must not drive that mode");
                expect (std::abs (orthogonal->contactShape) > 1.0e-4f,
                        "a nodal gesture must still drive the orthogonal member");
                const double nodalPressure = std::hypot (
                    nodal->micLeft, nodal->micLeftQuadrature);
                const double orthogonalPressure = std::hypot (
                    orthogonal->micLeft, orthogonal->micLeftQuadrature);
                expect (nodalPressure <= 1.0e-6 * std::max (orthogonalPressure, 1.0e-12),
                        "both pressure quadratures must share the microphone nodal diameter");
            }

            const auto angleZeroSecond = TaikoEngine::observeMode (
                drum, entryIndex, 1, strikeRadius, 0.0f);
            expect (angleZeroSecond.frequencyHz > 0.0f && angleZeroSecond.amplitude > 0.0f,
                    "authored angle zero must not suppress a rotated head's second branch");

            // Repositioning a finite collocated contact changes both its force
            // drive and displacement sensing by exactly the same factor.
            otherVoice->strikeAngle = 0.43f;
            engine.buildVoiceModes (*otherVoice, drum, stroke, 0.0f, false);
            for (int branch = 0; branch < 2; ++branch)
            {
                const auto* first = findMode (*voice, 2 * entryIndex + branch);
                const auto* second = findMode (*otherVoice, 2 * entryIndex + branch);
                if (first != nullptr && second != nullptr)
                    expect (near (static_cast<double> (first->drive) * second->contactShape,
                                  static_cast<double> (second->drive) * first->contactShape,
                                  2.0e-9),
                            "force and sensing projections must remain reciprocal as contact moves");
            }
        }
    }
};
}

int main()
{
    taikor::BachiModalIntegrationTestAccess::run();
    if (failures != 0)
        return 1;
    std::cout << "Live bachi footprint, principal axes and readout consistency passed.\n";
    return 0;
}
