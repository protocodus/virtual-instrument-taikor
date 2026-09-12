#include "DSP/TaikoEngine.h"
#include "DSP/ShellBoundary.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>

namespace taikor
{
struct TaikoEngineTestAccess
{
    // Deliberately calculate one entry at a time, independently of the render
    // loop's combined-head traversal. Shared cavity coordinates add coherently
    // before squaring; the two orientations of a non-axisymmetric shape add
    // their squared contributions. Require identical float bits because even
    // tiny strain differences can accumulate into a different nonlinear tail.
    static std::array<float, TaikoEngine::modeEntryCount> reference (
        const TaikoEngine::Voice& voice, bool rear)
    {
        std::array<float, TaikoEngine::modeEntryCount> result {};
        for (int entry = 0; entry < TaikoEngine::legacyModeEntryCount; ++entry)
        {
            const auto slot = static_cast<std::size_t> (entry);
            if (entry < TaikoEngine::axisymmetricEntryCount)
            {
                double amplitude = 0.0;
                float norm = 0.0f;
                for (int index = 0; index < voice.activeModeCount; ++index)
                {
                    const auto& mode = voice.modes[static_cast<std::size_t> (index)];
                    if (! mode.membrane || mode.modeEntry >= TaikoEngine::legacyModeEntryCount)
                        continue;
                    if (mode.sharedCavityMode)
                    {
                        amplitude += mode.cavityBasis[(rear ? cavity::radialModeCount : 0) + slot]
                                   * mode.resonator.y1;
                        norm = mode.cavityStretchNorm[slot];
                    }
                    else if (mode.modeEntry == entry && mode.circumferentialOrder == 0)
                    {
                        amplitude += static_cast<double> (
                            rear ? mode.resonantParticipation : mode.batterParticipation)
                            * mode.resonator.y1;
                        norm = mode.stretchNorm;
                    }
                }
                result[slot] = static_cast<float> (
                    static_cast<double> (norm) * amplitude * amplitude);
            }
            else
            {
                for (int index = 0; index < voice.activeModeCount; ++index)
                {
                    const auto& mode = voice.modes[static_cast<std::size_t> (index)];
                    if (mode.membrane && mode.modeEntry == entry
                        && ! mode.sharedCavityMode && mode.circumferentialOrder != 0
                        && mode.rearHeadMode == rear)
                        result[slot] += static_cast<float> (
                            static_cast<double> (mode.stretchNorm)
                                * mode.resonator.y1 * mode.resonator.y1);
                }
            }
        }
        return result;
    }

    static bool verify()
    {
        auto engine = std::make_unique<TaikoEngine>();
        std::uint32_t random = 0x716f39abu;
        const auto next = [&random]()
        {
            random ^= random << 13;
            random ^= random >> 17;
            random ^= random << 5;
            return random;
        };
        int checked = 0;
        for (const bool shared : { false, true })
        {
            TaikoEngine::setRealismFeatures (TaikoEngine::allRealismFeatures
                & ~(shared ? 0u : static_cast<std::uint32_t> (TaikoEngine::sharedCavity)));
            for (const double sampleRate : { 8000.0, 48000.0, 192000.0 })
                for (int family = 0; family < drumCount; ++family)
                {
                    EngineParameters parameters;
                    parameters.tension = family % 2 == 0 ? 0.2f : 0.9f;
                    engine->setParameters (parameters);
                    engine->prepare (sampleRate, 64);
                    const auto drum = TaikoEngine::resolveDrumFor (parameters, 0.0f, family);
                    engine->ensurePhysicalDrum (family, drum);
                    auto& voice = engine->physicalDrums_[static_cast<std::size_t> (family)];
                    // Different last norms expose an accidental first-write
                    // shortcut, even though a physical radial norm is shared.
                    for (int index = 0; index < voice.modeCount; ++index)
                    {
                        auto& mode = voice.modes[static_cast<std::size_t> (index)];
                        const float scale = 1.0f + 0.001f * static_cast<float> (index);
                        mode.stretchNorm *= scale;
                        for (auto& norm : mode.cavityStretchNorm)
                            norm *= scale;
                    }
                    for (int order = 0; order < 2; ++order)
                    {
                        if (order != 0)
                            std::reverse (voice.modes.begin(), voice.modes.begin() + voice.modeCount);
                        // Production rebuilds these after sorting. This test
                        // deliberately changes the order outside the builder.
                        TaikoEngine::rebuildModeTraversals (voice);
                        for (const int active : { voice.modeCount, voice.modeCount / 2, 0 })
                        {
                            voice.activeModeCount = active;
                            for (int state = 0; state < 4; ++state)
                            {
                                for (int index = 0; index < active; ++index)
                                {
                                    auto& mode = voice.modes[static_cast<std::size_t> (index)];
                                    const double amplitude = static_cast<double> (
                                        static_cast<std::int32_t> (next())) / 2147483648.0;
                                    mode.resonator.y1 = std::ldexp (amplitude,
                                        static_cast<int> (next() % 24u) - 20);
                                }
                                std::array<float, TaikoEngine::modeEntryCount> batter {}, rear {};
                                TaikoEngine::membraneSquaredSlopePerEntry (voice, batter, rear);
                                const auto expectedBatter = reference (voice, false);
                                const auto expectedRear = reference (voice, true);
                                for (std::size_t entry = 0; entry < batter.size(); ++entry)
                                    if (std::bit_cast<std::uint32_t> (batter[entry])
                                            != std::bit_cast<std::uint32_t> (expectedBatter[entry])
                                        || std::bit_cast<std::uint32_t> (rear[entry])
                                            != std::bit_cast<std::uint32_t> (expectedRear[entry]))
                                    {
                                        std::cerr << "Strain mismatch: shared=" << shared
                                                  << " sampleRate=" << sampleRate
                                                  << " family=" << family << " order=" << order
                                                  << " active=" << active << " state=" << state
                                                  << " entry=" << entry << '\n';
                                        return false;
                                    }
                                ++checked;
                                if (! verifyBoundary (voice))
                                    return false;
                            }
                        }
                    }
                }
        }
        TaikoEngine::setRealismFeatures (TaikoEngine::allRealismFeatures);
        std::cout << "Parametric strain is bit-identical across " << checked
                  << " bank/state combinations\n";
        return true;
    }

    static bool verifyBoundary (TaikoEngine::Voice& voice)
    {
        std::array<double, TaikoEngine::shellResonatorCount> relativeVelocity {};
        std::array<double, TaikoEngine::resonatorCount> expected {};
        for (int index = 0; index < voice.activeModeCount; ++index)
        {
            const auto& mode = voice.modes[static_cast<std::size_t> (index)];
            if (mode.shellBoundaryImpulseScale == 0.0)
                continue;
            const double velocity = shellboundary::modalVelocity (
                mode.liveOmega, static_cast<double> (mode.decayRate + mode.appliedPalmDecay),
                mode.quadratureFromCurrent, mode.quadratureFromPrevious,
                mode.resonator.y1, mode.resonator.y2);
            relativeVelocity[static_cast<std::size_t> (mode.shellBoundaryGroup)] +=
                static_cast<double> (mode.shellBoundaryProjection) * velocity;
        }
        for (int index = 0; index < voice.modeCount; ++index)
        {
            const auto& mode = voice.modes[static_cast<std::size_t> (index)];
            auto& value = expected[static_cast<std::size_t> (index)];
            value = mode.resonator.y2;
            if (index < voice.activeModeCount && mode.shellBoundaryImpulseScale != 0.0)
            {
                const double change = mode.shellBoundaryImpulseScale
                    * relativeVelocity[static_cast<std::size_t> (mode.shellBoundaryGroup)];
                value += change * mode.shellBoundaryVelocityToPrevious;
            }
        }
        TaikoEngine::applyShellBoundary (voice);
        for (int index = 0; index < voice.modeCount; ++index)
            if (std::bit_cast<std::uint64_t> (voice.modes[static_cast<std::size_t> (index)].resonator.y2)
                != std::bit_cast<std::uint64_t> (expected[static_cast<std::size_t> (index)]))
            {
                std::cerr << "Sparse shell boundary differs from full-bank traversal\n";
                return false;
            }
        return true;
    }
};
} // namespace taikor

int main()
{
    return taikor::TaikoEngineTestAccess::verify() ? 0 : 1;
}
