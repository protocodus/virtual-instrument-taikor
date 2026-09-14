#include "DSP/TaikoEngine.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <memory>

namespace taikor
{
struct TaikoEngineTestAccess
{
    static bool check()
    {
        auto engine = std::make_unique<TaikoEngine>();
        std::uint32_t random = 0x58acb124u;
        const auto next = [&random]()
        {
            random ^= random << 13; random ^= random >> 17; random ^= random << 5;
            return static_cast<float> (static_cast<std::int32_t> (random)) / 2147483648.0f;
        };
        int checked = 0;
        for (double rate : { 8000.0, 48000.0, 96000.0, 192000.0, 384000.0 })
        for (int drum = 0; drum < drumCount; ++drum)
        {
            EngineParameters parameters;
            engine->setParameters (parameters);
            engine->prepare (rate, 64);
            engine->ensurePhysicalDrum (drum, TaikoEngine::resolveDrumFor (parameters, 0.0f, drum));
            auto& voice = engine->physicalDrums_[static_cast<std::size_t> (drum)];
            for (int index = 0; index < voice.modeCount; ++index)
            {
                auto& mode = voice.modes[static_cast<std::size_t> (index)];
                mode.shellRingIndex = static_cast<std::int8_t> (index % TaikoEngine::shellResonatorCount);
                mode.shellRimCoupling = 0.125f;
                mode.rimDrive = 0.25f;
            }
            // Shared modes use radial weights regardless of their entry ID.
            for (auto& mode : voice.modes)
                if (mode.sharedCavityMode)
                {
                    mode.modeEntry = TaikoEngine::legacyModeEntryCount + 1;
                    break;
                }
            TaikoEngine::rebuildModeTraversals (voice);
            for (bool parametric : { false, true })
            for (bool rim : { false, true })
            for (int active : { voice.modeCount, voice.modeCount - 1,
                                voice.modeCount / 2, 3, 1, 0 })
            for (bool zero : { false, true })
            {
                voice.activeModeCount = active;
                std::array<float, TaikoEngine::modeEntryCount> ripple {}, rear {};
                // Upper entries are zero, exactly as in renderVoice().
                for (int entry = 0; entry < TaikoEngine::legacyModeEntryCount; ++entry)
                {
                    ripple[static_cast<std::size_t> (entry)] = zero ? 0.0f : next();
                    rear[static_cast<std::size_t> (entry)] = zero ? 0.0f : next();
                }
                for (auto& input : voice.modalInput) input = next() * 0.01f;
                for (auto& force : voice.shellRimForce) force = next();
                const auto previousRim = voice.shellRimForce;
                auto expectedRim = voice.shellRimForce;
                if (rim) expectedRim.fill (0.0f);
                auto expectedInput = voice.modalInput;
                std::array<double, TaikoEngine::resonatorCount> expectedY1 {}, expectedY2 {};
                const double scale = 0.25 / (rate * rate);
                for (int index = 0; index < voice.modeCount; ++index)
                {
                    auto& mode = voice.modes[static_cast<std::size_t> (index)];
                    auto& r = mode.resonator;
                    r.y1 = next() * 0.01; r.y2 = next() * 0.01;
                    double result = r.y1;
                    if (index < active)
                    {
                        if (rim && mode.shellRingIndex >= 0)
                        {
                            const auto ring = static_cast<std::size_t> (mode.shellRingIndex);
                            if (mode.membrane)
                                expectedRim[ring] += mode.shellRimCoupling * static_cast<float> (r.y1);
                            else
                                expectedInput[mode.physicalIndex] += mode.rimDrive * previousRim[ring];
                        }
                        float amount = parametric && mode.membrane
                            ? ripple[mode.modeEntry] * mode.batterTensionFraction
                                + rear[mode.modeEntry] * mode.rearTensionFraction : 0.0f;
                        if (parametric && mode.sharedCavityMode)
                        {
                            amount = 0.0f;
                            for (std::size_t radial = 0; radial < cavity::radialModeCount; ++radial)
                                amount += mode.cavityTensionWeights[radial] * ripple[radial]
                                        + mode.cavityRearTensionWeights[radial] * rear[radial];
                        }
                        result = r.b0 * static_cast<double> (expectedInput[mode.physicalIndex])
                               - r.a1 * r.y1 - r.a2 * r.y2;
                        if (amount != 0.0f)
                        {
                            const double increment = scale * static_cast<double> (amount)
                                * mode.liveOmega * mode.liveOmega * (1.0 + r.a2 - r.a1);
                            result /= std::max (1.0 + increment, 0.05);
                        }
                    }
                    expectedY1[static_cast<std::size_t> (index)] = result;
                    expectedY2[static_cast<std::size_t> (index)] = index < active ? r.y1 : r.y2;
                }
                TaikoEngine::advanceModalBank (voice, ripple, rear, scale, parametric, rim);
                for (int index = 0; index < voice.modeCount; ++index)
                {
                    const auto& r = voice.modes[static_cast<std::size_t> (index)].resonator;
                    if (std::bit_cast<std::uint64_t> (r.y1) != std::bit_cast<std::uint64_t> (expectedY1[static_cast<std::size_t> (index)])
                        || std::bit_cast<std::uint64_t> (r.y2) != std::bit_cast<std::uint64_t> (expectedY2[static_cast<std::size_t> (index)]))
                        return false;
                }
                for (std::size_t i = 0; i < expectedInput.size(); ++i)
                    if (std::bit_cast<std::uint32_t> (expectedInput[i]) != std::bit_cast<std::uint32_t> (voice.modalInput[i]))
                        return false;
                for (std::size_t i = 0; i < expectedRim.size(); ++i)
                    if (std::bit_cast<std::uint32_t> (expectedRim[i]) != std::bit_cast<std::uint32_t> (voice.shellRimForce[i]))
                        return false;
                ++checked;
            }
        }
        std::cout << "Split modal advance is bit-identical across " << checked << " cases\n";
        return true;
    }
};
} // namespace taikor

int main()
{
    if (taikor::TaikoEngineTestAccess::check()) return 0;
    std::cerr << "FAIL: modal advance differs from the original recurrence\n";
    return 1;
}
