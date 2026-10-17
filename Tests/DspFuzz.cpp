// Stateful coverage-guided fuzzing of the public JUCE-free DSP boundary, on
// the shared Tests/FuzzHarness.h. Every input starts from a fresh engine, so a
// saved libFuzzer crash replays without earlier corpus entries. Work is bounded
// independently of input length. Engine controls take raw IEEE floats and ints,
// out-of-range ones included, because setParameters promises to sanitise them;
// MIDI-derived values (velocity, bend, damping) stay inside what the plug-in
// processor can actually send; prepare also takes the nonsense sample rates
// the harness offers, since a host can pass anything.
// Run with -max_len=4096 -timeout=30 -rss_limit_mb=4096 and a corpus directory.
#include "DSP/EnsembleEngine.h"
#include "FuzzHarness.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>

namespace
{
using taikor::EngineParameters;
constexpr std::array floatControls {
    &EngineParameters::headDiameter, &EngineParameters::bodyDepth,
    &EngineParameters::tension, &EngineParameters::headMaterial,
    &EngineParameters::shellMaterial, &EngineParameters::resonantTension,
    &EngineParameters::cavityCoupling, &EngineParameters::headDamping,
    &EngineParameters::shellResonance, &EngineParameters::pitch,
    &EngineParameters::bachiHardness, &EngineParameters::strikePosition,
    &EngineParameters::strikeAzimuth, &EngineParameters::velocityDepth,
    &EngineParameters::velocityCurve, &EngineParameters::tensionModulation,
    &EngineParameters::strikeNoise, &EngineParameters::humanise,
    &EngineParameters::octaveBody, &EngineParameters::micDistance,
    &EngineParameters::micSpread, &EngineParameters::stereoWidth,
    &EngineParameters::drive, &EngineParameters::outputGain,
    &EngineParameters::outputHighPassHz, &EngineParameters::ensembleVariation,
    &EngineParameters::physicalFamilyMix, &EngineParameters::physicalRearTensionScale
};
constexpr std::array intControls {
    &EngineParameters::performer, &EngineParameters::ensembleSize,
    &EngineParameters::physicalFamily
};
} // namespace

extern "C" int LLVMFuzzerTestOneInput (const std::uint8_t* data, std::size_t size)
{
    using namespace protocodus::fuzz;
    if (size == 0)
        return 0;
    Input input { data, size };
    auto engine = std::make_unique<taikor::EnsembleEngine>();
    EngineParameters parameters;
    int maxBlock = 1 + input.Byte();
    engine->prepare (input.Rate(), maxBlock);
    int preparesLeft = 2; // Preparing rebuilds every voice bank; keep inputs fast.
    Sentinels<256> block;

    for (int command = 0; command < 64 && ! input.Empty(); ++command)
    {
        switch (input.Pick (14))
        {
            case 0:
                parameters.*floatControls[input.Pick (floatControls.size())] = input.Scalar();
                engine->setParameters (parameters);
                break;
            case 1:
                parameters.*intControls[input.Pick (intControls.size())] = static_cast<int> (input.Word());
                engine->setParameters (parameters);
                break;
            case 2:
                (void) engine->triggerMidi (input.Byte() % 128, static_cast<float> (1 + input.Byte() % 127) / 127.0f);
                break;
            case 3:
                engine->trigger (static_cast<taikor::Articulation> (input.Pick (taikor::articulationCount)),
                                 static_cast<int> (input.Pick (7)) - 3, input.Byte() / 255.0f);
                break;
            case 4: engine->setPitchBend (input.Byte() / 127.5f - 1.0f); break;
            case 5: engine->setHandDamping (input.Byte() / 255.0f); break;
            case 6: engine->setStrikePositionOverride (input.Byte() / 127.5f - 1.0f); break;
            case 7: engine->setStrikeAzimuthOverride ((input.Byte() / 127.5f - 1.0f) * 3.14159265f); break;
            case 8: engine->clearStrikeOverrides(); break;
            case 9: engine->setRearHeadStrike ((input.Byte() & 1) != 0); break;
            case 10: engine->reset(); break;
            case 11: engine->allSoundsOff(); break;
            case 12:
                if (preparesLeft-- > 0)
                {
                    maxBlock = 1 + input.Byte();
                    engine->prepare (input.Rate(), maxBlock);
                }
                break;
            default:
            {
                const auto frames = 1 + input.Pick (static_cast<std::size_t> (maxBlock));
                block.Arm();
                engine->process (block.Left(), block.Right(), static_cast<int> (frames));
                block.Check (frames);
                Require (engine->getActiveVoiceCount() >= 0);
                break;
            }
        }
    }
    return 0;
}
