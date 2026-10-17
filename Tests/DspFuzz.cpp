// Stateful coverage-guided fuzzing of the public JUCE-free DSP boundary,
// modelled on Electry's DspFuzz. Every input starts from a fresh engine, so a
// saved libFuzzer crash replays without earlier corpus entries. Work is bounded
// independently of input length. Engine controls take raw IEEE floats and ints,
// out-of-range ones included, because setParameters promises to sanitise them;
// MIDI-derived values (velocity, bend, damping) and lifecycle calls stay inside
// what the plug-in processor can actually send.
// Run with -max_len=4096 -timeout=30 -rss_limit_mb=4096 and a corpus directory.
#include "DSP/EnsembleEngine.h"

#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <span>

namespace
{
struct Input
{
    std::span<const std::uint8_t> bytes;
    std::uint8_t byte()
    {
        if (bytes.empty())
            return 0;
        const auto result = bytes.front();
        bytes = bytes.subspan (1);
        return result;
    }
    std::uint32_t word()
    {
        std::uint32_t result = 0;
        for (int shift = 0; shift != 32; shift += 8)
            result |= static_cast<std::uint32_t> (byte()) << shift;
        return result;
    }
    float unit() { return static_cast<float> (byte()) / 255.0f; }
    float bipolar() { return unit() * 2.0f - 1.0f; }
    float scalar()
    {
        switch (byte() % 8)
        {
            case 0: return std::bit_cast<float> (word());
            case 1: return std::numeric_limits<float>::max();
            case 2: return -std::numeric_limits<float>::max();
            case 3: return std::numeric_limits<float>::quiet_NaN();
            case 4: return std::numeric_limits<float>::infinity();
            case 5: return -std::numeric_limits<float>::infinity();
            default: return unit();
        }
    }
};

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
constexpr std::array rates { 22050.0, 44100.0, 48000.0, 88200.0, 96000.0, 176400.0, 192000.0 };

void require (bool value)
{
    if (! value)
        std::abort();
}
} // namespace

extern "C" int LLVMFuzzerTestOneInput (const std::uint8_t* data, std::size_t size)
{
    if (size == 0)
        return 0;
    Input input { { data, size } };
    auto engine = std::make_unique<taikor::EnsembleEngine>();
    EngineParameters parameters;
    int maxBlock = 1 + input.byte();
    engine->prepare (rates[input.byte() % rates.size()], maxBlock);
    int preparesLeft = 2; // Preparing rebuilds every voice bank; keep inputs fast.

    // Sentinels catch writes past the caller's frame count that stay inside the
    // allocation and would otherwise evade AddressSanitizer.
    std::array<float, 258> left {}, right {};
    constexpr float sentinel = 12345.0f;

    for (int command = 0; command < 64 && ! input.bytes.empty(); ++command)
    {
        switch (input.byte() % 14)
        {
            case 0:
                parameters.*floatControls[input.byte() % floatControls.size()] = input.scalar();
                engine->setParameters (parameters);
                break;
            case 1:
                parameters.*intControls[input.byte() % intControls.size()] =
                    static_cast<int> (input.word());
                engine->setParameters (parameters);
                break;
            case 2:
                (void) engine->triggerMidi (input.byte() % 128, static_cast<float> (1 + input.byte() % 127) / 127.0f);
                break;
            case 3:
                engine->trigger (static_cast<taikor::Articulation> (input.byte() % taikor::articulationCount),
                                 static_cast<int> (input.byte() % 7) - 3, input.unit());
                break;
            case 4: engine->setPitchBend (input.bipolar()); break;
            case 5: engine->setHandDamping (input.unit()); break;
            case 6: engine->setStrikePositionOverride (input.bipolar()); break;
            case 7: engine->setStrikeAzimuthOverride (input.bipolar() * 3.14159265f); break;
            case 8: engine->clearStrikeOverrides(); break;
            case 9: engine->setRearHeadStrike ((input.byte() & 1) != 0); break;
            case 10: engine->reset(); break;
            case 11: engine->allSoundsOff(); break;
            case 12:
                if (preparesLeft-- > 0)
                {
                    maxBlock = 1 + input.byte();
                    engine->prepare (rates[input.byte() % rates.size()], maxBlock);
                }
                break;
            default:
            {
                const int count = 1 + input.byte() % maxBlock;
                left.fill (sentinel);
                right.fill (sentinel);
                engine->process (left.data() + 1, right.data() + 1, count);
                require (left.front() == sentinel && right.front() == sentinel);
                for (std::size_t i = 1; i < left.size(); ++i)
                {
                    if (i <= static_cast<std::size_t> (count))
                        require (std::isfinite (left[i]) && std::isfinite (right[i]));
                    else
                        require (left[i] == sentinel && right[i] == sentinel);
                }
                const int voices = engine->getActiveVoiceCount();
                require (voices >= 0);
                break;
            }
        }
    }
    return 0;
}
