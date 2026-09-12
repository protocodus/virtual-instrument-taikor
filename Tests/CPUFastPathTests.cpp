#include "PluginProcessor.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace
{
void require (bool condition, const char* message)
{
    if (! condition) throw std::runtime_error (message);
}

bool same (float a, float b)
{
    return std::bit_cast<std::uint32_t> (a) == std::bit_cast<std::uint32_t> (b);
}

void parameter (TaikorAudioProcessor& processor, const char* id, float value)
{
    auto* target = processor.parameters.getParameter (id);
    require (target != nullptr, "Missing parameter");
    target->setValueNotifyingHost (target->convertTo0to1 (value));
}

void checkMeterRecurrence()
{
    const juce::ScopedNoDenormals noDenormals;
    for (double rate : { 48000.0, 96000.0 })
        for (int block : { 64, 257 })
        {
            auto processor = std::make_unique<TaikorAudioProcessor>();
            processor->prepareToPlay (rate, block);
            const float release = static_cast<float> (std::exp (-1.0 / (0.22 * rate)));
            std::array<float, 2> expected {};
            juce::AudioBuffer<float> audio (2, block);
            juce::MidiBuffer midi;
            for (int index = 0; index < 48; ++index)
            {
                midi.clear();
                if (index == 2 || index == 38)
                    midi.addEvent (juce::MidiMessage::noteOn (1, 84, 0.8f), 17);
                if (index == 10)
                {
                    parameter (*processor, taikor::parameters::reverbRoom, 3.0f);
                    parameter (*processor, taikor::parameters::reverbMix, 0.6f);
                }
                if (index == 28)
                {
                    midi.addEvent (juce::MidiMessage::allSoundOff (1), 31);
                    parameter (*processor, taikor::parameters::reverbRoom, 0.0f);
                }
                processor->processBlock (audio, midi);
                for (int sample = 0; sample < block; ++sample)
                    for (int channel = 0; channel < 2; ++channel)
                        expected[static_cast<std::size_t> (channel)] = std::max (
                            std::abs (audio.getSample (channel, sample)),
                            expected[static_cast<std::size_t> (channel)] * release);
                for (int channel = 0; channel < 2; ++channel)
                    require (same (expected[static_cast<std::size_t> (channel)],
                                   processor->getOutputLevel (channel)),
                             "Meter must retain its exact sample recurrence through MIDI, reverb and panic");
            }
        }

    // The block starts with a hit and ends after the physical bank freezes.
    // Checking frozen state only after rendering incorrectly drops its peaks.
    auto processor = std::make_unique<TaikorAudioProcessor>();
    constexpr int rate = 48000, samples = 13 * rate;
    processor->prepareToPlay (rate, 256);
    juce::AudioBuffer<float> audio (2, samples);
    juce::MidiBuffer midi;
    midi.addEvent (juce::MidiMessage::noteOn (1, 84, 0.8f), 0);
    processor->processBlock (audio, midi);
    const float release = static_cast<float> (std::exp (-1.0 / (0.22 * rate)));
    for (int channel = 0; channel < 2; ++channel)
    {
        float expected = 0.0f;
        for (int sample = 0; sample < samples; ++sample)
            expected = std::max (std::abs (audio.getSample (channel, sample)), expected * release);
        require (expected > 0.0f && same (expected, processor->getOutputLevel (channel)),
                 "A block ending in frozen output must retain earlier peaks");
    }
}

void checkReverbLifecycle()
{
    const juce::ScopedNoDenormals noDenormals;
    auto reused = std::make_unique<taikor::IRReverb>();
    for (double rate : { 48000.0, 96000.0, 48000.0 })
        for (int room : { 0, 1, 2, 3 })
        {
            auto fresh = std::make_unique<taikor::IRReverb>();
            const float mix = room == 0 ? 0.0f : 0.6f;
            reused->prepare (rate, 257, room, mix);
            fresh->prepare (rate, 257, room, mix);
            require (reused->isBypassed() == (room == 0), "Prepared bypass cache must match the room");
            for (int block = 0; block < 12; ++block)
            {
                if (block == 4)
                {
                    reused->setParameters (0, 0.0f);
                    fresh->setParameters (0, 0.0f);
                }
                // Calling the same setter every block must not restart fades.
                reused->setParameters (block < 4 ? room : 0, block < 4 ? mix : 0.0f);
                std::array<float, 257> a {}, b {}, c {}, d {};
                if (block == 0)
                    a[3] = b[3] = c[3] = d[3] = 0.25f;
                reused->process (a.data(), b.data(), static_cast<int> (a.size()));
                fresh->process (c.data(), d.data(), static_cast<int> (c.size()));
                for (std::size_t i = 0; i < a.size(); ++i)
                    require (same (a[i], c[i]) && same (b[i], d[i]),
                             "Cached reverb targets must survive prepare and repeated setters without changing audio");
            }
            reused->reset();
            require (reused->isBypassed(), "Reset must finish the fade to dry and clear its cached activity");
        }
}
} // namespace

int main()
{
    const juce::ScopedJuceInitialiser_GUI gui;
    try
    {
        checkMeterRecurrence();
        checkReverbLifecycle();
        std::cout << "CPU fast-path audio and meter checks passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
