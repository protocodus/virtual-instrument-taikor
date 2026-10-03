#include "PluginProcessor.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace
{
constexpr double sampleRate = 48000.0;
constexpr int totalSamples = 768;
constexpr int referenceNote = 48;

int failureCount = 0;

void expect (bool condition, const std::string& message)
{
    if (! condition)
    {
        ++failureCount;
        std::cerr << "FAIL: " << message << '\n';
    }
}

struct MidiEvent
{
    int sample = 0;
    juce::MidiMessage message;
};

MidiEvent at (int sample, juce::MidiMessage message)
{
    return { sample, std::move (message) };
}

juce::MidiMessage note (float velocity = 0.9f)
{
    return juce::MidiMessage::noteOn (1, referenceNote, velocity);
}

juce::MidiMessage controller (int number, int value)
{
    return juce::MidiMessage::controllerEvent (1, number, value);
}

std::vector<int> uniformPartitions (int blockSize)
{
    std::vector<int> result;
    for (int remaining = totalSamples; remaining > 0; remaining -= blockSize)
        result.push_back (std::min (blockSize, remaining));
    return result;
}

std::vector<float> render (const std::vector<MidiEvent>& events,
                           const std::vector<int>& partitions)
{
    TaikorAudioProcessor processor;
    processor.prepareToPlay (sampleRate, totalSamples);

    std::vector<float> output;
    output.reserve (static_cast<std::size_t> (totalSamples) * 2u);

    int blockStart = 0;
    for (const auto blockSize : partitions)
    {
        expect (blockSize > 0 && blockStart + blockSize <= totalSamples,
                "partition exceeds the rendered sample range");
        if (blockSize <= 0 || blockStart + blockSize > totalSamples)
            break;

        juce::AudioBuffer<float> audio (2, blockSize);
        audio.clear();
        juce::MidiBuffer midi;
        const auto blockEnd = blockStart + blockSize;

        // Reinsert in authored order. JUCE's same-sample ordering is part of
        // the host contract being tested here.
        for (const auto& event : events)
            if (event.sample >= blockStart && event.sample < blockEnd)
                midi.addEvent (event.message, event.sample - blockStart);

        processor.processBlock (audio, midi);
        for (int sample = 0; sample < blockSize; ++sample)
        {
            output.push_back (audio.getSample (0, sample));
            output.push_back (audio.getSample (1, sample));
        }
        blockStart = blockEnd;
    }

    expect (blockStart == totalSamples, "partitions do not cover the render");
    return output;
}

float peak (const std::vector<float>& audio)
{
    float result = 0.0f;
    for (const auto sample : audio)
        result = std::max (result, std::abs (sample));
    return result;
}

float maximumDifference (const std::vector<float>& first,
                         const std::vector<float>& second)
{
    expect (first.size() == second.size(), "audio captures have different lengths");
    const auto count = std::min (first.size(), second.size());
    float result = 0.0f;
    for (std::size_t index = 0; index < count; ++index)
        result = std::max (result, std::abs (first[index] - second[index]));
    return result;
}

void expectExact (const std::vector<float>& first, const std::vector<float>& second,
                  const std::string& message)
{
    if (first.size() != second.size())
    {
        expect (false, message + " (length)");
        return;
    }

    for (std::size_t index = 0; index < first.size(); ++index)
        if (first[index] != second[index])
        {
            expect (false, message + " (sample " + std::to_string (index / 2u) + ")");
            return;
        }
}

void expectOrderDifference (const std::vector<MidiEvent>& first,
                            const std::vector<MidiEvent>& second,
                            const std::string& message)
{
    const auto firstAudio = render (first, uniformPartitions (totalSamples));
    const auto secondAudio = render (second, uniformPartitions (totalSamples));
    expect (maximumDifference (firstAudio, secondAudio) > 1.0e-7f, message);
}

void testHostStableOrder()
{
    juce::MidiBuffer midi;
    midi.addEvent (controller (16, 127), 11);
    midi.addEvent (controller (17, 0), 11);
    midi.addEvent (controller (18, 100), 11);
    midi.addEvent (juce::MidiMessage::pitchWheel (1, 16383), 11);
    midi.addEvent (note(), 11);

    const std::vector<int> expectedControllers { 16, 17, 18 };
    std::vector<int> controllers;
    std::vector<bool> notePositions;
    juce::MidiMessage message;
    int position = 0;
    juce::MidiBuffer::Iterator iterator (midi);
    while (iterator.getNextEvent (message, position))
    {
        expect (position == 11, "same-sample MIDI event moved to another timestamp");
        if (message.isController())
            controllers.push_back (message.getControllerNumber());
        notePositions.push_back (message.isNoteOn());
    }

    expect (controllers == expectedControllers,
            "host MIDI controller insertion order was not retained");
    expect (notePositions.size() == 5u && ! notePositions[0] && ! notePositions[1]
                && ! notePositions[2] && ! notePositions[3] && notePositions[4],
            "host MIDI note order was not retained at a shared timestamp");
}

void testCapturedControllerOrder()
{
    const auto cc16Before = std::vector<MidiEvent> {
        at (0, controller (16, 127)), at (0, note())
    };
    const auto cc16After = std::vector<MidiEvent> {
        at (0, note()), at (0, controller (16, 127))
    };
    expectOrderDifference (cc16Before, cc16After,
                           "CC16 before/after note did not change captured azimuth");

    const auto cc17Before = std::vector<MidiEvent> {
        at (0, controller (17, 0)), at (0, note())
    };
    const auto cc17After = std::vector<MidiEvent> {
        at (0, note()), at (0, controller (17, 0))
    };
    expectOrderDifference (cc17Before, cc17After,
                           "CC17 before/after note did not change captured position");

    const auto cc18Before = std::vector<MidiEvent> {
        at (0, controller (18, 100)), at (0, note())
    };
    const auto cc18After = std::vector<MidiEvent> {
        at (0, note()), at (0, controller (18, 100))
    };
    expectOrderDifference (cc18Before, cc18After,
                           "CC18 before/after note did not capture the struck head");
}

void testSmoothedControllerOrder()
{
    const auto partitions = uniformPartitions (totalSamples);
    const auto baseline = render ({ at (0, note()) }, partitions);
    constexpr int delayedSample = 256;

    for (const auto& change : { juce::MidiMessage::pitchWheel (1, 16383),
                               controller (1, 127) })
    {
        const std::string name = change.isPitchWheel() ? "pitch wheel" : "CC1";
        const auto before = render ({ at (0, change), at (0, note()) }, partitions);
        const auto after = render ({ at (0, note()), at (0, change) }, partitions);

        // These setters change only smoothing targets. Unlike CC16/17/18's
        // captured strike attributes, they commute with a note at the same
        // timestamp: neither the dispatch nor the note advances a sample.
        expectExact (before, after,
                     name + " changed audio when reordered around a same-sample note");
        expect (maximumDifference (before, baseline) > 1.0e-7f,
                name + " at note onset had no audible effect");

        // Keep the note fixed and move only the controller. The ringing head
        // must respond once samples advance, without changing earlier audio.
        const std::vector<MidiEvent> delayedStream {
            at (0, note()), at (delayedSample, change)
        };
        const auto delayed = render (delayedStream, partitions);
        expect (std::equal (delayed.begin(), delayed.begin() + 2 * delayedSample,
                            baseline.begin()),
                name + " affected audio before its scheduled sample");
        expect (maximumDifference (delayed, baseline) > 1.0e-7f,
                name + " after note onset did not affect the ringing drum");
        expect (maximumDifference (before, delayed) > 1.0e-7f,
                name + " ignored the positive sample offset from note onset");
        expectExact (delayed, render (delayedStream, uniformPartitions (64)),
                     name + " timing changed when its event fell on a block boundary");
    }
}

void testResetAndPanicOrder()
{
    const std::vector<MidiEvent> configured {
        at (0, controller (16, 127)),
        at (0, controller (17, 0)),
        at (0, controller (18, 100)),
        at (0, controller (1, 127)),
        at (0, juce::MidiMessage::pitchWheel (1, 16383)),
    };

    auto resetBefore = configured;
    resetBefore.push_back (at (0, controller (121, 0)));
    resetBefore.push_back (at (0, note()));

    auto resetAfter = configured;
    resetAfter.push_back (at (0, note()));
    resetAfter.push_back (at (0, controller (121, 0)));
    expectOrderDifference (resetBefore, resetAfter,
                           "CC121 before/after note did not restore controllers in order");

    const auto panicBefore = std::vector<MidiEvent> {
        at (0, juce::MidiMessage::allSoundOff (1)), at (0, note())
    };
    const auto panicAfter = std::vector<MidiEvent> {
        at (0, note()), at (0, juce::MidiMessage::allSoundOff (1))
    };
    const auto audible = render (panicBefore, uniformPartitions (totalSamples));
    const auto silenced = render (panicAfter, uniformPartitions (totalSamples));
    expect (peak (audible) > 1.0e-7f, "note after panic was unexpectedly silent");
    expect (peak (silenced) == 0.0f,
            "panic after note did not clear the percussive voice immediately");
}

void testPercussiveNoteOffPolicy()
{
    const auto noteOnly = std::vector<MidiEvent> { at (0, note()) };
    const auto noteOff = std::vector<MidiEvent> {
        at (0, note()), at (0, juce::MidiMessage::noteOff (1, referenceNote))
    };
    const auto zeroVelocity = std::vector<MidiEvent> {
        at (0, note()), at (0, juce::MidiMessage::noteOn (1, referenceNote, 0.0f))
    };

    const auto baseline = render (noteOnly, uniformPartitions (totalSamples));
    expectExact (baseline, render (noteOff, uniformPartitions (totalSamples)),
                 "note-off truncated a percussive tail");
    expectExact (baseline, render (zeroVelocity, uniformPartitions (totalSamples)),
                 "zero-velocity note-on truncated a percussive tail");
}

std::vector<MidiEvent> boundaryStream()
{
    return {
        at (0, controller (16, 127)), at (0, note (0.72f)),
        at (64, controller (17, 0)), at (64, note (0.81f)),
        at (128, controller (18, 100)), at (128, note (0.63f)),
        at (192, juce::MidiMessage::pitchWheel (1, 0)),
        at (256, controller (1, 127)), at (256, note (0.88f)),
        at (384, controller (121, 0)), at (384, note (0.54f)),
        at (512, note (0.94f)), at (512, juce::MidiMessage::allSoundOff (1)),
        at (640, controller (18, 0)), at (640, note (0.67f)),
    };
}

void testBlockPartitionInvariance()
{
    const auto stream = boundaryStream();
    const auto whole = render (stream, uniformPartitions (totalSamples));
    const auto fixed = render (stream, uniformPartitions (64));
    const std::vector<int> irregular { 1, 7, 32, 64, 3, 127, 5, 256, 128, 145 };

    expectExact (whole, fixed,
                 "same MIDI stream changed when split into equal host blocks");
    expectExact (whole, render (stream, irregular),
                 "same MIDI stream changed at irregular block boundaries");
}

class FixedRandom
{
public:
    explicit FixedRandom (std::uint32_t seed) : state (seed) {}

    std::uint32_t next() noexcept
    {
        state = state * 1664525u + 1013904223u;
        return state;
    }

private:
    std::uint32_t state;
};

std::vector<MidiEvent> seededStream (std::uint32_t seed)
{
    FixedRandom random (seed);
    std::vector<MidiEvent> result;
    int sample = 0;

    for (int index = 0; index < 32; ++index)
    {
        sample += static_cast<int> (random.next() % 24u);
        if (sample >= totalSamples)
            break;

        const auto value = static_cast<int> (random.next() & 0x7fu);
        switch (index % 8)
        {
            case 0:
                result.push_back (at (sample, note (0.35f + 0.6f * (value / 127.0f))));
                break;
            case 1:
                result.push_back (at (sample, controller (1, value)));
                break;
            case 2:
                result.push_back (at (sample, controller (16, value)));
                break;
            case 3:
                result.push_back (at (sample, controller (17, value)));
                break;
            case 4:
                result.push_back (at (sample, controller (18, value)));
                break;
            case 5:
                result.push_back (at (sample, juce::MidiMessage::pitchWheel (
                                            1, value * 129)));
                break;
            case 6:
                result.push_back (at (sample, controller (121, 0)));
                break;
            default:
                result.push_back (at (sample, juce::MidiMessage::allSoundOff (1)));
                break;
        }
    }

    // Guarantee that every replay contains a same-sample note/controller pair
    // even if a future RNG implementation changes the generated spacing.
    result.push_back (at (totalSamples / 2, controller (18, 64)));
    result.push_back (at (totalSamples / 2, note (0.8f)));
    return result;
}

void testSeededReplay()
{
    constexpr std::uint32_t seed = 0x5EED1234u;
    const auto stream = seededStream (seed);
    const auto first = render (stream, uniformPartitions (64));
    const auto second = render (stream, uniformPartitions (64));
    expectExact (first, second,
                 "seeded valid MIDI replay was not deterministic for seed 0x5EED1234");

    const std::vector<int> replayPartitions { 17, 31, 64, 2, 127, 53, 256, 218 };
    expectExact (first, render (stream, replayPartitions),
                 "seeded valid MIDI replay was not partition invariant");
}
} // namespace

int main()
{
    testHostStableOrder();
    testCapturedControllerOrder();
    testSmoothedControllerOrder();
    testResetAndPanicOrder();
    testPercussiveNoteOffPolicy();
    testBlockPartitionInvariance();
    testSeededReplay();

    if (failureCount != 0)
    {
        std::cerr << failureCount << " MIDI ordering test(s) failed\n";
        return 1;
    }

    std::cout << "Taikor MIDI ordering tests passed (seed 0x5EED1234)\n";
    return 0;
}
