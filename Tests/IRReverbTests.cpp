#include "DSP/IRReverb.h"
#include "DSP/OutputLimiter.h"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
using taikor::IRReverb;
int failures = 0;

void expect (bool passed, const std::string& message)
{
    if (! passed)
    {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

struct Audio
{
    explicit Audio (int samples)
        : left (static_cast<std::size_t> (samples)), right (left.size()) {}

    std::vector<float> left, right;

    double energy (int begin = 0, int end = -1) const
    {
        if (end < 0)
            end = static_cast<int> (left.size());
        double result = 0.0;
        for (int index = begin; index < end; ++index)
        {
            const auto slot = static_cast<std::size_t> (index);
            result += static_cast<double> (left[slot]) * left[slot]
                    + static_cast<double> (right[slot]) * right[slot];
        }
        return result;
    }

    bool silent (int begin = 0) const
    {
        const auto zero = [] (float value) { return value == 0.0f; };
        return std::all_of (left.begin() + begin, left.end(), zero)
            && std::all_of (right.begin() + begin, right.end(), zero);
    }

    bool protectedOutput() const
    {
        const auto safe = [] (float value)
        {
            return std::isfinite (value)
                && std::abs (value) <= taikor::OutputLimiter::ceiling;
        };
        return std::all_of (left.begin(), left.end(), safe)
            && std::all_of (right.begin(), right.end(), safe);
    }

    bool sameBits (const Audio& other) const
    {
        const auto same = [] (float a, float b)
        { return std::bit_cast<std::uint32_t> (a) == std::bit_cast<std::uint32_t> (b); };
        return left.size() == other.left.size()
            && std::equal (left.begin(), left.end(), other.left.begin(), same)
            && std::equal (right.begin(), right.end(), other.right.begin(), same);
    }
};

Audio signal (int samples, double rate)
{
    Audio audio (samples);
    for (int index = 0; index < samples; ++index)
    {
        const double time = static_cast<double> (index) / rate;
        const auto slot = static_cast<std::size_t> (index);
        audio.left[slot] = static_cast<float> (
            0.11 * std::sin (6.283185307179586 * 293.0 * time)
            + 0.06 * std::cos (6.283185307179586 * 1733.0 * time));
        audio.right[slot] = static_cast<float> (
            0.09 * std::cos (6.283185307179586 * 397.0 * time)
            - 0.04 * std::sin (6.283185307179586 * 2137.0 * time));
        if (index % 197 == 0)
            audio.left[slot] += 0.19f;
        if (index % 311 == 0)
            audio.right[slot] -= 0.17f;
    }
    return audio;
}

void process (IRReverb& reverb, Audio& audio, int blockSize)
{
    const auto samples = static_cast<int> (audio.left.size());
    for (int offset = 0; offset < samples; offset += blockSize)
        reverb.process (audio.left.data() + offset, audio.right.data() + offset,
                        std::min (blockSize, samples - offset));
}

double maximumDifference (const Audio& a, const Audio& b)
{
    double result = 0.0;
    for (std::size_t sample = 0; sample < a.left.size(); ++sample)
        result = std::max ({ result,
            std::abs (static_cast<double> (a.left[sample]) - b.left[sample]),
            std::abs (static_cast<double> (a.right[sample]) - b.right[sample]) });
    return result;
}

void checkExactBypass()
{
    auto reverb = std::make_unique<IRReverb>();
    for (const double rate : { 48000.0, 96000.0 })
        for (const int room : { 0, 1, 2, 3 })
        {
            reverb->prepare (rate, 512, room, room == 0 ? 0.73f : 0.0f);
            expect (reverb->getLatency() == 0, "IR reverb must not add host latency");
            const auto input = [&]
            {
                auto audio = signal (8193, rate);
                // Bypass must also preserve signed zero and dry samples above
                // the wet-path limiter ceiling. A dry setting is an identity.
                audio.left[0] = -0.0f;
                audio.right[0] = 0.0f;
                audio.left[1] = 1.25f;
                audio.right[1] = -1.5f;
                return audio;
            }();
            for (const int blockSize : { 1, 64, 257, 512, 4097 })
            {
                reverb->reset();
                auto actual = input;
                process (*reverb, actual, blockSize);
                expect (actual.sameBits (input),
                        "Off and 0% must preserve every dry sample bit at rate "
                            + std::to_string (rate) + ", room " + std::to_string (room)
                            + ", block " + std::to_string (blockSize));
            }
        }
}

void checkRoomsAndEndpoints()
{
    for (const double rate : { 48000.0, 96000.0 })
    {
        auto reverb = std::make_unique<IRReverb>();
        std::vector<Audio> rooms;
        const int length = static_cast<int> (std::ceil (IRReverb::maximumTailSeconds * rate))
                         + 2048;
        expect (IRReverb::maximumTailSeconds > 0.1,
                "IR reverb must declare a nonzero room-tail bound");
        for (int room = 1; room <= 3; ++room)
        {
            reverb->prepare (rate, 512, room, 1.0f);
            Audio response (length);
            response.left[0] = response.right[0] = 0.25f;
            process (*reverb, response, 257);
            const auto label = "room " + std::to_string (room)
                             + " at " + std::to_string (rate) + " Hz";
            expect (response.protectedOutput() && response.energy() > 1.0e-10,
                    label + " must produce a finite, audible, protected impulse response");
            expect (response.energy (static_cast<int> (0.05 * rate)) > 1.0e-12,
                    label + " must retain a room tail after the input becomes silent");
            expect (response.energy (length - 1024) < 1.0e-16,
                    label + " must finish within its declared tail bound");
            double stereoDifference = 0.0;
            for (std::size_t index = 0; index < response.left.size(); ++index)
            {
                const double difference = static_cast<double> (response.left[index])
                                        - response.right[index];
                stereoDifference += difference * difference;
            }
            expect (stereoDifference > response.energy() * 1.0e-6,
                    label + " must retain the captured stereo room response");
            for (const auto& previous : rooms)
                expect (maximumDifference (response, previous) > 1.0e-6,
                        "the three rooms must have distinct responses");
            rooms.push_back (std::move (response));

            const auto input = signal (16387, rate);
            reverb->prepare (rate, 512, room, 1.0f);
            auto wet = input;
            process (*reverb, wet, 64);
            reverb->prepare (rate, 512, room, 0.5f);
            auto halfway = input;
            process (*reverb, halfway, 64);
            auto expected = input;
            for (std::size_t index = 0; index < input.left.size(); ++index)
            {
                expected.left[index] = 0.5f * (input.left[index] + wet.left[index]);
                expected.right[index] = 0.5f * (input.right[index] + wet.right[index]);
            }
            expect (maximumDifference (wet, input) > 1.0e-5,
                    label + " at 100% must contain the room response");
            expect (maximumDifference (halfway, expected) <= 2.0e-6,
                    label + " at 50% must be the linear dry/wet blend");

            reverb->reset();
            Audio afterReset (8193);
            process (*reverb, afterReset, 512);
            expect (afterReset.silent(), "reset must immediately clear " + label + " tails");

            reverb->prepare (rate, 512, room, 1.0f);
            auto overloaded = signal (8193, rate);
            for (auto& value : overloaded.left)
                value *= 10000.0f;
            for (auto& value : overloaded.right)
                value *= 10000.0f;
            process (*reverb, overloaded, 4097);
            expect (overloaded.protectedOutput(),
                    label + " wet overload must remain finite and limited");
        }
    }
}

void checkTailRetirement()
{
    auto reverb = std::make_unique<IRReverb>();
    for (const double rate : { 48000.0, 96000.0 })
        for (const bool roomOff : { false, true })
        {
            reverb->prepare (rate, 512, 2, 0.8f);
            auto initial = signal (4096, rate);
            process (*reverb, initial, 257);
            reverb->setParameters (roomOff ? 0 : 2, roomOff ? 0.8f : 0.0f);
            Audio release (static_cast<int> (0.5 * rate));
            process (*reverb, release, 257);
            expect (release.silent (static_cast<int> (0.4 * rate)),
                    "Off and 0% must settle to exact silence while the old room decays");

            auto dry = signal (4097, rate);
            const auto expected = dry;
            process (*reverb, dry, 4097);
            expect (dry.sameBits (expected),
                    "settled Off and 0% must preserve the dry signal exactly");
            reverb->setParameters (2, 0.8f);
            Audio returned (static_cast<int> (0.5 * rate));
            process (*reverb, returned, 64);
            expect (returned.silent(),
                    "reenabling a room must not replay tails or audio received during bypass");
        }
}

struct Event
{
    int sample, room;
    float mix;
};

Audio automated (IRReverb& reverb, double rate, int blockSize)
{
    constexpr std::array<Event, 11> events {{
        { 0, 1, 0.32f }, { 3, 2, 0.8f }, { 97, 3, 1.0f },
        { 991, 0, 0.7f }, { 2049, 2, 0.0f }, { 3077, 1, 0.65f },
        { 5003, 3, 0.2f }, { 7319, 2, 1.0f }, { 11987, 0, 0.8f },
        { 14701, 3, 0.5f }, { 19661, 1, 1.0f }
    }};
    reverb.setParameters (1, 0.32f);
    reverb.reset();
    auto audio = signal (static_cast<int> (0.5 * rate), rate);
    // Start and later return from silence while room weights are moving, so
    // the idle shortcut's smoothing is also checked against small host blocks.
    for (const auto range : { std::array<int, 2> { 0, 256 },
                              std::array<int, 2> { 12001, 15000 } })
    {
        std::fill (audio.left.begin() + range[0], audio.left.begin() + range[1], 0.0f);
        std::fill (audio.right.begin() + range[0], audio.right.begin() + range[1], 0.0f);
    }
    const int samples = static_cast<int> (audio.left.size());
    std::size_t nextEvent = 0;
    for (int offset = 0; offset < samples;)
    {
        while (nextEvent < events.size() && events[nextEvent].sample == offset)
        {
            reverb.setParameters (events[nextEvent].room, events[nextEvent].mix);
            ++nextEvent;
        }
        int count = std::min (blockSize, samples - offset);
        if (nextEvent < events.size())
            count = std::min (count, events[nextEvent].sample - offset);
        reverb.process (audio.left.data() + offset, audio.right.data() + offset, count);
        offset += count;
    }
    return audio;
}

void checkPartitionAndAutomation()
{
    double worstDifference = 0.0;
    for (const double rate : { 48000.0, 96000.0 })
    {
        auto reverb = std::make_unique<IRReverb>();
        reverb->prepare (rate, 512, 1, 0.32f);
        const auto reference = automated (*reverb, rate, 64);
        expect (reference.energy() > 1.0e-5 && reference.protectedOutput(),
                "rapid room and mix changes must remain audible, finite and limited");
        for (const int blockSize : { 1, 257, 512, 4097 })
        {
            const auto candidate = automated (*reverb, rate, blockSize);
            const double difference = maximumDifference (reference, candidate);
            worstDifference = std::max (worstDifference, difference);
            expect (candidate.protectedOutput() && difference <= 2.0e-6,
                    "room/mix automation must be independent of host block partition at "
                        + std::to_string (rate) + " Hz, block " + std::to_string (blockSize)
                        + " (maximum sample difference " + std::to_string (difference) + ")");
        }
    }
    std::cout << "IR block-partition maximum sample difference: "
              << std::scientific << worstDifference << '\n';
}

void reportCpu()
{
    // Informational only: wall-clock timings vary with CPU scheduling. Measure
    // the same five seconds of stereo input after prepare has loaded all IRs;
    // exclude allocation, input generation and result inspection from timing.
    constexpr double rate = 48000.0;
    constexpr int samples = static_cast<int> (5.0 * rate);
    const auto input = signal (samples, rate);
    auto reverb = std::make_unique<IRReverb>();
    for (int room = 0; room <= 3; ++room)
    {
        reverb->prepare (rate, 512, room, 0.2f);
        auto rendered = input;
        const auto start = std::chrono::steady_clock::now();
        process (*reverb, rendered, 256);
        const auto end = std::chrono::steady_clock::now();
        const double milliseconds = std::chrono::duration<double, std::milli> (end - start).count();
        expect (room == 0 ? rendered.sameBits (input) : rendered.protectedOutput(),
                "the CPU probe must produce valid dry or reverberant audio");
        std::cout << "IR room " << room << ": " << std::fixed << std::setprecision (3)
                  << milliseconds << " ms / 5 s audio (" << milliseconds / 50.0
                  << "% of one real-time core)\n";
        reverb->reset();
        Audio silence (samples);
        const auto idleStart = std::chrono::steady_clock::now();
        process (*reverb, silence, 256);
        const auto idleEnd = std::chrono::steady_clock::now();
        const double idleMilliseconds =
            std::chrono::duration<double, std::milli> (idleEnd - idleStart).count();
        expect (silence.silent(), "the idle CPU probe must remain exactly silent");
        std::cout << "IR room " << room << " idle: " << idleMilliseconds
                  << " ms / 5 s silence\n";
    }
}
} // namespace

int main()
{
    checkExactBypass();
    checkRoomsAndEndpoints();
    checkTailRetirement();
    checkPartitionAndAutomation();
    reportCpu();
    if (failures == 0)
        std::cout << "IR reverb checks passed\n";
    return failures == 0 ? 0 : 1;
}
