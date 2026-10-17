#include "IRReverb.h"
#include <ImpulseData.h>

#include <algorithm>
#include <cmath>

namespace taikor
{
IRReverb::IRReverb()
{
    for (auto& room : rooms)
        room = std::make_unique<juce::dsp::Convolution> (
            juce::dsp::Convolution::NonUniform { 256 }, loader);
    weights[0].setCurrentAndTargetValue (1.0f);
}

void IRReverb::prepare (double sampleRate, int maximumBlockSize, int room, float mix)
{
    prepared = false;
    const double rate = std::isfinite (sampleRate)
        ? std::clamp (sampleRate, 8000.0, 384000.0) : 48000.0;
    chunkCapacity = std::clamp (maximumBlockSize, 1, blockCapacity);
    const juce::dsp::ProcessSpec spec { rate, static_cast<juce::uint32> (chunkCapacity), 2 };
    const std::array<const char*, roomCount> data {
        TaikorImpulseData::hall_wav, TaikorImpulseData::theater_wav, TaikorImpulseData::opera_wav };
    const std::array<int, roomCount> sizes {
        TaikorImpulseData::hall_wavSize, TaikorImpulseData::theater_wavSize,
        TaikorImpulseData::opera_wavSize };
    tailLengthSamples = 0;
    for (int index = 0; index < roomCount; ++index)
    {
        auto& convolution = *rooms[static_cast<std::size_t> (index)];
        // Full original length and stereo channels, including natural arrival
        // delay. JUCE's energy normalization supplies a useful wet gain without
        // shortening or otherwise editing the supplied room recordings.
        convolution.loadImpulseResponse (data[index], static_cast<std::size_t> (sizes[index]),
            juce::dsp::Convolution::Stereo::yes, juce::dsp::Convolution::Trim::no, 0,
            juce::dsp::Convolution::Normalise::yes);
        // JUCE guarantees the preceding IR is ready synchronously on prepare.
        convolution.prepare (spec);
        jassert (convolution.getLatency() == 0);
        tailLengthSamples = std::max (tailLengthSamples, convolution.getCurrentIRSize() + 1);
    }
    limiter.prepare (rate);
    for (auto& weight : weights)
        weight.reset (rate, 0.020);
    setParameters (room, mix);
    prepared = true;
    reset();
}

void IRReverb::setParameters (int room, float mix) noexcept
{
    room = std::clamp (room, 0, roomCount);
    mix = std::isfinite (mix) ? std::clamp (mix, 0.0f, 1.0f) : 0.0f;
    if (room == selectedRoom && mix == selectedMix)
        return;
    selectedRoom = room;
    selectedMix = mix;
    const float amount = selectedRoom == 0 ? 0.0f : selectedMix;
    weights[0].setTargetValue (1.0f - amount);
    for (int index = 0; index < roomCount; ++index)
        weights[static_cast<std::size_t> (index + 1)].setTargetValue (
            selectedRoom == index + 1 ? amount : 0.0f);
    updateActivity();
}

void IRReverb::updateActivity() noexcept
{
    activeRoomCount = 0;
    smoothing = weights[0].isSmoothing();
    for (int index = 0; index < roomCount; ++index)
    {
        const auto& weight = weights[static_cast<std::size_t> (index + 1)];
        smoothing = smoothing || weight.isSmoothing();
        if (weight.getCurrentValue() != 0.0f || weight.getTargetValue() != 0.0f)
            activeRooms[static_cast<std::size_t> (activeRoomCount++)] = index;
    }
    bypassed = prepared && ! smoothing && activeRoomCount == 0 && tailRemaining == 0
        && weights[0].getCurrentValue() == 1.0f
        && std::none_of (dirty.begin(), dirty.end(), [] (bool value) { return value; });
}

void IRReverb::clearHistory() noexcept
{
    for (int index = 0; index < roomCount; ++index)
        if (dirty[static_cast<std::size_t> (index)])
        {
            rooms[static_cast<std::size_t> (index)]->reset();
            dirty[static_cast<std::size_t> (index)] = false;
        }
    tailRemaining = 0;
    limiter.reset();
}

void IRReverb::reset() noexcept
{
    clearHistory();
    for (auto& weight : weights)
        weight.setCurrentAndTargetValue (weight.getTargetValue());
    updateActivity();
}

void IRReverb::process (float* left, float* right, int numSamples) noexcept
{
    if (bypassed || ! prepared || left == nullptr || right == nullptr || numSamples <= 0)
        return;

    for (int offset = 0; offset < numSamples;)
    {
        const int count = std::min (chunkCapacity, numSamples - offset);
        if (weights[0].getCurrentValue() == 1.0f && ! weights[0].isSmoothing())
        {
            // Off and settled 0% are exact dry bypass, with no FFT work and no
            // wet history that could reappear when the effect is enabled later.
            clearHistory();
            updateActivity();
            return;
        }

        // Ordinary audio ends in a nonzero sample: locating the last input
        // backwards avoids walking a whole sounding block for its tail clock.
        int lastInput = count - 1;
        while (lastInput >= 0 && left[offset + lastInput] == 0.0f && right[offset + lastInput] == 0.0f)
            --lastInput;
        if (lastInput < 0 && tailRemaining == 0)
        {
            if (smoothing)
            {
                for (auto& weight : weights)
                    weight.skip (count);
                updateActivity();
            }
            offset += count;
            continue;
        }

        const bool ramping = smoothing;
        const float dryGain = weights[0].getCurrentValue();
        if (ramping)
            for (std::size_t index = 0; index < weights.size(); ++index)
                for (int sample = 0; sample < count; ++sample)
                    gains[index][static_cast<std::size_t> (sample)] = weights[index].getNextValue();

        const float* inputs[] { left + offset, right + offset };
        const juce::dsp::AudioBlock<const float> input { inputs, 2, static_cast<std::size_t> (count) };
        for (int roomIndex = 0; roomIndex < activeRoomCount; ++roomIndex)
        {
            const int index = activeRooms[static_cast<std::size_t> (roomIndex)];
            auto* wetLeft = wet[static_cast<std::size_t> (2 * index)].data();
            auto* wetRight = wet[static_cast<std::size_t> (2 * index + 1)].data();
            float* channels[] { wetLeft, wetRight };
            juce::dsp::AudioBlock<float> block { channels, 2, static_cast<std::size_t> (count) };
            rooms[static_cast<std::size_t> (index)]->process (
                juce::dsp::ProcessContextNonReplacing<float> { input, block });
            dirty[static_cast<std::size_t> (index)] = true;
        }
        for (int sample = 0; sample < count; ++sample)
        {
            const auto position = static_cast<std::size_t> (sample);
            // Once the crossfade reaches dry, even its remaining samples must
            // bypass the extra limiter to retain the original dry waveform.
            const float dry = ramping ? gains[0][position] : dryGain;
            if (dry == 1.0f)
                continue;
            float l = left[offset + sample] * dry;
            float r = right[offset + sample] * dry;
            for (int roomIndex = 0; roomIndex < activeRoomCount; ++roomIndex)
            {
                const int index = activeRooms[static_cast<std::size_t> (roomIndex)];
                const auto weightIndex = static_cast<std::size_t> (index + 1);
                const float wetGain = ramping ? gains[weightIndex][position] : weights[weightIndex].getCurrentValue();
                l += wet[static_cast<std::size_t> (2 * index)][position] * wetGain;
                r += wet[static_cast<std::size_t> (2 * index + 1)][position] * wetGain;
            }
            limiter.process (l, r);
            left[offset + sample] = l;
            right[offset + sample] = r;
        }
        for (int index = 0; index < roomCount; ++index)
            if (dirty[static_cast<std::size_t> (index)]
                && weights[static_cast<std::size_t> (index + 1)].getCurrentValue() == 0.0f
                && ! weights[static_cast<std::size_t> (index + 1)].isSmoothing())
            {
                rooms[static_cast<std::size_t> (index)]->reset();
                dirty[static_cast<std::size_t> (index)] = false;
            }
        tailRemaining = lastInput >= 0 ? std::max (0, tailLengthSamples - (count - 1 - lastInput))
                                      : std::max (0, tailRemaining - count);
        if (tailRemaining == 0)
            clearHistory();
        // With settled weights the active set cannot change. Tail expiry
        // clears convolution history but does not change any room's weight.
        if (ramping)
            updateActivity();
        offset += count;
    }
}

bool IRReverb::isBypassed() const noexcept
{
    return bypassed;
}
} // namespace taikor
