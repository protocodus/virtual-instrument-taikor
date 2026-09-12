#pragma once

#include <juce_dsp/juce_dsp.h>

#include "OutputLimiter.h"

#include <array>
#include <memory>

namespace taikor
{
// One shared room after the complete ensemble. The original stereo IRs are
// prepared up front; changing rooms never decodes files or builds FFT plans
// on the audio thread. Only audible rooms (including a short crossfade) run.
class IRReverb
{
public:
    IRReverb();
    void prepare (double sampleRate, int maximumBlockSize,
                  int room = 0, float mix = 0.2f);
    void setParameters (int room, float mix) noexcept;
    void reset() noexcept;
    void process (float* left, float* right, int numSamples) noexcept;
    [[nodiscard]] bool isBypassed() const noexcept;

    // Opera is the longest original: 88,594 frames at 44,100 Hz. Round up
    // enough to cover resampling to any supported host clock as well.
    static constexpr double maximumTailSeconds = 2.01;
    [[nodiscard]] int getLatency() const noexcept { return 0; }

private:
    static constexpr int roomCount = 3;
    static constexpr int blockCapacity = 256;
    void clearHistory() noexcept;
    void updateActivity() noexcept;

    // The shared loader must outlive all three convolution instances.
    juce::dsp::ConvolutionMessageQueue loader;
    std::array<std::unique_ptr<juce::dsp::Convolution>, roomCount> rooms;
    std::array<juce::SmoothedValue<float>, roomCount + 1> weights;
    std::array<std::array<float, blockCapacity>, roomCount + 1> gains {};
    std::array<std::array<float, blockCapacity>, roomCount * 2> wet {};
    std::array<bool, roomCount> dirty {};
    std::array<int, roomCount> activeRooms {};
    int activeRoomCount = 0;
    bool bypassed = false;
    bool smoothing = false;
    OutputLimiter limiter;
    int selectedRoom = 0;
    float selectedMix = 0.2f;
    int chunkCapacity = blockCapacity;
    int tailLengthSamples = 0;
    int tailRemaining = 0;
    bool prepared = false;
};
} // namespace taikor
