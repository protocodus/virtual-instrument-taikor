#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>

namespace taikor
{
// Elapsed callback time relative to the audio deadline, not whole-machine CPU.
// Filter state belongs to the audio thread; the UI reads lock-free atomics.
class AudioCpuMeter
{
public:
    struct Snapshot
    {
        float average = 0.0f, peak = 0.0f;
        std::uint32_t overruns = 0;
    };
    void reset (double sampleRate = 48000.0) noexcept
    {
        const auto rate = std::isfinite (sampleRate)
            ? std::clamp (sampleRate, 8000.0, 768000.0) : 48000.0;
        inverseRate = 1.0 / rate;
        average = peak = 0.0;
        initialized = false;
        publishedAverage.store (0.0f, std::memory_order_relaxed);
        publishedPeak.store (0.0f, std::memory_order_relaxed);
        overruns.store (0, std::memory_order_relaxed);
    }
    void record (double elapsedSeconds, int samples) noexcept
    {
        if (samples <= 0 || ! std::isfinite (elapsedSeconds) || elapsedSeconds < 0.0)
            return;
        const double duration = static_cast<double> (samples) * inverseRate;
        const double load = 100.0 * elapsedSeconds / duration;
        if (! std::isfinite (load) || load > 1.0e9)
            return;
        average = initialized
            ? average + std::min (duration / 0.25, 1.0) * (load - average) : load;
        peak = std::max (load, peak * (1.0 - std::min (duration / 2.0, 1.0)));
        initialized = true;
        if (load > 100.0)
        {
            auto count = overruns.load (std::memory_order_relaxed);
            while (count != std::numeric_limits<std::uint32_t>::max()
                   && ! overruns.compare_exchange_weak (count, count + 1,
                                                        std::memory_order_relaxed,
                                                        std::memory_order_relaxed))
            {
            }
        }
        publishedAverage.store (static_cast<float> (average), std::memory_order_relaxed);
        publishedPeak.store (static_cast<float> (peak), std::memory_order_relaxed);
    }
    [[nodiscard]] Snapshot snapshot() const noexcept
    {
        return { publishedAverage.load (std::memory_order_relaxed),
                 publishedPeak.load (std::memory_order_relaxed),
                 overruns.load (std::memory_order_relaxed) };
    }
private:
    static_assert (std::atomic<float>::is_always_lock_free);
    static_assert (std::atomic<std::uint32_t>::is_always_lock_free);
    double inverseRate = 1.0 / 48000.0, average = 0.0, peak = 0.0;
    bool initialized = false;
    std::atomic<float> publishedAverage { 0.0f }, publishedPeak { 0.0f };
    std::atomic<std::uint32_t> overruns { 0 };
};
} // namespace taikor
