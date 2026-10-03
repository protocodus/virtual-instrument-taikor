#include "DSP/EnsembleEngine.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>

namespace
{
int failures = 0;
std::uint64_t compared = 0;

void expect (bool condition, const char* message)
{
    if (! condition)
    {
        ++failures;
        if (failures <= 20)
            std::cerr << "FAIL: " << message << '\n';
    }
}

void compare (taikor::EnsembleEngine& serial, taikor::EnsembleEngine& parallel, int count)
{
    std::array<float, 1026> left {}, right {}, parallelLeft {}, parallelRight {};
    left.front() = left.back() = right.front() = right.back() = 1234.0f;
    parallelLeft.front() = parallelLeft.back() = parallelRight.front() = parallelRight.back() = 1234.0f;
    serial.process (left.data() + 1, right.data() + 1, count);
    parallel.process (parallelLeft.data() + 1, parallelRight.data() + 1, count);
    for (int sample = 1; sample <= count; ++sample)
    {
        expect (std::isfinite (left[sample]) && std::isfinite (right[sample]), "finite ensemble output");
        expect (std::abs (left[sample]) <= taikor::OutputLimiter::ceiling
                    && std::abs (right[sample]) <= taikor::OutputLimiter::ceiling,
                "shared output limiter remains active");
        expect (std::bit_cast<std::uint32_t> (left[sample])
                    == std::bit_cast<std::uint32_t> (parallelLeft[sample])
                && std::bit_cast<std::uint32_t> (right[sample])
                    == std::bit_cast<std::uint32_t> (parallelRight[sample]),
                "parallel rendering must preserve exact sample order and values");
        ++compared;
    }
    expect (left.front() == 1234.0f && left.back() == 1234.0f
                && right.front() == 1234.0f && right.back() == 1234.0f
                && parallelLeft.front() == 1234.0f && parallelLeft.back() == 1234.0f
                && parallelRight.front() == 1234.0f && parallelRight.back() == 1234.0f,
            "render buffers keep their guard samples");
}

void exercise (double rate, int workers)
{
    auto serial = std::make_unique<taikor::EnsembleEngine>();
    auto parallel = std::make_unique<taikor::EnsembleEngine>();
    taikor::EngineParameters parameters;
    parameters.ensembleSize = taikor::maximumEnsembleSize;
    parameters.ensembleVariation = 0.7f;
    serial->setParameters (parameters);
    parallel->setParameters (parameters);
    serial->prepare (rate, 64);
    parallel->prepare (rate, 64);
    parallel->prepareOfflineRendering (workers);
    parallel->setOfflineRendering (true);
    expect (parallel->getOfflineWorkerCount() >= 0 && parallel->getOfflineWorkerCount() <= 3,
            "worker count is bounded, including resource-exhaustion fallback");
    constexpr std::array sizes { 0, 1, 7, 63, 64, 127, 256, 257, 1024 };
    for (int event = 0; event < 36; ++event)
    {
        if (event % 6 == 0)
        {
            parameters.ensembleSize = event % 12 == 0 ? 8 : 3;
            parameters.ensembleVariation = event % 12 == 0 ? 0.7f : 0.0f;
            serial->setParameters (parameters);
            parallel->setParameters (parameters);
        }
        for (auto* engine : { serial.get(), parallel.get() })
        {
            engine->setRearHeadStrike (event % 3 == 0);
            engine->setHandDamping (event % 5 == 0 ? 0.8f : 0.0f);
            engine->setPitchBend (event % 4 == 0 ? 0.2f : 0.0f);
            engine->setStrikePositionOverride (0.05f * static_cast<float> (event % 9 - 4));
            engine->setStrikeAzimuthOverride (0.1f * static_cast<float> (event % 7 - 3));
            engine->trigger (static_cast<taikor::Articulation> (event % 4), event % 4, 0.8f);
        }
        compare (*serial, *parallel, sizes[static_cast<std::size_t> (event) % sizes.size()]);
        if (event == 17)
        {
            serial->allSoundsOff();
            parallel->allSoundsOff();
        }
        if (event == 25)
            parallel->setOfflineRendering (false);
        if (event == 27)
            parallel->setOfflineRendering (true);
    }
    serial->reset();
    parallel->reset();
    compare (*serial, *parallel, 1024);
    parallel->releaseOfflineRendering();
    expect (parallel->getOfflineWorkerCount() == 0, "release joins all worker threads");
    compare (*serial, *parallel, 1024);
}

void hostileParameters()
{
    auto engine = std::make_unique<taikor::EnsembleEngine>();
    std::array<float, 257> left {}, right {};
    engine->process (left.data(), right.data(), 257);
    expect (std::all_of (left.begin(), left.end(), [] (float value) { return value == 0.0f; }),
            "unprepared rendering is silent");
    engine->process (nullptr, right.data(), 257);
    engine->process (left.data(), nullptr, 257);
    engine->process (left.data(), right.data(), -1);
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float inf = std::numeric_limits<float>::infinity();
    for (const double rate : { 0.0, -1.0, 48000.0, static_cast<double> (nan), static_cast<double> (inf) })
    {
        engine->prepare (rate, -1);
        for (const float value : { nan, inf, -inf, -1.0e30f, 1.0e30f, 0.0f })
        {
            taikor::EngineParameters parameters;
            parameters.ensembleSize = value == 0.0f ? 0 : 999;
            parameters.ensembleVariation = value;
            parameters.tension = parameters.headMaterial = parameters.shellMaterial = value;
            parameters.headDamping = parameters.strikePosition = parameters.strikeAzimuth = value;
            engine->setParameters (parameters);
            engine->setStrikePositionOverride (value);
            engine->setStrikeAzimuthOverride (value);
            engine->setHandDamping (value);
            engine->setPitchBend (value);
            expect (! engine->triggerMidi (std::numeric_limits<int>::min(), 0.8f)
                        && ! engine->triggerMidi (std::numeric_limits<int>::max(), 0.8f),
                    "hostile MIDI note integers are rejected");
            engine->trigger (static_cast<taikor::Articulation> (255), 0, 0.8f);
            engine->trigger (taikor::Articulation::Don, 0, nan);
            engine->trigger (taikor::Articulation::Don, 0, 0.8f);
            engine->process (left.data(), right.data(), 257);
            const auto finite = [] (float sample) { return std::isfinite (sample); };
            expect (std::all_of (left.begin(), left.end(), finite)
                        && std::all_of (right.begin(), right.end(), finite),
                    "hostile ensemble parameters cannot poison output");
            engine->allSoundsOff();
        }
    }
}
} // namespace

int main()
{
    hostileParameters();
    for (const int workers : { -1, 1, 3, 999 })
        exercise (48000.0, workers);
    exercise (96000.0, 3);
    exercise (384000.0, 3);
    std::cout << "Ensemble robustness: " << compared << " stereo sample pairs compared; "
              << failures << " failures\n";
    return failures == 0 ? 0 : 1;
}
