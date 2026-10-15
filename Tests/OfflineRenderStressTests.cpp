#include "DSP/EnsembleEngine.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cfenv>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <latch>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#if defined (__SSE__)
 #include <xmmintrin.h>
#endif

namespace taikor
{
// Read-only observations ensure the pending-event and parallel-eligible paths
// were reached. Do not manufacture private state or bypass host lifecycle rules.
struct EnsembleEngineTestAccess
{
    static std::size_t pendingHits (const EnsembleEngine& engine)
    { return engine.pendingCount; }
    static int size (const EnsembleEngine& engine)
    { return engine.parameters.ensembleSize; }
    static float variation (const EnsembleEngine& engine)
    { return engine.parameters.ensembleVariation; }
    static bool offline (const EnsembleEngine& engine)
    { return engine.offlineRendering; }
    static int activeCompanions (const EnsembleEngine& engine)
    {
        int count = 0;
        for (int member = 1; member < maximumEnsembleSize; ++member)
            count += engine.players[member]->getActiveVoiceCount() > 0 ? 1 : 0;
        return count;
    }
};
}

namespace
{
using Engine = taikor::EnsembleEngine;
using Access = taikor::EnsembleEngineTestAccess;
constexpr int hostThreads = 2;
constexpr int defaultIterations = 2;
constexpr std::uint32_t defaultSeed = 0x7461696bu;
constexpr float guard = 1234.0f;

// Explicit PRNG, independent of scheduling and standard-library distributions.
struct Random
{
    std::uint32_t state;
    std::uint32_t next()
    {
        state += 0x9e3779b9u;
        auto value = state;
        value = (value ^ (value >> 16)) * 0x21f0aaadu;
        value = (value ^ (value >> 15)) * 0x735a2d97u;
        return value ^ (value >> 15);
    }
};

std::uint64_t denormalControl()
{
#if defined (__SSE__)
    return _mm_getcsr() & 0x8040u; // FTZ and DAZ.
#elif defined (__aarch64__) && (defined (__clang__) || defined (__GNUC__))
    std::uint64_t value = 0;
    asm volatile ("mrs %0, fpcr" : "=r" (value));
    return value & (std::uint64_t { 1 } << 24); // FZ.
#else
    return 0;
#endif
}

void setDenormalControl (bool flush)
{
#if defined (__SSE__)
    _mm_setcsr ((_mm_getcsr() & ~0x8040u) | (flush ? 0x8040u : 0u));
#elif defined (__aarch64__) && (defined (__clang__) || defined (__GNUC__))
    std::uint64_t value = 0;
    asm volatile ("mrs %0, fpcr" : "=r" (value));
    constexpr auto mask = std::uint64_t { 1 } << 24;
    value = (value & ~mask) | (flush ? mask : 0);
    asm volatile ("msr fpcr, %0" : : "r" (value));
#else
    (void) flush;
#endif
}

struct FloatingEnvironment
{
    std::fenv_t saved {};
    bool valid = std::fegetenv (&saved) == 0;
    ~FloatingEnvironment()
    {
        if (valid)
            (void) std::fesetenv (&saved);
    }
};

struct Result
{
    std::uint64_t failures = 0, samples = 0, eligibleBlocks = 0, fallbackPools = 0;
    std::uint64_t digest = 14695981039346656037ull;
    int lane = 0, iteration = 0, phase = -1, block = 0;
    std::vector<std::string> messages;

    void expect (bool condition, const std::string& message)
    {
        if (condition)
            return;
        ++failures;
        if (messages.size() < 12)
        {
            std::ostringstream text;
            text << "lane=" << lane << " iteration=" << iteration
                 << " phase=" << phase << " block=" << block << ": " << message;
            messages.push_back (text.str());
        }
    }
};

void compare (Engine& serial, Engine& candidate, int count, Result& result,
              bool silent = false)
{
    std::array<float, 1026> left, right, otherLeft, otherRight;
    for (auto* buffer : { &left, &right, &otherLeft, &otherRight })
        buffer->fill (guard);

    // Workers were prepared under a different environment. Change controls
    // between successive jobs, including jobs with independent host callers.
    constexpr std::array modes { FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO };
    const auto selector = result.block + result.lane;
    const int rounding = modes[static_cast<std::size_t> (selector) % modes.size()];
    FloatingEnvironment restore;
    result.expect (restore.valid, "capture host floating-point environment");
    result.expect (std::fesetround (rounding) == 0, "set host rounding mode");
    setDenormalControl ((selector / 4) % 2 != 0);
    const auto denormals = denormalControl();
    std::fenv_t renderEnvironment {};
    const bool captured = std::fegetenv (&renderEnvironment) == 0;
    result.expect (captured, "capture render environment");

    if (count >= 64 && candidate.getOfflineWorkerCount() > 0
        && Access::offline (candidate) && Access::pendingHits (candidate) == 0
        && Access::activeCompanions (candidate) >= 2)
        ++result.eligibleBlocks;

    serial.process (left.data() + 1, right.data() + 1, count);
    result.expect (std::fegetround() == rounding && denormalControl() == denormals,
                   "serial render preserves caller floating-point controls");
    if (captured)
        result.expect (std::fesetenv (&renderEnvironment) == 0, "restore comparison environment");
    std::this_thread::yield();
    candidate.process (otherLeft.data() + 1, otherRight.data() + 1, count);
    result.expect (std::fegetround() == rounding && denormalControl() == denormals,
                   "offline render preserves caller floating-point controls");

    int firstMismatch = -1;
    bool finite = true, zero = true;
    for (int sample = 1; sample <= count; ++sample)
    {
        const auto a = std::bit_cast<std::uint32_t> (left[sample]);
        const auto b = std::bit_cast<std::uint32_t> (right[sample]);
        if ((a != std::bit_cast<std::uint32_t> (otherLeft[sample])
             || b != std::bit_cast<std::uint32_t> (otherRight[sample])) && firstMismatch < 0)
            firstMismatch = sample - 1;
        finite = finite && std::isfinite (left[sample]) && std::isfinite (right[sample])
                        && std::isfinite (otherLeft[sample]) && std::isfinite (otherRight[sample]);
        zero = zero && left[sample] == 0.0f && right[sample] == 0.0f
                    && otherLeft[sample] == 0.0f && otherRight[sample] == 0.0f;
        result.digest = (result.digest ^ a) * 1099511628211ull;
        result.digest = (result.digest ^ b) * 1099511628211ull;
        ++result.samples;
    }
    result.expect (firstMismatch < 0, "serial/offline bit mismatch at sample " + std::to_string (firstMismatch));
    result.expect (finite, "all four channels remain finite");
    result.expect (! silent || zero, "panic/reset/unprepared output is exactly silent");
    for (const auto* buffer : { &left, &right, &otherLeft, &otherRight })
        result.expect (buffer->front() == guard
                           && std::all_of (buffer->begin() + count + 1, buffer->end(),
                                          [] (float sample) { return sample == guard; }),
                       "render preserves prefix and entire unused buffer suffix");
    result.expect (serial.getActiveVoiceCount() == candidate.getActiveVoiceCount(),
                   "published active voice counts match");
    result.expect (Access::pendingHits (serial) == Access::pendingHits (candidate),
                   "pending event counts match");
    ++result.block;
}

void preparePool (Engine& engine, int requested, Result& result)
{
    engine.prepareOfflineRendering (requested);
    const auto cores = std::thread::hardware_concurrency();
    const int cap = cores == 0 ? 3 : static_cast<int> (std::min (cores - 1u, 3u));
    const int expected = std::clamp (requested, 0, cap);
    const int actual = engine.getOfflineWorkerCount();
    result.expect (actual == 0 || actual == expected,
                   "worker count equals hardware/request clamp or complete serial fallback");
    result.expect (! Access::offline (engine), "pool replacement disables offline selection");
    if (expected > 0 && actual == 0)
        ++result.fallbackPools;
}

void lifecycle (Random& random, Result& result)
{
    auto serial = std::make_unique<Engine>();
    auto candidate = std::make_unique<Engine>();
    // Pool creation before DSP prepare, including silent unprepared processing.
    preparePool (*candidate, 3, result);
    candidate->setOfflineRendering (true);
    compare (*serial, *candidate, 65, result, true);
    taikor::EngineParameters parameters;
    parameters.ensembleSize = taikor::maximumEnsembleSize;
    parameters.ensembleVariation = 0.0f;
    for (auto* engine : { serial.get(), candidate.get() })
    {
        engine->setParameters (parameters);
        engine->prepare (48000.0, 64);
    }

    constexpr std::array workers { std::numeric_limits<int>::min(), 0, 1, 2, 3,
                                  std::numeric_limits<int>::max() };
    constexpr std::array sizes { std::numeric_limits<int>::min(), 0, 1, 2,
                                taikor::maximumEnsembleSize, std::numeric_limits<int>::max() };
    const std::array variations { -1.0f, 0.0f, 1.0f, 2.0f,
                                 std::numeric_limits<float>::quiet_NaN(),
                                 std::numeric_limits<float>::infinity() };
    constexpr std::array blocks { 0, 1, 63, 64, 65, 255, 256, 257 };
    for (std::size_t phase = 0; phase < workers.size(); ++phase)
    {
        result.phase = static_cast<int> (phase);
        preparePool (*candidate, workers[phase], result);
        candidate->setOfflineRendering (true);
        parameters.ensembleSize = sizes[phase];
        parameters.ensembleVariation = variations[phase];
        for (auto* engine : { serial.get(), candidate.get() })
        {
            engine->setParameters (parameters);
            result.expect (Access::size (*engine) == std::clamp (sizes[phase], 1, taikor::maximumEnsembleSize),
                           "ensemble size clamps at integer boundaries");
            result.expect (std::isfinite (Access::variation (*engine))
                               && Access::variation (*engine) >= 0.0f && Access::variation (*engine) <= 1.0f,
                           "variation sanitises nonfinite/out-of-range boundaries");
        }
        compare (*serial, *candidate, 64, result);
        parameters.ensembleSize = taikor::maximumEnsembleSize;
        parameters.ensembleVariation = 0.0f;
        const auto articulation = static_cast<taikor::Articulation> (random.next() % 4u);
        const int octave = static_cast<int> (random.next() % 4u);
        const float velocity = 0.4f + static_cast<float> (random.next() % 512u) / 1024.0f;
        for (auto* engine : { serial.get(), candidate.get() })
        {
            engine->allSoundsOff();
            engine->setParameters (parameters);
            engine->trigger (articulation, octave, velocity);
        }
        result.expect (Access::activeCompanions (*candidate) >= 2,
                       "immediate ensemble hit activates parallel-eligible companions");
        for (std::size_t block = 0; block < blocks.size(); ++block)
        {
            candidate->setOfflineRendering (block % 2 != 0);
            compare (*serial, *candidate, blocks[block], result);
        }

        // Queue delayed hits, shrink while pending, then grow before they fire.
        parameters.ensembleVariation = 1.0f;
        for (auto* engine : { serial.get(), candidate.get() })
        {
            engine->setParameters (parameters);
            engine->trigger (taikor::Articulation::Don, octave, velocity);
            result.expect (Access::pendingHits (*engine) > 0, "delayed hits were actually queued");
        }
        parameters.ensembleSize = 1;
        for (auto* engine : { serial.get(), candidate.get() })
            engine->setParameters (parameters);
        compare (*serial, *candidate, 65, result);
        parameters.ensembleSize = taikor::maximumEnsembleSize;
        for (auto* engine : { serial.get(), candidate.get() })
            engine->setParameters (parameters);
        for (const int count : { 257, 1024, 257 })
            compare (*serial, *candidate, count, result);
        result.expect (Access::pendingHits (*candidate) == 0, "queued hits drain past maximum delay");

        for (auto* engine : { serial.get(), candidate.get() })
        {
            engine->trigger (taikor::Articulation::Don, octave, velocity);
            result.expect (Access::pendingHits (*engine) > 0, "panic has pending work to cancel");
            engine->allSoundsOff();
            result.expect (Access::pendingHits (*engine) == 0 && engine->getActiveVoiceCount() == 0,
                           "panic immediately clears pending hits and active voices");
        }
        // Longer than 30 ms at 48 kHz: canceled companions must never reappear.
        compare (*serial, *candidate, 1024, result, true);
        compare (*serial, *candidate, 513, result, true);

        if (phase == 2 || phase == 4)
        {
            const int previousWorkers = candidate->getOfflineWorkerCount();
            for (auto* engine : { serial.get(), candidate.get() })
            {
                engine->trigger (taikor::Articulation::Don, octave, velocity);
                engine->prepare (48000.0, phase == 2 ? 1 : 1024);
                result.expect (Access::pendingHits (*engine) == 0, "reprepare cancels queued events");
            }
            result.expect (candidate->getOfflineWorkerCount() == previousWorkers,
                           "DSP reprepare retains idle pool");
            compare (*serial, *candidate, 257, result, true);
        }
        if (phase % 2 == 0)
        {
            candidate->releaseOfflineRendering();
            candidate->releaseOfflineRendering();
            result.expect (candidate->getOfflineWorkerCount() == 0 && ! Access::offline (*candidate),
                           "repeated release joins helpers and clears offline selection");
            candidate->setOfflineRendering (true);
            compare (*serial, *candidate, 65, result, true);
        }
    }
    // Destruction with an unreleased pool, ringing lead and queued companions.
    // All process calls are complete; destruction never races its own render.
    candidate->trigger (taikor::Articulation::Don, 0, 0.7f);
    result.expect (Access::pendingHits (*candidate) > 0, "destruction covers queued events");
}

bool parseUnsigned (std::string_view text, std::uint32_t& value)
{
    const auto result = std::from_chars (text.data(), text.data() + text.size(), value);
    return result.ec == std::errc {} && result.ptr == text.data() + text.size();
}
}

int main (int argc, char** argv)
{
    std::uint32_t iterations = defaultIterations, seed = defaultSeed;
    for (int argument = 1; argument < argc; ++argument)
    {
        const std::string_view option = argv[argument];
        if (option == "--help")
        {
            std::cout << "Usage: OfflineRenderStressTests [--iterations 1..100] [--seed uint32]\n"
                         "Defaults: 2 iterations per host thread, decimal seed 1952541035; 2 host threads.\n";
            return 0;
        }
        std::uint32_t value = 0;
        if ((option != "--iterations" && option != "--seed") || argument + 1 >= argc
            || ! parseUnsigned (argv[++argument], value)
            || (option == "--iterations" && (value == 0 || value > 100)))
        {
            std::cerr << "Invalid arguments; use --help.\n";
            return 2;
        }
        (option == "--iterations" ? iterations : seed) = value;
    }

    // Realism feature flags are process-global. Leave them fixed at their
    // defaults throughout construction, concurrent rendering and destruction.
    std::array<Result, hostThreads> results;
    std::array<std::jthread, hostThreads> hosts;
    std::latch start { 1 };
    try
    {
        for (int lane = 0; lane < hostThreads; ++lane)
        {
            results[lane].lane = lane;
            hosts[lane] = std::jthread ([&, lane]
            {
                start.wait();
                auto& result = results[lane];
                try
                {
                    Random random { seed ^ (0x9e3779b9u * static_cast<std::uint32_t> (lane + 1)) };
                    for (std::uint32_t iteration = 0; iteration < iterations; ++iteration)
                    {
                        result.iteration = static_cast<int> (iteration);
                        result.phase = -1;
                        lifecycle (random, result);
                        std::this_thread::yield();
                    }
                }
                catch (const std::exception& error)
                { result.expect (false, std::string ("exception: ") + error.what()); }
                catch (...)
                { result.expect (false, "unknown exception"); }
            });
        }
    }
    catch (const std::exception& error)
    {
        start.count_down();
        for (auto& host : hosts)
            if (host.joinable())
                host.join();
        std::cerr << "Cannot create stress host threads: " << error.what() << '\n';
        return 1;
    }
    start.count_down();
    for (auto& host : hosts)
        host.join();

    std::uint64_t failures = 0;
    std::cout << "Offline render stress: seed=" << seed << " iterations=" << iterations
              << " hostThreads=" << hostThreads << '\n';
    for (const auto& result : results)
    {
        failures += result.failures;
        for (const auto& message : result.messages)
            std::cerr << "FAIL: " << message << '\n';
        std::cout << "lane=" << result.lane << " stereoSamples=" << result.samples
                  << " parallelEligibleBlocks=" << result.eligibleBlocks
                  << " fallbackPools=" << result.fallbackPools
                  << " serialDigest=" << result.digest << " failures=" << result.failures << '\n';
    }
    return failures == 0 ? 0 : 1;
}
