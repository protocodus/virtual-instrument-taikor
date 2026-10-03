#include "DSP/TaikoEngine.h"

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>

#if defined(_MSC_VER)
 #include <malloc.h>
#endif

// Allocation failure is injected only in this standalone test executable.
// Cover aligned allocation too: the optional modal storage can be over-aligned.
namespace allocation_failure
{
thread_local int remaining = -1;
thread_local unsigned failures = 0;

void beforeAllocation()
{
    if (remaining == 0)
    {
        ++failures;
        throw std::bad_alloc {};
    }
    if (remaining > 0)
        --remaining;
}

void freeAligned (void* pointer) noexcept
{
#if defined(_MSC_VER)
    _aligned_free (pointer);
#else
    std::free (pointer);
#endif
}
} // namespace allocation_failure

void* operator new (std::size_t bytes)
{
    allocation_failure::beforeAllocation();
    if (void* pointer = std::malloc (bytes == 0 ? 1 : bytes))
        return pointer;
    throw std::bad_alloc {};
}

void* operator new[] (std::size_t bytes) { return ::operator new (bytes); }
void operator delete (void* pointer) noexcept { std::free (pointer); }
void operator delete[] (void* pointer) noexcept { std::free (pointer); }

void* operator new (std::size_t bytes, std::align_val_t alignment)
{
    allocation_failure::beforeAllocation();
    const auto align = static_cast<std::size_t> (alignment);
    const auto size = bytes == 0 ? std::size_t { 1 } : bytes;
#if defined(_MSC_VER)
    if (void* pointer = _aligned_malloc (size, align))
        return pointer;
#else
    void* pointer = nullptr;
    if (posix_memalign (&pointer, align, size) == 0)
        return pointer;
#endif
    throw std::bad_alloc {};
}

void* operator new[] (std::size_t bytes, std::align_val_t alignment)
{ return ::operator new (bytes, alignment); }
void operator delete (void* pointer, std::align_val_t) noexcept
{ allocation_failure::freeAligned (pointer); }
void operator delete[] (void* pointer, std::align_val_t) noexcept
{ allocation_failure::freeAligned (pointer); }

#if defined(__cpp_sized_deallocation)
void operator delete (void* pointer, std::size_t) noexcept { std::free (pointer); }
void operator delete[] (void* pointer, std::size_t) noexcept { std::free (pointer); }
void operator delete (void* pointer, std::size_t, std::align_val_t) noexcept
{ allocation_failure::freeAligned (pointer); }
void operator delete[] (void* pointer, std::size_t, std::align_val_t) noexcept
{ allocation_failure::freeAligned (pointer); }
#endif

// Standalone deterministic stress regression, not an unbounded fuzz target.
// Link with TaikorDSP and run under ASan/UBSan (including float-cast-overflow).
// The seed and operation printed on failure reproduce the failing stream.
namespace
{
using Engine = taikor::TaikoEngine;
using Parameters = taikor::EngineParameters;
constexpr int capacity = 257;
constexpr float guard = 12345.0f;
std::uint32_t activeSeed = 0;
int activeStep = 0;

void require (bool condition, const char* message)
{
    if (! condition)
    {
        std::cerr << "Engine fuzz failure: seed=" << activeSeed
                  << " step=" << activeStep << " " << message << '\n';
        std::exit (1);
    }
}

struct Random
{
    std::uint32_t state;
    std::uint32_t next() noexcept
    {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return state;
    }
    int integer() noexcept { return std::bit_cast<std::int32_t> (next()); }
    float scalar() noexcept
    {
        constexpr std::array<std::uint32_t, 14> edges {{
            0u, 0x80000000u, 1u, 0x80000001u, 0x7f7fffffu, 0xff7fffffu,
            0x7f800000u, 0xff800000u, 0x7fc00000u, 0x3f800000u,
            0xbf800000u, 0x3f000000u, 0x40000000u, 0x00800000u
        }};
        const auto choice = next();
        return std::bit_cast<float> ((choice & 3u) == 0u
            ? next() : edges[(choice >> 2) % edges.size()]);
    }
};

constexpr auto fields = std::to_array<float Parameters::*> ({
    &Parameters::headDiameter, &Parameters::bodyDepth, &Parameters::tension,
    &Parameters::headMaterial, &Parameters::shellMaterial, &Parameters::resonantTension,
    &Parameters::cavityCoupling, &Parameters::headDamping, &Parameters::shellResonance,
    &Parameters::pitch, &Parameters::bachiHardness, &Parameters::strikePosition,
    &Parameters::strikeAzimuth, &Parameters::velocityDepth, &Parameters::velocityCurve,
    &Parameters::tensionModulation, &Parameters::strikeNoise, &Parameters::humanise,
    &Parameters::octaveBody, &Parameters::micDistance, &Parameters::micSpread,
    &Parameters::stereoWidth, &Parameters::drive, &Parameters::outputGain,
    &Parameters::outputHighPassHz, &Parameters::ensembleVariation,
    &Parameters::physicalFamilyMix, &Parameters::physicalRearTensionScale
});

void checkVisual (const Engine& engine)
{
    taikor::DrumVisualState state;
    engine.getVisualState (state);
    require (std::isfinite (state.strikeRadius) && state.strikeRadius >= 0.0f
                 && state.strikeRadius <= 1.0f, "strike radius");
    require (std::isfinite (state.strikeAngle) && std::isfinite (state.strikeLevel)
                 && std::isfinite (state.fundamentalHz), "nonfinite visual state");
    require (state.activeVoices >= 0 && state.activeVoices <= 16, "voice bound");
    require (state.lastOctaveOffset >= 0 && state.lastOctaveOffset < taikor::drumCount,
             "visual octave");
    require (static_cast<std::size_t> (state.lastArticulation) < taikor::articulationCount,
             "visual articulation");
    require (std::isfinite (engine.getOutputLevel (std::numeric_limits<int>::min()))
                 && std::isfinite (engine.getOutputLevel (std::numeric_limits<int>::max())),
             "output meters");
}

std::uint64_t stream (std::uint32_t seed, int iterations)
{
    activeSeed = seed;
    Random random { seed };
    Engine engine;
    Parameters parameters;
    std::uint64_t hash = 1469598103934665603ull;
    std::array<float, capacity + 2> left, right;
    std::array<float, capacity> extraLeft {}, extraRight {}, gain {};
    std::array<taikor::StereoPan, capacity> pan {};
    constexpr std::array<double, 10> rates {{
        std::numeric_limits<double>::quiet_NaN(), 8000.0, 44100.0, 48000.0,
        96000.0, 384000.0, 0.0, -1.0,
        std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()
    }};
    constexpr std::array<int, 10> blocks {{ -1, 0, 1, 2, 7, 31, 32, 33, 127, capacity }};

    // Exhaust all underlying enum values, including Count and invalid values
    // that could otherwise accidentally map onto another drum's MIDI note.
    for (int value = 0; value < 256; ++value)
    {
        const auto articulation = static_cast<taikor::Articulation> (value);
        const int note = taikor::midiNoteFor (articulation, random.integer());
        if (value < static_cast<int> (taikor::articulationCount))
            require (taikor::articulationForMidiNote (note) == articulation, "MIDI round trip");
        else
        {
            require (note == -1, "unsupported articulation MIDI mapping");
            engine.trigger (articulation, random.integer(), 1.0f);
        }
        require (! taikor::getArticulationDisplayName (articulation).empty(), "metadata fallback");
    }
    require (engine.getActiveVoiceCount() == 0, "invalid enums started voices");

    const int ratePeriod = std::max (1, iterations / static_cast<int> (rates.size()));
    for (activeStep = 0; activeStep < iterations; ++activeStep)
    {
        // Cover every listed rate and parameter field independently of PRNG
        // luck; rebuild only every fourth operation to keep CI cost bounded.
        if (activeStep % ratePeriod == 0)
            engine.prepare (rates[static_cast<std::size_t> (activeStep / ratePeriod)
                                   % rates.size()], random.integer());
        if (activeStep % 4 == 0)
        {
            if (activeStep % 20 == 0)
            {
                parameters = {};
                for (const auto field : fields)
                    parameters.*field = random.scalar();
            }
            parameters.*fields[static_cast<std::size_t> (activeStep / 4) % fields.size()]
                = random.scalar();
            parameters.performer = random.integer();
            parameters.ensembleSize = random.integer();
            parameters.physicalFamily = random.integer();
            const auto safe = Engine::sanitiseParameters (parameters);
            require (Engine::sanitiseParameters (safe) == safe, "sanitization not idempotent");
            for (const auto field : fields)
                require (std::isfinite (safe.*field), "sanitization left a nonfinite field");
            require (safe.performer >= 0 && safe.performer <= 3
                         && safe.ensembleSize >= 1 && safe.ensembleSize <= taikor::maximumEnsembleSize,
                     "sanitized player or ensemble range");
            require (safe.physicalFamily == 0 && safe.physicalFamilyMix == 0.0f
                         && safe.physicalRearTensionScale == 1.0f, "internal coordinates not reset");
            engine.setParameters (parameters);
        }
        switch (activeStep % 12)
        {
            case 0: engine.setEnsembleMember (random.integer()); break;
            case 1: engine.setHandDamping (random.scalar()); break;
            case 2: engine.setPitchBend (random.scalar()); break;
            case 3: engine.setStrikeAzimuthOverride (random.scalar()); break;
            case 4: engine.setStrikePositionOverride (random.scalar()); break;
            case 5: engine.setRearHeadStrike ((random.next() & 1u) != 0u); break;
            case 6: engine.clearStrikeOverrides(); break;
            case 7: engine.allSoundsOff(); break;
            case 8: engine.reset(); break;
            default: break;
        }
        const int midi = activeStep % 2 == 0 ? random.integer()
            : taikor::referenceNote + (activeStep % 4) * 12 + activeStep % 3;
        const bool mapped = taikor::articulationForMidiNote (midi).has_value();
        require (engine.triggerMidi (midi, random.scalar()) == mapped, "MIDI return contract");
        engine.trigger (static_cast<taikor::Articulation> (activeStep % 4),
                        random.integer(), activeStep % 3 == 0 ? 0.85f : random.scalar(),
                        random.scalar(), random.scalar());
        if (activeStep % 30 == 0)
        {
            const auto contact = Engine::measureContact (parameters,
                static_cast<taikor::Articulation> (random.next() & 255u),
                random.integer(), random.scalar());
            require (std::isfinite (contact) && contact > 0.0f, "contact readout domain");
        }
        // Dense simultaneous contacts exceed the voice pool without adding
        // a long audio render or unbounded event stream.
        if (activeStep == 39 || activeStep == 119)
            for (int hit = 0; hit < 24; ++hit)
                engine.trigger (static_cast<taikor::Articulation> (hit % 4), hit % 4, 0.9f);

        left.fill (guard);
        right.fill (guard);
        const int count = blocks[static_cast<std::size_t> (activeStep) % blocks.size()];
        for (int index = 0; index < capacity; ++index)
        {
            const auto slot = static_cast<std::size_t> (index);
            extraLeft[slot] = random.scalar();
            extraRight[slot] = random.scalar();
            gain[slot] = random.scalar();
            pan[slot] = { random.scalar(), random.scalar(), random.scalar(), random.scalar() };
        }
        engine.process (nullptr, right.data() + 1, capacity);
        engine.processRaw (left.data() + 1, nullptr, capacity);
        engine.processWithEnsemble (nullptr, nullptr, std::numeric_limits<int>::max(),
                                    nullptr, nullptr, nullptr, true);
        const bool raw = activeStep % 3 == 1;
        if (activeStep % 3 == 0)
            engine.process (left.data() + 1, right.data() + 1, count);
        else if (raw)
            engine.processRaw (left.data() + 1, right.data() + 1, count);
        else
            engine.processWithEnsemble (left.data() + 1, right.data() + 1, count,
                extraLeft.data(), activeStep % 5 == 0 ? nullptr : extraRight.data(),
                gain.data(), true, pan.data());

        require (left.front() == guard && right.front() == guard, "leading buffer guard");
        for (int index = 0; index < capacity + 1; ++index)
        {
            const auto slot = static_cast<std::size_t> (index + 1);
            if (index >= std::max (count, 0))
                require (left[slot] == guard && right[slot] == guard, "trailing buffer guard");
            else
            {
                require (std::isfinite (left[slot]) && std::isfinite (right[slot]), "nonfinite audio");
                if (! raw)
                    require (std::abs (left[slot]) <= 1.0f && std::abs (right[slot]) <= 1.0f,
                             "processed audio exceeded full scale");
                hash = (hash ^ std::bit_cast<std::uint32_t> (left[slot])) * 1099511628211ull;
                hash = (hash ^ std::bit_cast<std::uint32_t> (right[slot])) * 1099511628211ull;
            }
        }
        checkVisual (engine);
    }
    engine.allSoundsOff();
    engine.process (left.data(), right.data(), capacity);
    require (engine.getActiveVoiceCount() == 0, "panic retained contacts");
    for (int index = 0; index < capacity; ++index)
        require (left[static_cast<std::size_t> (index)] == 0.0f
                     && right[static_cast<std::size_t> (index)] == 0.0f, "panic not silent");
    return hash;
}

void allocationFailureRegression()
{
    activeSeed = 0;
    activeStep = 0;
    // Run before any successful measure() call initializes its thread-local
    // engines: fail the contact audit first, then the bank builder, then its
    // retry. No allocations are allowed for diagnostics while failure is armed.
    for (const int successfulAllocations : { 0, 1, 0 })
    {
        const unsigned before = allocation_failure::failures;
        allocation_failure::remaining = successfulAllocations;
        const auto result = Engine::measure ({}, 0);
        allocation_failure::remaining = -1;
        require (allocation_failure::failures == before + 1,
                 "measurement scratch allocation failure was not exercised");
        require (result.idealFundamentalHz == 0.0f && result.loadedFundamentalHz == 0.0f
                     && result.breathingModeHz == 0.0f && result.soundingHz == 0.0f
                     && result.tailSeconds == 0.0f, "allocation failure invented a pitch or tail");
        require (std::isfinite (result.radiusMetres) && std::isfinite (result.depthMetres)
                     && std::isfinite (result.tensionNewtonsPerMetre)
                     && std::isfinite (result.arealDensityKgPerSquareMetre)
                     && std::isfinite (result.waveSpeedMetresPerSecond)
                     && std::isfinite (result.headStiffnessParameter)
                     && std::isfinite (result.cavityStiffnessFactor)
                     && std::isfinite (result.cavityEnergyFraction),
                 "allocation failure returned nonfinite physical fields");
    }
    const auto recovered = Engine::measure ({}, 0);
    require (std::isfinite (recovered.soundingHz) && recovered.soundingHz > 0.0f
                 && recovered.idealFundamentalHz > 0.0f,
             "measurement scratch did not recover after allocation failure");
    // Both scratch engines are now initialized. Even with further allocations
    // forbidden, a repeat must return exactly the successful measurement.
    const unsigned before = allocation_failure::failures;
    allocation_failure::remaining = 0;
    const auto repeated = Engine::measure ({}, 0);
    allocation_failure::remaining = -1;
    require (allocation_failure::failures == before
                 && repeated.soundingHz == recovered.soundingHz
                 && repeated.idealFundamentalHz == recovered.idealFundamentalHz
                 && repeated.tailSeconds == recovered.tailSeconds,
             "initialized measurement allocated or changed its result");
    std::cout << "measurement_allocation_failures=" << allocation_failure::failures
              << " recovery=passed\n";
}

void partitionRegression()
{
    activeSeed = 0;
    activeStep = 0;
    Engine whole, split;
    whole.trigger (taikor::Articulation::Don, 0, 0.9f);
    split.trigger (taikor::Articulation::Don, 0, 0.9f);
    whole.setPitchBend (0.7f);
    split.setPitchBend (0.7f);
    std::array<float, 1024> left {}, right {}, splitLeft {}, splitRight {};
    whole.process (left.data(), right.data(), static_cast<int> (left.size()));
    int offset = 0;
    constexpr std::array<int, 7> sizes {{ 1, 31, 32, 7, 127, 33, 257 }};
    for (int block = 0; offset < static_cast<int> (left.size()); ++block)
    {
        const int count = std::min (sizes[static_cast<std::size_t> (block) % sizes.size()],
                                     static_cast<int> (left.size()) - offset);
        split.process (splitLeft.data() + offset, splitRight.data() + offset, 0);
        split.processRaw (splitLeft.data() + offset, splitRight.data() + offset, -1);
        split.process (splitLeft.data() + offset, nullptr, count);
        split.process (splitLeft.data() + offset, splitRight.data() + offset, count);
        offset += count;
    }
    require (left == splitLeft && right == splitRight, "block partition or empty call changed sound");
}
} // namespace

int main (int argc, char** argv)
{
    int iterations = 60;
    bool singleSeed = false;
    std::uint32_t selectedSeed = 0;
    for (int argument = 1; argument < argc; ++argument)
    {
        const std::string_view option = argv[argument];
        if ((option != "--seed" && option != "--iterations") || argument + 1 >= argc)
        {
            std::cerr << "Usage: EngineFuzzTests [--seed <nonzero uint32, decimal or 0xhex>]"
                         " [--iterations <1..1000000>]\n";
            return 2;
        }
        std::string_view value = argv[++argument];
        int base = 10;
        if (value.starts_with ("0x") || value.starts_with ("0X"))
        {
            value.remove_prefix (2);
            base = 16;
        }
        std::uint32_t parsed = 0;
        const auto result = std::from_chars (value.data(), value.data() + value.size(), parsed, base);
        if (result.ec != std::errc {} || result.ptr != value.data() + value.size() || parsed == 0u
            || (option == "--iterations" && parsed > 1000000u))
        {
            std::cerr << "Invalid value for " << option << '\n';
            return 2;
        }
        if (option == "--seed")
        {
            selectedSeed = parsed;
            singleSeed = true;
        }
        else
            iterations = static_cast<int> (parsed);
    }
    allocationFailureRegression();
    partitionRegression();
    const std::array<std::uint32_t, 3> seeds {{
        singleSeed ? selectedSeed : 0x12345678u, 0x9e3779b9u, 0xd15ea5e5u
    }};
    const int seedCount = singleSeed ? 1 : static_cast<int> (seeds.size());
    for (int index = 0; index < seedCount; ++index)
    {
        const auto seed = seeds[static_cast<std::size_t> (index)];
        const auto first = stream (seed, iterations);
        const auto replay = stream (seed, iterations);
        require (first == replay, "seed replay changed audio");
        std::cout << "seed=" << seed << " iterations=" << iterations
                  << " replays=2 digest=" << first << " passed\n";
    }
    std::cout << "Engine boundary fuzz passed: seeds=" << seedCount
              << " stream_operations=" << seedCount * 2 * iterations
              << " enum_probes=" << seedCount * 2 * 256
              << " partition_samples=1024\n";
    return 0;
}
