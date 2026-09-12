// CPU/audio regression protocol: build Release for the host architecture, save
// this executable before changing DSP, and run baseline/candidate alternately
// on an otherwise idle machine with identical arguments. Compare median CPU
// percentages (elapsed DSP seconds / rendered seconds * 100), including trigger
// and live control costs. Preparation, allocation, validation and file I/O are
// excluded. Each timed host block includes all sample-accurate event splits.
// Block percentiles describe the repetition whose total time is the median;
// the maximum includes initial strikes and is not a scheduling guarantee.
// A short discarded warmup precedes fresh, identically prepared repetitions.
// Fresh instances prevent a previous control sweep's retained caches or poles
// from influencing the next take; construction/preparation remain untimed.
//
// Defaults render 2 s of playing plus 1 s of tail. Solo is one sustained open
// stroke; dense/ensembles begin with all four drums, then mix all four gestures
// at 16 hits/s. Controls additionally changes pitch bend, palm damping, strike
// position, microphone/drive controls, drum tension and ensemble size while
// sound is running, after changing controls during an initial silent interval.
// Engine reset fixes its internal gesture/noise sequence; --seed fixes score
// velocities and the performer identity. Every repetition must be finite,
// audible (except idle) and byte-identical to the first. Capture files contain headerless,
// interleaved stereo IEEE-754 float32, little-endian, without normalisation.
// Compare captures with cmp, or decode them as <f4 to quantify differences.
// Example:
//   TaikorBenchmarkCPU --case dense --repeats 5 --capture /tmp/taikor-before
//   TaikorBenchmarkCPU --case controls --seconds 14 --tail 12 --repeats 1

#include "DSP/EnsembleEngine.h"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#if defined (__SSE2__) || defined (_M_X64)
 #include <xmmintrin.h>
#endif

namespace
{
using Clock = std::chrono::steady_clock;

// Match the plug-in's ScopedNoDenormals on x86/SSE and GCC/Clang AArch64.
// Other targets retain their default floating-point environment.
class ScopedNoDenormals
{
public:
    ScopedNoDenormals()
    {
#if defined (__SSE2__) || defined (_M_X64)
        previous = _mm_getcsr();
        _mm_setcsr (static_cast<unsigned int> (previous) | 0x8040u);
#elif defined (__aarch64__) && (defined (__GNUC__) || defined (__clang__))
        __asm__ __volatile__ ("mrs %0, fpcr" : "=r" (previous));
        const auto next = previous | (std::uint64_t { 1 } << 24);
        __asm__ __volatile__ ("msr fpcr, %0" : : "r" (next));
#endif
    }
    ~ScopedNoDenormals()
    {
#if defined (__SSE2__) || defined (_M_X64)
        _mm_setcsr (static_cast<unsigned int> (previous));
#elif defined (__aarch64__) && (defined (__GNUC__) || defined (__clang__))
        __asm__ __volatile__ ("msr fpcr, %0" : : "r" (previous));
#endif
    }
private:
    [[maybe_unused]] std::uint64_t previous { 0 };
};

struct Options
{
    std::string selectedCase { "all" };
    int rate { 48000 }, block { 256 }, repeats { 3 };
    double seconds { 3.0 }, tail { 1.0 }, warmup { 0.25 };
    std::uint32_t seed { 1 };
    std::filesystem::path capture;
};

void usage()
{
    std::cout << "Usage: TaikorBenchmarkCPU [options]\n"
        "  --case idle|solo|dense|ensemble4|ensemble8|controls|all (default all)\n"
        "  --rate 48000|96000       Sample rate (default 48000)\n"
        "  --block 64|256          Host block size (default 256)\n"
        "  --seconds N             Total duration including tail (default 3)\n"
        "  --tail N                Final unstruck duration (default 1)\n"
        "  --repeats N             Timed repetitions (default 3)\n"
        "  --warmup N              Discarded warmup seconds (default 0.25)\n"
        "  --seed N                Unsigned score/performer seed (default 1)\n"
        "  --capture DIRECTORY     Save first repetition as stereo .f32le\n"
        "  --help                  Show this help\n"
        "CSV output; CPU is percent of one real-time audio thread. Captures\n"
        "are headerless interleaved float32 LE at --rate, without gain trims.\n"
        "Match every option across baseline/candidate and alternate runs.\n";
}

double number (const std::string& value)
{
    std::size_t consumed = 0;
    const double parsed = std::stod (value, &consumed);
    if (consumed != value.size() || ! std::isfinite (parsed))
        throw std::runtime_error ("Invalid number: " + value);
    return parsed;
}

int integer (const std::string& value)
{
    const double parsed = number (value);
    if (parsed < 0.0 || parsed > 1000000.0 || std::floor (parsed) != parsed)
        throw std::runtime_error ("Invalid integer: " + value);
    return static_cast<int> (parsed);
}

Options parse (int argc, char** argv)
{
    Options result;
    for (int index = 1; index < argc; ++index)
    {
        const std::string key = argv[index];
        if (key == "--help")
        {
            usage();
            std::exit (0);
        }
        if (index + 1 >= argc)
            throw std::runtime_error ("Missing value for " + key);
        const std::string value = argv[++index];
        if (key == "--case") result.selectedCase = value;
        else if (key == "--rate") result.rate = integer (value);
        else if (key == "--block") result.block = integer (value);
        else if (key == "--repeats") result.repeats = integer (value);
        else if (key == "--seconds") result.seconds = number (value);
        else if (key == "--tail") result.tail = number (value);
        else if (key == "--warmup") result.warmup = number (value);
        else if (key == "--capture") result.capture = value;
        else if (key == "--seed")
        {
            const double parsed = number (value);
            if (parsed < 0 || parsed > std::numeric_limits<std::uint32_t>::max()
                || std::floor (parsed) != parsed)
                throw std::runtime_error ("Seed must be an unsigned 32-bit integer");
            result.seed = static_cast<std::uint32_t> (parsed);
        }
        else throw std::runtime_error ("Unknown option: " + key);
    }
    if ((result.rate != 48000 && result.rate != 96000)
        || (result.block != 64 && result.block != 256)
        || result.repeats < 1 || result.repeats > 100
        || result.seconds <= 0.0 || result.seconds > 600.0
        || std::llround (result.seconds * result.rate) < 1
        || result.tail < 0.0 || result.tail >= result.seconds
        || result.warmup < 0.0 || result.warmup > 60.0
        || (result.warmup > 0.0
            && std::llround (std::min (result.seconds, result.warmup) * result.rate) < 1))
        throw std::runtime_error ("Options out of range; see --help");
    return result;
}

enum class EventKind { hit, bend, palm, position, mix, tension, idleMix, size };
struct Event
{
    std::int64_t frame;
    EventKind kind;
    int note { 48 };
    float value { 0.0f };
};

std::vector<Event> score (const Options& options, const std::string& name)
{
    std::vector<Event> events;
    if (name == "idle") return events;
    const double playing = options.seconds - options.tail;
    const double onset = name == "controls" ? std::min (0.125, playing * 0.1) : 0.0;
    const auto frame = [&] (double time)
    { return static_cast<std::int64_t> (std::llround (time * options.rate)); };
    std::uint32_t random = options.seed;
    const auto velocity = [&random] ()
    {
        random = random * 1664525u + 1013904223u;
        return 0.55f + 0.4f * static_cast<float> (random >> 8) / 16777215.0f;
    };
    events.push_back ({ frame (onset), EventKind::hit, 48, velocity() });
    if (name != "solo")
    {
        for (int drum = 1; drum < 4; ++drum)
            events.push_back ({ frame (onset), EventKind::hit, 48 + 12 * drum, velocity() });
        int stroke = 0;
        for (double time = onset + 0.0625; time < playing; time += 0.0625, ++stroke)
            events.push_back ({ frame (time), EventKind::hit,
                48 + 12 * (stroke % 4) + ((stroke / 4 + stroke) % 4), velocity() });
    }
    if (name == "controls")
    {
        events.push_back ({ 0, EventKind::bend, 0, -0.2f });
        events.push_back ({ 0, EventKind::palm, 0, 0.03f });
        events.push_back ({ 0, EventKind::idleMix, 0, 0.0f });
        events.push_back ({ frame (playing * 0.15), EventKind::bend, 0, 0.13f });
        events.push_back ({ frame (playing * 0.25), EventKind::position, 0, -0.35f });
        events.push_back ({ frame (playing * 0.35), EventKind::palm, 0, 0.07f });
        events.push_back ({ frame (playing * 0.40), EventKind::size, 0, 4.0f });
        events.push_back ({ frame (playing * 0.45), EventKind::mix, 0, 0.63f });
        events.push_back ({ frame (playing * 0.55), EventKind::bend, 0, -0.08f });
        events.push_back ({ frame (playing * 0.65), EventKind::tension, 0, 0.66f });
        events.push_back ({ frame (playing * 0.75), EventKind::palm, 0, 0.0f });
        events.push_back ({ frame (playing * 0.80), EventKind::size, 0, 1.0f });
        events.push_back ({ frame (playing * 0.85), EventKind::bend, 0, 0.0f });
    }
    std::stable_sort (events.begin(), events.end(), [] (const auto& a, const auto& b)
    { return a.frame < b.frame; });
    return events;
}

struct Take
{
    std::vector<float> left, right;
    std::vector<double> blockSeconds;
    double elapsed { 0.0 }, maxBlockPercent { 0.0 }, peak { 0.0 }, rms { 0.0 };
    std::uint64_t hash { 14695981039346656037ull };
};

Take render (taikor::EnsembleEngine& engine, taikor::EngineParameters parameters,
             const Options& options, const std::vector<Event>& events,
             double seconds, bool expectSilence = false)
{
    engine.setParameters (parameters);
    engine.setHandDamping (0.0f);
    engine.setPitchBend (0.0f);
    engine.reset();
    Take take;
    const auto frames = static_cast<std::int64_t> (std::llround (seconds * options.rate));
    take.left.resize (static_cast<std::size_t> (frames));
    take.right.resize (static_cast<std::size_t> (frames));
    take.blockSeconds.reserve (static_cast<std::size_t> (
        (frames + options.block - 1) / options.block));
    std::size_t event = 0;
    for (std::int64_t start = 0; start < frames; start += options.block)
    {
        const auto end = std::min (frames, start + options.block);
        auto cursor = start;
        const auto before = Clock::now();
        while (cursor < end)
        {
            while (event < events.size() && events[event].frame <= cursor)
            {
                const auto& next = events[event++];
                switch (next.kind)
                {
                    case EventKind::hit:
                        if (! engine.triggerMidi (next.note, next.value))
                            throw std::runtime_error ("Benchmark scheduled an invalid note");
                        break;
                    case EventKind::bend: engine.setPitchBend (next.value); break;
                    case EventKind::palm: engine.setHandDamping (next.value); break;
                    case EventKind::position:
                        engine.setStrikePositionOverride (next.value);
                        engine.setStrikeAzimuthOverride (0.7f);
                        break;
                    case EventKind::mix:
                        parameters.drive = next.value;
                        parameters.outputHighPassHz = 120.0f;
                        parameters.micSpread = 0.7f;
                        engine.setParameters (parameters);
                        break;
                    case EventKind::tension:
                        parameters.tension = next.value;
                        engine.setParameters (parameters);
                        break;
                    case EventKind::idleMix:
                        parameters.outputGain = 0.05f;
                        parameters.stereoWidth = 0.85f;
                        parameters.drive = 0.15f;
                        engine.setParameters (parameters);
                        break;
                    case EventKind::size:
                        parameters.ensembleSize = static_cast<int> (next.value);
                        engine.setParameters (parameters);
                        break;
                }
            }
            const auto until = event < events.size()
                ? std::min (end, events[event].frame) : end;
            engine.process (take.left.data() + cursor, take.right.data() + cursor,
                            static_cast<int> (until - cursor));
            cursor = until;
        }
        const double elapsed = std::chrono::duration<double> (Clock::now() - before).count();
        take.elapsed += elapsed;
        take.blockSeconds.push_back (elapsed);
        take.maxBlockPercent = std::max (take.maxBlockPercent,
            elapsed * options.rate / static_cast<double> (end - start) * 100.0);
    }
    double energy = 0.0;
    for (std::size_t index = 0; index < take.left.size(); ++index)
        for (float value : { take.left[index], take.right[index] })
        {
            if (! std::isfinite (value))
                throw std::runtime_error ("Non-finite audio");
            take.peak = std::max (take.peak, std::abs (static_cast<double> (value)));
            energy += static_cast<double> (value) * value;
            const auto bits = std::bit_cast<std::uint32_t> (value);
            for (int byte = 0; byte < 4; ++byte)
            {
                take.hash ^= (bits >> (8 * byte)) & 0xffu;
                take.hash *= 1099511628211ull;
            }
        }
    take.rms = std::sqrt (energy / (2.0 * static_cast<double> (frames)));
    if (! expectSilence && take.peak <= 1.0e-12)
        throw std::runtime_error ("Audio was unexpectedly silent");
    if (expectSilence && take.peak != 0.0)
        throw std::runtime_error ("Idle audio was unexpectedly nonzero");
    return take;
}

void capture (const std::filesystem::path& path, const Take& take)
{
    static_assert (sizeof (float) == 4 && std::numeric_limits<float>::is_iec559);
    std::ofstream output (path, std::ios::binary);
    if (! output) throw std::runtime_error ("Cannot create " + path.string());
    for (std::size_t frame = 0; frame < take.left.size(); ++frame)
        for (float value : { take.left[frame], take.right[frame] })
        {
            const auto bits = std::bit_cast<std::uint32_t> (value);
            for (int byte = 0; byte < 4; ++byte)
                output.put (static_cast<char> ((bits >> (8 * byte)) & 0xffu));
        }
    output.close();
    if (! output) throw std::runtime_error ("Cannot finish " + path.string());
}

double percentile (const std::vector<double>& sorted, double fraction)
{
    return sorted[static_cast<std::size_t> (std::ceil (
        fraction * static_cast<double> (sorted.size() - 1)))];
}

void benchmark (const Options& options, const std::string& name)
{
    taikor::EngineParameters parameters;
    parameters.performer = static_cast<int> (options.seed % 4);
    parameters.ensembleSize = name == "ensemble8" ? 8 : name == "ensemble4" ? 4 : 1;
    const auto makeEngine = [&]
    {
        auto engine = std::make_unique<taikor::EnsembleEngine>();
        engine->setParameters (parameters);
        engine->prepare (options.rate, options.block);
        return engine;
    };
    const auto events = score (options, name);
    if (options.warmup > 0.0)
    {
        const auto firstHit = std::find_if (events.begin(), events.end(), [] (const auto& event)
        { return event.kind == EventKind::hit; });
        const double warmupSeconds = std::min (options.seconds, options.warmup);
        const bool warmupSilent = firstHit == events.end()
            || firstHit->frame >= std::llround (warmupSeconds * options.rate);
        auto engine = makeEngine();
        (void) render (*engine, parameters, options, events,
                       warmupSeconds, warmupSilent);
    }
    std::vector<Take> repetitions;
    for (int repeat = 0; repeat < options.repeats; ++repeat)
    {
        auto engine = makeEngine();
        auto take = render (*engine, parameters, options, events, options.seconds, name == "idle");
        if (! repetitions.empty() && (take.hash != repetitions.front().hash
            || take.left != repetitions.front().left || take.right != repetitions.front().right))
            throw std::runtime_error ("Audio differs between reset repetitions of " + name);
        repetitions.push_back (std::move (take));
    }
    if (! options.capture.empty())
    {
        std::filesystem::create_directories (options.capture);
        capture (options.capture / (name + "-" + std::to_string (options.rate)
            + "-b" + std::to_string (options.block) + "-s" + std::to_string (options.seed)
            + ".f32le"), repetitions.front());
    }
    std::sort (repetitions.begin(), repetitions.end(), [] (const auto& a, const auto& b)
    { return a.elapsed < b.elapsed; });
    auto& middle = repetitions[repetitions.size() / 2];
    const double duration = static_cast<double> (middle.left.size()) / options.rate;
    const auto cpu = [duration] (double elapsed) { return elapsed / duration * 100.0; };
    const double medianElapsed = repetitions.size() % 2 == 0
        ? 0.5 * (middle.elapsed + repetitions[repetitions.size() / 2 - 1].elapsed)
        : middle.elapsed;
    std::sort (middle.blockSeconds.begin(), middle.blockSeconds.end());
    std::cout << name << ',' << options.rate << ',' << options.block << ','
        << middle.left.size() << ',' << options.repeats << ',' << options.seed << ','
        << std::fixed << std::setprecision (6)
        << cpu (medianElapsed) << ',' << cpu (repetitions.front().elapsed) << ','
        << cpu (repetitions.back().elapsed) << ','
        << percentile (middle.blockSeconds, 0.5) * 1.0e6 << ','
        << percentile (middle.blockSeconds, 0.95) * 1.0e6 << ','
        << middle.blockSeconds.back() * 1.0e6 << ',' << middle.maxBlockPercent << ','
        << middle.peak << ',' << middle.rms << ',' << std::hex << middle.hash
        << std::dec << '\n' << std::flush;
}
} // namespace

int main (int argc, char** argv)
{
    try
    {
        const auto options = parse (argc, argv);
        const ScopedNoDenormals noDenormals;
        const std::vector<std::string> cases { "idle", "solo", "dense", "ensemble4", "ensemble8", "controls" };
        if (options.selectedCase != "all"
            && std::find (cases.begin(), cases.end(), options.selectedCase) == cases.end())
            throw std::runtime_error ("Unknown case: " + options.selectedCase);
        std::cout << "case,rate_hz,block,frames,repeats,seed,median_cpu_percent,"
            "min_cpu_percent,max_cpu_percent,median_block_us,p95_block_us,max_block_us,"
            "max_block_cpu_percent,audio_peak,audio_rms,audio_fnv1a64\n";
        for (const auto& name : cases)
            if (options.selectedCase == "all" || options.selectedCase == name)
                benchmark (options, name);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "TaikorBenchmarkCPU: " << error.what() << '\n';
        return 1;
    }
}
