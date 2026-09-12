// Full callback and isolated convolution probes. Use the same options and
// compare captures with BenchmarkCompare.py before accepting a CPU change.
#include "PluginProcessor.h"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using Clock = std::chrono::steady_clock;
struct Options
{
    std::string selected = "all";
    int rate = 48000, block = 256, repeats = 3;
    double seconds = 3.0, tail = 1.0, warmup = 0.25;
    std::uint32_t seed = 1;
    std::filesystem::path capture;
};

Options parse (int argc, char** argv)
{
    Options result;
    for (int i = 1; i < argc; ++i)
    {
        const std::string key = argv[i];
        if (key == "--help")
        {
            std::cout << "TaikorBenchmarkAudioPath --case all|plugin-idle|plugin-dense|"
                "plugin-controls|room-off|room-hall|room-theater|room-opera|"
                "room-switch|room-tail|room-silence --rate 48000|96000 --block 64|256 "
                "--seconds N --tail N --repeats N --warmup N --seed N --capture DIR\n";
            std::exit (0);
        }
        if (++i == argc) throw std::runtime_error ("Missing value for " + key);
        const std::string value = argv[i];
        if (key == "--case") result.selected = value;
        else if (key == "--capture") result.capture = value;
        else
        {
            std::size_t consumed = 0;
            const double n = std::stod (value, &consumed);
            if (consumed != value.size() || ! std::isfinite (n) || n < 0.0)
                throw std::runtime_error ("Invalid number for " + key);
            if (key == "--seconds") result.seconds = n;
            else if (key == "--tail") result.tail = n;
            else if (key == "--warmup") result.warmup = n;
            else
            {
                if (n > 1000000.0 || n != std::floor (n))
                    throw std::runtime_error ("Invalid integer for " + key);
                const int number = static_cast<int> (n);
                if (key == "--rate") result.rate = number;
                else if (key == "--block") result.block = number;
                else if (key == "--repeats") result.repeats = number;
                else if (key == "--seed") result.seed = static_cast<std::uint32_t> (number);
                else throw std::runtime_error ("Unknown option: " + key);
            }
        }
    }
    if ((result.rate != 48000 && result.rate != 96000)
        || (result.block != 64 && result.block != 256)
        || result.repeats < 1 || result.repeats > 100
        || result.seconds > 600.0 || std::llround (result.seconds * result.rate) < 1
        || result.tail >= result.seconds || result.warmup > 60.0)
        throw std::runtime_error ("Options out of range");
    return result;
}

void parameter (TaikorAudioProcessor& processor, const char* id, float value)
{
    auto* p = processor.parameters.getParameter (id);
    if (p == nullptr) throw std::runtime_error (std::string ("Missing parameter ") + id);
    p->setValueNotifyingHost (p->convertTo0to1 (value));
}

struct Take
{
    std::vector<float> audio;
    std::vector<double> blocks;
    double elapsed = 0.0, peak = 0.0, rms = 0.0, maxBlockPercent = 0.0;
    std::uint64_t hash = 14695981039346656037ull;
    int processCalls = 0, eventBlocks = 0;
};

Take render (const Options& options, const std::string& name, double seconds)
{
    const juce::ScopedNoDenormals noDenormals;
    const bool plugin = name.starts_with ("plugin-");
    const bool controls = name == "plugin-controls";
    const bool silent = name == "plugin-idle" || name == "room-silence";
    auto processor = plugin ? std::make_unique<TaikorAudioProcessor>() : nullptr;
    auto reverb = plugin ? nullptr : std::make_unique<taikor::IRReverb>();
    if (plugin)
    {
        parameter (*processor, taikor::parameters::performer,
                   static_cast<float> (options.seed % 4));
        parameter (*processor, taikor::parameters::reverbRoom, controls ? 1.0f : 0.0f);
        processor->prepareToPlay (options.rate, options.block);
    }
    else
    {
        const int room = name == "room-off" ? 0 : name == "room-theater" ? 2
            : name == "room-opera" ? 3 : 1;
        reverb->prepare (options.rate, options.block, room, 0.35f);
    }

    const int frames = static_cast<int> (std::llround (seconds * options.rate));
    const int playing = static_cast<int> (std::llround (
        (options.seconds - options.tail) * options.rate));
    Take take;
    take.audio.resize (static_cast<std::size_t> (2 * frames));
    take.blocks.reserve (static_cast<std::size_t> ((frames + options.block - 1) / options.block));
    juce::AudioBuffer<float> buffer (2, options.block);
    juce::MidiBuffer midi;
    midi.ensureSize (8192);
    std::uint32_t random = options.seed;
    int stroke = 0, phase = -1;
    const int spacing = options.rate / 16;
    const int onset = controls ? std::min (options.rate / 8, playing / 10) : 0;
    for (int start = 0; start < frames; start += options.block)
    {
        const int count = std::min (options.block, frames - start);
        buffer.setSize (2, count, false, false, true);
        buffer.clear();
        midi.clear();
        if (plugin && ! silent)
        {
            for (int sample = start; sample < start + count && sample < playing; ++sample)
                if (sample >= onset && (sample - onset) % spacing == 0)
                {
                    random = random * 1664525u + 1013904223u;
                    const int velocity = 70 + static_cast<int> ((random >> 16) % 51);
                    if (stroke == 0)
                        for (int drum = 0; drum < 4; ++drum)
                            midi.addEvent (juce::MidiMessage::noteOn (
                                1, 48 + 12 * drum, static_cast<juce::uint8> (velocity)), sample - start);
                    else
                        midi.addEvent (juce::MidiMessage::noteOn (
                            1, 48 + 12 * (stroke % 4) + stroke % 12,
                            static_cast<juce::uint8> (velocity)), sample - start);
                    ++stroke;
                }
        }
        if (! plugin && ! silent)
            for (int sample = 0; sample < count; ++sample)
            {
                const int absolute = start + sample;
                if (absolute >= playing || (name == "room-tail" && absolute != 0)) continue;
                random = random * 1664525u + 1013904223u;
                const float noise = static_cast<float> (random >> 8) / 16777215.0f - 0.5f;
                const float envelope = std::exp (-40.0f * static_cast<float> (absolute % spacing)
                                                / static_cast<float> (options.rate));
                buffer.setSample (0, sample, 0.3f * noise * envelope);
                buffer.setSample (1, sample, -0.21f * noise * envelope);
            }
        const int nextPhase = std::min (9, static_cast<int> (
            10.0 * start / std::max (1, playing)));
        const bool change = nextPhase != phase;
        const auto before = Clock::now();
        if (controls && change)
        {
            parameter (*processor, taikor::parameters::reverbRoom, static_cast<float> (nextPhase % 4));
            parameter (*processor, taikor::parameters::reverbMix, nextPhase == 5 ? 0.0f : 0.6f);
            parameter (*processor, taikor::parameters::drive, nextPhase % 3 == 0 ? 0.4f : 0.0f);
            parameter (*processor, taikor::parameters::ensembleSize, nextPhase >= 3 && nextPhase < 7 ? 4.0f : 1.0f);
            if (nextPhase == 6)
                parameter (*processor, taikor::parameters::tension, 0.66f);
            midi.addEvent (juce::MidiMessage::pitchWheel (1, nextPhase % 2 ? 9200 : 8192), count / 3);
            midi.addEvent (juce::MidiMessage::controllerEvent (1, 1, nextPhase % 2 ? 5 : 0), count / 2);
            if (nextPhase == 8)
                midi.addEvent (juce::MidiMessage::allSoundOff (1), count / 2);
        }
        if (plugin)
            processor->processBlock (buffer, midi);
        else
        {
            if (name == "room-switch" && change)
                reverb->setParameters (nextPhase % 4, nextPhase == 5 ? 0.0f : 0.6f);
            else if (name != "room-switch")
                reverb->setParameters (name == "room-off" ? 0 : name == "room-theater" ? 2
                    : name == "room-opera" ? 3 : 1, 0.35f);
            reverb->process (buffer.getWritePointer (0), buffer.getWritePointer (1), count);
        }
        const double elapsed = std::chrono::duration<double> (Clock::now() - before).count();
        take.elapsed += elapsed;
        take.blocks.push_back (elapsed);
        take.maxBlockPercent = std::max (take.maxBlockPercent, elapsed * options.rate / count * 100.0);
        ++take.processCalls;
        if (! midi.isEmpty() || ((controls || name == "room-switch") && change)) ++take.eventBlocks;
        phase = nextPhase;
        for (int sample = 0; sample < count; ++sample)
            for (int channel = 0; channel < 2; ++channel)
                take.audio[static_cast<std::size_t> (2 * (start + sample) + channel)] =
                    buffer.getSample (channel, sample);
    }
    double energy = 0.0;
    for (float value : take.audio)
    {
        if (! std::isfinite (value)) throw std::runtime_error ("Non-finite audio in " + name);
        take.peak = std::max (take.peak, std::abs (static_cast<double> (value)));
        energy += static_cast<double> (value) * value;
        const auto bits = std::bit_cast<std::uint32_t> (value);
        for (int byte = 0; byte < 4; ++byte)
        {
            take.hash ^= (bits >> (8 * byte)) & 0xffu;
            take.hash *= 1099511628211ull;
        }
    }
    take.rms = std::sqrt (energy / static_cast<double> (take.audio.size()));
    if (silent && take.peak != 0.0) throw std::runtime_error ("Nonzero idle output");
    if (! silent && frames > onset && take.peak == 0.0)
        throw std::runtime_error ("Unexpectedly silent " + name);
    return take;
}

void benchmark (const Options& options, const std::string& name)
{
    if (options.warmup > 0.0)
        (void) render (options, name, std::min (options.seconds, options.warmup));
    std::vector<Take> takes;
    for (int repeat = 0; repeat < options.repeats; ++repeat)
    {
        auto take = render (options, name, options.seconds);
        if (! takes.empty() && (take.hash != takes.front().hash || take.audio != takes.front().audio))
            throw std::runtime_error ("Nondeterministic " + name);
        takes.push_back (std::move (take));
    }
    if (! options.capture.empty())
    {
        std::filesystem::create_directories (options.capture);
        const auto path = options.capture / (name + "-" + std::to_string (options.rate)
            + "-b" + std::to_string (options.block) + "-s" + std::to_string (options.seed) + ".f32le");
        std::ofstream stream (path, std::ios::binary);
        for (float value : takes.front().audio)
        {
            const auto bits = std::bit_cast<std::uint32_t> (value);
            for (int byte = 0; byte < 4; ++byte)
                stream.put (static_cast<char> ((bits >> (8 * byte)) & 0xffu));
        }
        stream.close();
        if (! stream) throw std::runtime_error ("Cannot write " + path.string());
    }
    std::sort (takes.begin(), takes.end(), [] (const auto& a, const auto& b) { return a.elapsed < b.elapsed; });
    auto& median = takes[takes.size() / 2];
    const double duration = static_cast<double> (median.audio.size()) / (2.0 * options.rate);
    const double elapsed = takes.size() % 2 ? median.elapsed
        : 0.5 * (median.elapsed + takes[takes.size() / 2 - 1].elapsed);
    std::sort (median.blocks.begin(), median.blocks.end());
    std::cout << name << ',' << options.rate << ',' << options.block << ','
        << median.audio.size() / 2 << ',' << options.repeats << ',' << options.seed << ','
        << std::fixed << std::setprecision (6) << elapsed / duration * 100.0 << ','
        << takes.front().elapsed / duration * 100.0 << ',' << takes.back().elapsed / duration * 100.0 << ','
        << median.blocks[median.blocks.size() / 2] * 1.0e6 << ','
        << median.blocks[static_cast<std::size_t> (std::ceil (0.95 * (median.blocks.size() - 1)))] * 1.0e6 << ','
        << median.blocks.back() * 1.0e6 << ',' << median.maxBlockPercent << ','
        << median.peak << ',' << median.rms << ',' << std::hex << median.hash << std::dec << ','
        << median.processCalls << ',' << median.eventBlocks << '\n' << std::flush;
}
} // namespace

int main (int argc, char** argv)
{
    const juce::ScopedJuceInitialiser_GUI juce;
    try
    {
        const auto options = parse (argc, argv);
        const std::vector<std::string> cases { "plugin-idle", "plugin-dense", "plugin-controls",
            "room-off", "room-hall", "room-theater", "room-opera", "room-switch", "room-tail", "room-silence" };
        if (options.selected != "all" && std::find (cases.begin(), cases.end(), options.selected) == cases.end())
            throw std::runtime_error ("Unknown case: " + options.selected);
        std::cout << "case,rate_hz,block,frames,repeats,seed,median_cpu_percent,min_cpu_percent,max_cpu_percent,"
            "median_block_us,p95_block_us,max_block_us,max_block_cpu_percent,audio_peak,audio_rms,audio_fnv1a64,"
            "process_calls,event_blocks\n";
        for (const auto& name : cases)
            if (options.selected == "all" || options.selected == name) benchmark (options, name);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "TaikorBenchmarkAudioPath: " << error.what() << '\n';
        return 1;
    }
}
