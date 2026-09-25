#include "PluginProcessor.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// Standalone, deterministic, bounded regression corpus. Run with --verbose to
// print the reproducer immediately before each operation (also useful under
// ASan/UBSan). No wall-clock or audio-performance assertions are made here.
// All state lengths describe real owned storage; MIDI uses JUCE's public raw
// insertion API, never a fabricated MidiBuffer header or dangling payload.
namespace
{
constexpr std::uint32_t defaultSeed = 0x5441494bu;
std::uint32_t seed = defaultSeed;
unsigned iterations = 16;
constexpr int preparedBlockSize = 128;
constexpr std::size_t maximumCorpusBytes = 4u * 1024u * 1024u + 1u;
constexpr float guardValue = 12345.25f;
int failures = 0;
int stateOperations = 0;
int audioOperations = 0;
int midiEvents = 0;
std::size_t prefixCases = 0;
unsigned randomStateCases = 0;
unsigned valueCases = 0;
unsigned treeCases = 0;
unsigned lifecycleCycles = 0;
unsigned malformedMidiCases = 0;
unsigned boundaryCases = 0;
unsigned randomMidiBlocks = 0;
bool verbose = false;
std::string context;

void beginCase (const std::string& label)
{
    context = label;
    if (verbose)
        std::cout << "CASE " << context << std::endl;
}

void expect (bool condition, const std::string& detail)
{
    if (! condition)
    {
        ++failures;
        if (failures <= 40)
            std::cerr << "FAIL [" << context << "]: " << detail << '\n';
    }
}

struct Random
{
    std::uint32_t state = seed;
    std::uint32_t next() noexcept
    {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return state;
    }
};

std::vector<juce::RangedAudioParameter*> parameterList (TaikorAudioProcessor& processor)
{
    std::vector<juce::RangedAudioParameter*> result;
    for (auto* parameter : processor.getParameters())
    {
        auto* ranged = dynamic_cast<juce::RangedAudioParameter*> (parameter);
        expect (ranged != nullptr, "non-ranged host parameter");
        if (ranged != nullptr)
            result.push_back (ranged);
    }
    expect (! result.empty(), "empty parameter layout");
    return result;
}

std::vector<float> snapshot (TaikorAudioProcessor& processor)
{
    std::vector<float> result;
    for (const auto* parameter : parameterList (processor))
    {
        const auto* raw = processor.parameters.getRawParameterValue (parameter->paramID);
        expect (raw != nullptr, "missing raw parameter");
        result.push_back (raw != nullptr ? raw->load (std::memory_order_relaxed)
                                       : std::numeric_limits<float>::quiet_NaN());
    }
    return result;
}

bool sameValues (const std::vector<float>& first, const std::vector<float>& second)
{
    if (first.size() != second.size())
        return false;
    for (std::size_t index = 0; index < first.size(); ++index)
        if (! std::isfinite (first[index]) || ! std::isfinite (second[index])
            || std::abs (first[index] - second[index]) > 0.0001f)
            return false;
    return true;
}

void checkParameters (TaikorAudioProcessor& processor)
{
    for (const auto* parameter : parameterList (processor))
    {
        const auto id = parameter->paramID.toStdString();
        const auto* raw = processor.parameters.getRawParameterValue (parameter->paramID);
        expect (raw != nullptr, id + " has no raw value");
        if (raw == nullptr)
            continue;
        const float value = raw->load (std::memory_order_relaxed);
        const auto& range = parameter->getNormalisableRange();
        const bool legal = std::isfinite (value) && value >= range.start
                                                  && value <= range.end;
        expect (legal, id + " raw value is nonfinite/out of range");
        const float normalised = parameter->getValue();
        expect (std::isfinite (normalised) && normalised >= 0.0f && normalised <= 1.0f,
                id + " host value is nonfinite/out of range");
        if (legal)
            expect (std::abs (range.snapToLegalValue (value) - value) <= 0.0001f,
                    id + " raw value is not snapped to its legal step");
    }
}

juce::MemoryBlock encode (const juce::ValueTree& tree)
{
    juce::MemoryBlock bytes;
    const auto xml = tree.createXml();
    expect (xml != nullptr, "corpus XML creation failed");
    if (xml != nullptr)
        juce::AudioProcessor::copyXmlToBinary (*xml, bytes);
    return bytes;
}

juce::MemoryBlock save (TaikorAudioProcessor& processor)
{
    juce::MemoryBlock bytes;
    processor.getStateInformation (bytes);
    expect (bytes.getSize() > 8u && bytes.getSize() <= maximumCorpusBytes,
            "saved state has implausible size");
    return bytes;
}

void restore (TaikorAudioProcessor& processor, const void* bytes, std::size_t length)
{
    expect (length <= maximumCorpusBytes, "corpus allocation exceeded bound");
    if (length > maximumCorpusBytes)
        return;
    ++stateOperations;
    processor.setStateInformation (bytes, static_cast<int> (length));
    checkParameters (processor);
}

void restore (TaikorAudioProcessor& processor, const juce::MemoryBlock& bytes)
{
    restore (processor, bytes.getData(), bytes.getSize());
}

void checkSavedState (TaikorAudioProcessor& processor)
{
    const auto before = snapshot (processor);
    const auto bytes = save (processor);
    const auto xml = juce::AudioProcessor::getXmlFromBinary (
        bytes.getData(), static_cast<int> (bytes.getSize()));
    expect (xml != nullptr, "saved state is not valid binary XML");
    if (xml == nullptr)
        return;
    const auto tree = juce::ValueTree::fromXml (*xml);
    expect (tree.hasType ("TAIKOR_STATE"), "saved state has wrong root");
    for (const auto* parameter : parameterList (processor))
    {
        int matches = 0;
        for (const auto& child : tree)
        {
            if (! child.hasType ("PARAM")
                || child.getProperty ("id").toString() != parameter->paramID)
                continue;
            ++matches;
            const auto valueText = child.getProperty ("value").toString();
            const char* start = valueText.toRawUTF8();
            char* end = nullptr;
            const double value = std::strtod (start, &end);
            const auto& range = parameter->getNormalisableRange();
            expect (end != start && *end == '\0' && std::isfinite (value)
                        && value >= static_cast<double> (range.start) - 0.0001
                        && value <= static_cast<double> (range.end) + 0.0001,
                    parameter->paramID.toStdString() + " saved an invalid value");
            const auto* raw = processor.parameters.getRawParameterValue (parameter->paramID);
            if (raw != nullptr)
                expect (std::abs (value - static_cast<double> (raw->load())) <= 0.0001,
                        parameter->paramID.toStdString() + " saved/raw values disagree");
        }
        expect (matches == 1, parameter->paramID.toStdString()
                                  + " saved missing/duplicate parameter nodes");
    }
    restore (processor, bytes);
    expect (sameValues (before, snapshot (processor)), "save/reload changed parameters");
}

// Both channels own a prefix and suffix guard. A zero-sample block still has
// valid channel pointers. Only the supported stereo output layout is used.
std::vector<float> render (TaikorAudioProcessor& processor, int samples,
                           juce::MidiBuffer midi)
{
    constexpr int guardSamples = 8;
    std::array<std::vector<float>, 2> storage;
    std::array<float*, 2> channels {};
    for (std::size_t channel = 0; channel < storage.size(); ++channel)
    {
        storage[channel].assign (static_cast<std::size_t> (samples + 2 * guardSamples),
                                 guardValue);
        channels[channel] = storage[channel].data() + guardSamples;
        std::fill_n (channels[channel], samples, 0.125f);
    }
    juce::AudioBuffer<float> buffer { channels.data(), 2, samples };
    ++audioOperations;
    midiEvents += midi.getNumEvents();
    processor.processBlock (buffer, midi);
    std::vector<float> result;
    for (std::size_t channel = 0; channel < storage.size(); ++channel)
    {
        for (int index = 0; index < guardSamples; ++index)
        {
            expect (storage[channel][static_cast<std::size_t> (index)] == guardValue,
                    "audio prefix guard overwritten");
            expect (storage[channel][static_cast<std::size_t> (guardSamples + samples + index)]
                        == guardValue, "audio suffix guard overwritten");
        }
        for (int index = 0; index < samples; ++index)
        {
            const float value = channels[channel][index];
            expect (std::isfinite (value), "nonfinite audio");
            result.push_back (value);
        }
        expect (std::isfinite (processor.getOutputLevel (static_cast<int> (channel))),
                "nonfinite output meter");
    }
    expect (processor.getActiveVoiceCount() >= 0, "negative active voice count");
    return result;
}

void setNormalised (TaikorAudioProcessor& processor, const juce::String& id, float value)
{
    auto* parameter = processor.parameters.getParameter (id);
    expect (parameter != nullptr, "missing setup parameter");
    if (parameter != nullptr)
        parameter->setValueNotifyingHost (value);
}

void configureMidi (TaikorAudioProcessor& processor)
{
    setNormalised (processor, taikor::parameters::humanise, 0.0f);
    setNormalised (processor, taikor::parameters::ensembleSize, 0.0f);
    setNormalised (processor, taikor::parameters::reverbRoom, 0.0f);
    setNormalised (processor, taikor::parameters::reverbMix, 0.0f);
    processor.prepareToPlay (48000.0, preparedBlockSize);
}

void stateBytes()
{
    std::cout << "State bytes: seed=" << seed << std::endl;
    auto processor = std::make_unique<TaikorAudioProcessor>();
    beginCase ("state baseline");
    setNormalised (*processor, taikor::parameters::pitch, 0.75f);
    const auto baseline = save (*processor);
    const auto baselineValues = snapshot (*processor);
    beginCase ("null empty state");
    restore (*processor, nullptr, 0u);
    expect (sameValues (baselineValues, snapshot (*processor)), "empty state changed parameters");
    const auto* source = static_cast<const std::uint8_t*> (baseline.getData());
    for (std::size_t length = 0; length < baseline.getSize(); ++length)
    {
        ++prefixCases;
        beginCase ("state prefix length=" + std::to_string (length));
        // Exact-sized ownership makes ASan useful at the actual truncation edge.
        std::vector<std::uint8_t> prefix (source, source + length);
        restore (*processor, prefix.data(), prefix.size());
    }
    Random random;
    for (unsigned index = 0; index < iterations; ++index)
    {
        randomStateCases += 2u;
        beginCase ("state random index=" + std::to_string (index));
        std::vector<std::uint8_t> bytes (1u + random.next() % 2048u);
        for (auto& byte : bytes)
            byte = static_cast<std::uint8_t> (random.next());
        restore (*processor, bytes.data(), bytes.size());
        beginCase ("state mutation index=" + std::to_string (index));
        bytes.assign (source, source + baseline.getSize());
        for (unsigned flip = 0; flip < 1u + index % 4u; ++flip)
            bytes[random.next() % bytes.size()] ^= static_cast<std::uint8_t> (
                1u << (random.next() % 8u));
        restore (*processor, bytes.data(), bytes.size());
    }
    for (const std::uint32_t claimed : { 0u, 1u, 0x7fffffffu, 0xffffffffu })
    {
        beginCase ("state header declared length=" + std::to_string (claimed));
        std::vector<std::uint8_t> bytes (source, source + baseline.getSize());
        if (bytes.size() >= 8u)
            for (unsigned byte = 0; byte < 4u; ++byte)
                bytes[4u + byte] = static_cast<std::uint8_t> (claimed >> (8u * byte));
        restore (*processor, bytes.data(), bytes.size());
    }
    beginCase ("oversized owned garbage");
    std::vector<std::uint8_t> oversized (maximumCorpusBytes, 0xa5u);
    restore (*processor, oversized.data(), oversized.size());
    beginCase ("oversized valid XML with unknown payload");
    juce::ValueTree oversizedTree { "TAIKOR_STATE" };
    oversizedTree.setProperty ("unknownPayload",
                              juce::String::repeatedString ("x", 4 * 1024 * 1024 - 512),
                              nullptr);
    const auto encoded = encode (oversizedTree);
    expect (encoded.getSize() <= maximumCorpusBytes, "oversized XML corpus overflow");
    restore (*processor, encoded);
    beginCase ("known-good recovery after byte corpus");
    restore (*processor, baseline);
    expect (sameValues (baselineValues, snapshot (*processor)), "valid restore failed after garbage");
    checkSavedState (*processor);
}

void stateValues()
{
    std::cout << "State values and tree structure" << std::endl;
    auto processor = std::make_unique<TaikorAudioProcessor>();
    beginCase ("value corpus baseline");
    const auto defaults = snapshot (*processor);
    const auto baselineTree = processor->parameters.copyState();
    const auto baseline = encode (baselineTree);
    processor->prepareToPlay (48000.0, preparedBlockSize);
    const std::array<const char*, 19> values {
        "nan", "NaN", "inf", "-inf", "1e309", "-1e309", "1.7976931348623157e308",
        "-1.7976931348623157e308", "3.4028234663852886e38", "-3.4028234663852886e38",
        "1e-300", "-1e-300", "", "garbage", "0.5junk", " 0.5 ", "-0", "0.37", "2.5"
    };
    for (const auto* parameter : parameterList (*processor))
    {
        for (std::size_t index = 0; index < values.size(); ++index)
        {
            ++valueCases;
            beginCase ("state value id=" + parameter->paramID.toStdString()
                       + " index=" + std::to_string (index) + " text=" + values[index]);
            auto tree = baselineTree.createCopy();
            auto child = tree.getChildWithProperty ("id", parameter->paramID);
            expect (child.isValid(), "missing corpus parameter child");
            if (! child.isValid())
                continue;
            child.setProperty ("value", juce::String (values[index]), nullptr);
            restore (*processor, encode (tree));
            checkSavedState (*processor);
            // Exercise raw-value consumers, including reverb integer conversion.
            render (*processor, 1, {});
        }
    }
    for (int shape = 0; shape < 8; ++shape)
    {
        ++treeCases;
        beginCase ("state tree shape=" + std::to_string (shape));
        auto tree = baselineTree.createCopy();
        auto first = tree.getChild (0);
        if (shape == 0)
            first.removeProperty ("value", nullptr);
        else if (shape == 1)
            first.removeProperty ("id", nullptr);
        else if (shape == 2)
        {
            auto duplicate = first.createCopy();
            duplicate.setProperty ("value", "nan", nullptr);
            tree.appendChild (duplicate, nullptr);
        }
        else if (shape == 3)
        {
            juce::ValueTree unknown { "PARAM" };
            unknown.setProperty ("id", "futureParameter", nullptr);
            unknown.setProperty ("value", "nan", nullptr);
            tree.appendChild (unknown, nullptr);
        }
        else if (shape == 4)
        {
            tree = juce::ValueTree { "TAIKOR_STATE" };
            auto parent = tree;
            for (int depth = 0; depth < 64; ++depth)
            {
                juce::ValueTree nested { "NESTED" };
                parent.appendChild (nested, nullptr);
                parent = nested;
            }
        }
        else if (shape == 5)
        {
            tree = juce::ValueTree { "WRONG_ROOT" };
            tree.appendChild (first.createCopy(), nullptr);
        }
        else if (shape == 6)
            first.appendChild (juce::ValueTree { "NESTED" }, nullptr);
        else
            for (int index = 0; index < 512; ++index)
                tree.appendChild (first.createCopy(), nullptr);

        const auto bytes = encode (tree);
        restore (*processor, baseline);
        restore (*processor, bytes);
        const auto firstResult = snapshot (*processor);
        checkSavedState (*processor);
        restore (*processor, baseline);
        restore (*processor, bytes);
        expect (sameValues (firstResult, snapshot (*processor)), "tree restore was nondeterministic");
        render (*processor, 17, {});
    }
    beginCase ("legacy empty state fills defaults");
    for (auto* parameter : parameterList (*processor))
        parameter->setValueNotifyingHost (1.0f);
    restore (*processor, encode (juce::ValueTree { "TAIKOR_STATE" }));
    expect (sameValues (defaults, snapshot (*processor)), "missing parameters inherited previous preset");
    checkSavedState (*processor);
    processor->releaseResources();
}

void lifecycle()
{
    std::cout << "Restore/save/reprepare lifecycle" << std::endl;
    auto processor = std::make_unique<TaikorAudioProcessor>();
    Random random;
    beginCase ("save before prepare");
    checkSavedState (*processor);
    constexpr std::array<double, 3> rates { 44100.0, 48000.0, 96000.0 };
    constexpr std::array<int, 6> sizes { 0, 1, 2, 17, 127, 128 };
    for (unsigned iteration = 0; iteration < std::min (iterations, 12u); ++iteration)
    {
        ++lifecycleCycles;
        beginCase ("lifecycle iteration=" + std::to_string (iteration));
        for (auto* parameter : parameterList (*processor))
        {
            const float normalised = iteration % 3 == 0 ? 0.0f
                                   : iteration % 3 == 1 ? 1.0f
                                   : static_cast<float> (random.next() % 1001u) / 1000.0f;
            parameter->setValueNotifyingHost (normalised);
        }
        checkParameters (*processor);
        const auto expected = snapshot (*processor);
        const auto bytes = save (*processor);
        restore (*processor, bytes);
        processor->prepareToPlay (rates[static_cast<std::size_t> (iteration) % rates.size()],
                                  preparedBlockSize);
        expect (processor->isEngineReady(), "prepare did not publish readiness");
        for (const int size : sizes)
        {
            juce::MidiBuffer midi;
            if (size > 0)
                midi.addEvent (juce::MidiMessage::noteOn (
                    1, taikor::midiNoteFor (taikor::Articulation::Don, 0), 0.75f), size - 1);
            render (*processor, size, midi);
        }
        checkSavedState (*processor);
        processor->releaseResources();
        expect (! processor->isEngineReady(), "release left engine ready");
        restore (*processor, bytes);
        expect (sameValues (expected, snapshot (*processor)), "lifecycle changed saved parameters");
        checkSavedState (*processor);
    }
}

void addRaw (juce::MidiBuffer& midi, const std::vector<juce::uint8>& bytes, int position)
{
    expect (! bytes.empty(), "empty raw MIDI corpus entry");
    if (! bytes.empty())
        expect (midi.addEvent (bytes.data(), static_cast<int> (bytes.size()), position),
                "JUCE rejected raw MIDI insertion");
}

void midiCorpus()
{
    std::cout << "Malformed/mixed MIDI and guarded buffer boundaries" << std::endl;
    const auto note = static_cast<juce::uint8> (
        taikor::midiNoteFor (taikor::Articulation::Don, 0));
    const std::vector<std::vector<juce::uint8>> ignored {
        { 0x90 }, { 0x90, note }, { 0x80 }, { 0x80, note },
        { 0xb0 }, { 0xb0, 120 }, { 0xe0 }, { 0xe0, 127 },
        { 0x90, static_cast<juce::uint8> (note | 0x80u), 127 },
        { 0x90, note, 0xff }, { 0x90, note, 0x80 },
        { 0xb0, 0x90, 127 }, { 0xb0, 16, 0xff },
        { 0xb0, 17, 0xff }, { 0xb0, 18, 0xff },
        { 0xe0, 0xff, 127 }, { 0xe0, 0, 0xff },
        { 0x90, note, 0 }, { 0x80, note, 127 },
        { 0xf0, 0x7d, 0x01, 0xf7 }, { 0xf0, 0x7d, 0x01 },
        { 0xf8 }, { 0xfa }, { 0xfb }, { 0xfc }, { 0xfe },
        { 0xc0, 1 }, { 0xd0, 127 }, { 0xa0, note, 127 }
    };
    auto processor = std::make_unique<TaikorAudioProcessor>();
    beginCase ("MIDI preparation");
    configureMidi (*processor);
    for (std::size_t index = 0; index < ignored.size(); ++index)
    {
        ++malformedMidiCases;
        beginCase ("ignored MIDI packet=" + std::to_string (index));
        processor->requestPanic();
        const auto before = processor->getTriggerCounter (taikor::Articulation::Don);
        juce::MidiBuffer midi;
        addRaw (midi, ignored[index], 0);
        const auto audio = render (*processor, 17, midi);
        expect (processor->getTriggerCounter (taikor::Articulation::Don) == before,
                "ignored/malformed packet triggered a note");
        expect (std::all_of (audio.begin(), audio.end(), [] (float value) { return value == 0.0f; }),
                "ignored/malformed packet produced audio from silence");
    }
    // Compare identical fresh engines with and without inserted junk. This
    // catches ignored packets changing controller state even when silent.
    processor->releaseResources();
    processor.reset();
    auto clean = std::make_unique<TaikorAudioProcessor>();
    auto dirty = std::make_unique<TaikorAudioProcessor>();
    beginCase ("mixed MIDI deterministic comparison");
    configureMidi (*clean);
    configureMidi (*dirty);
    juce::MidiBuffer cleanMidi;
    juce::MidiBuffer dirtyMidi;
    for (const auto& packet : ignored)
        addRaw (dirtyMidi, packet, 0);
    for (auto* midi : { &cleanMidi, &dirtyMidi })
        midi->addEvent (juce::MidiMessage::noteOn (1, note, 0.75f), 17);
    auto reference = render (*clean, 128, cleanMidi);
    auto candidate = render (*dirty, 128, dirtyMidi);
    bool audible = false;
    for (int block = 0; block < 8; ++block)
    {
        expect (reference.size() == candidate.size(), "mixed MIDI changed buffer size");
        for (std::size_t index = 0; index < std::min (reference.size(), candidate.size()); ++index)
        {
            expect (std::abs (reference[index] - candidate[index]) <= 1.0e-7f,
                    "ignored packets changed subsequent valid audio");
            audible = audible || std::abs (reference[index]) > 1.0e-8f;
        }
        reference = render (*clean, 128, {});
        candidate = render (*dirty, 128, {});
    }
    expect (audible, "mixed MIDI comparison never rendered a valid stroke");
    dirty->releaseResources();
    dirty.reset();

    for (const int size : { 0, 1, 2, 17, 127, 128 })
    {
        for (const int position : { -1, 0, std::max (0, size - 1), size, size + 1,
                                    std::numeric_limits<int>::max() })
        {
            ++boundaryCases;
            beginCase ("MIDI boundary samples=" + std::to_string (size)
                       + " position=" + std::to_string (position));
            clean->requestPanic();
            juce::MidiBuffer midi;
            midi.addEvent (juce::MidiMessage::noteOn (16, note, 0.5f), position);
            render (*clean, size, midi);
            // Out-of-block timestamp policy is deliberately unspecified here;
            // safe iteration and buffer ownership are the regression contract.
        }
    }
    beginCase ("in-block MIDI must not sound before its timestamp");
    clean->requestPanic();
    juce::MidiBuffer delayed;
    delayed.addEvent (juce::MidiMessage::noteOn (1, note, 0.8f), 63);
    const auto delayedAudio = render (*clean, 128, delayed);
    for (std::size_t channel = 0; channel < 2u; ++channel)
        for (std::size_t sample = 0; sample < 63u; ++sample)
            expect (delayedAudio[channel * 128u + sample] == 0.0f,
                    "note sounded before its event timestamp");

    beginCase ("same-sample note/panic ordering");
    juce::MidiBuffer noteThenPanic;
    noteThenPanic.addEvent (juce::MidiMessage::noteOn (1, note, 1.0f), 0);
    noteThenPanic.addEvent (juce::MidiMessage::allSoundOff (1), 0);
    const auto stopped = render (*clean, 128, noteThenPanic);
    expect (std::all_of (stopped.begin(), stopped.end(), [] (float value) { return value == 0.0f; }),
            "same-sample all-sound-off did not cancel preceding note");

    Random random;
    for (unsigned block = 0; block < iterations; ++block)
    {
        ++randomMidiBlocks;
        beginCase ("seeded mixed MIDI block=" + std::to_string (block));
        const int size = static_cast<int> (random.next() % 129u);
        juce::MidiBuffer midi;
        for (unsigned event = 0; event < 16u; ++event)
        {
            const int position = static_cast<int> (random.next() % static_cast<unsigned> (size + 3)) - 1;
            const int channel = 1 + static_cast<int> (random.next() % 16u);
            switch (random.next() % 4u)
            {
                case 0: midi.addEvent (juce::MidiMessage::noteOn (channel, note, 0.5f), position); break;
                case 1: midi.addEvent (juce::MidiMessage::controllerEvent (
                            channel, static_cast<int> (random.next() % 128u),
                            static_cast<int> (random.next() % 128u)), position); break;
                case 2: midi.addEvent (juce::MidiMessage::pitchWheel (
                            channel, static_cast<int> (random.next() % 16384u)), position); break;
                default: addRaw (midi, ignored[random.next() % ignored.size()], position); break;
            }
        }
        render (*clean, size, midi);
    }
    clean->releaseResources();
}
} // namespace

int main (int argc, char** argv)
{
    const auto usage = []
    {
        std::cerr << "Usage: TaikorStateFuzzTests [--verbose] "
                     "[--case all|state-bytes|state-values|lifecycle|midi] "
                     "[--iterations 1..256] [--seed 1..4294967295]\n"
                     "Defaults: 16 iterations, seed 0x5441494b. Hex seeds accepted. "
                     "Fixed regression cases always run; lifecycle caps at 12 cycles.\n";
    };
    std::string selected = "all";
    for (int index = 1; index < argc; ++index)
    {
        const std::string argument { argv[index] };
        if (argument == "--verbose")
            verbose = true;
        else if (argument == "--case" && index + 1 < argc)
            selected = argv[++index];
        else if ((argument == "--iterations" || argument == "--seed") && index + 1 < argc)
        {
            const std::string text { argv[++index] };
            const bool hex = text.size() > 2u && text[0] == '0'
                             && (text[1] == 'x' || text[1] == 'X');
            const char* start = text.c_str();
            char* end = nullptr;
            errno = 0;
            const auto number = std::strtoull (start, &end, hex ? 16 : 10);
            const auto maximum = argument == "--iterations" ? 256ull : 4294967295ull;
            if (text.empty() || text[0] < '0' || text[0] > '9' || end == start
                || *end != '\0' || errno == ERANGE || number == 0ull || number > maximum)
            {
                usage();
                return 2;
            }
            if (argument == "--iterations")
                iterations = static_cast<unsigned> (number);
            else
                seed = static_cast<std::uint32_t> (number);
        }
        else
        {
            usage();
            return 2;
        }
    }
    if (selected != "all" && selected != "state-bytes" && selected != "state-values"
        && selected != "lifecycle" && selected != "midi")
    {
        std::cerr << "Unknown case: " << selected << '\n';
        return 2;
    }
    juce::ScopedJuceInitialiser_GUI guiInitialiser;
    if (selected == "all" || selected == "state-bytes") stateBytes();
    if (selected == "all" || selected == "state-values") stateValues();
    if (selected == "all" || selected == "lifecycle") lifecycle();
    if (selected == "all" || selected == "midi") midiCorpus();
    std::cout << "seed=" << seed << " restores=" << stateOperations
              << " iterations=" << iterations << " audio-blocks=" << audioOperations
              << " midi-events=" << midiEvents << " prefixes=" << prefixCases
              << " random-state-cases=" << randomStateCases << " parameter-value-cases=" << valueCases
              << " tree-shapes=" << treeCases << " lifecycle-cycles=" << lifecycleCycles
              << " malformed-midi-packets=" << malformedMidiCases << " boundary-cases=" << boundaryCases
              << " random-midi-blocks=" << randomMidiBlocks << " failures=" << failures << '\n';
    return failures == 0 ? 0 : 1;
}
