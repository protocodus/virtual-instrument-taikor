#include "PluginProcessor.h"
#include "PluginEditor.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <locale>
#include <memory>
#include <sstream>
#include <string_view>
#include <vector>

namespace
{
namespace ids = taikor::parameters;

// Index into the processor's raw-pointer cache. Kept in the same order as the
// layout below so the two cannot drift apart.
enum Slot
{
    slotHeadDiameter = 0,
    slotBodyDepth,
    slotTension,
    slotHeadMaterial,
    slotShellMaterial,
    slotResonantTension,
    slotCavityCoupling,
    slotHeadDamping,
    slotShellResonance,
    slotPitch,
    slotBachiHardness,
    slotStrikePosition,
    slotVelocityDepth,
    slotTensionModulation,
    slotStrikeNoise,
    slotHumanise,
    slotOctaveBody,
    slotMicDistance,
    slotMicSpread,
    slotStereoWidth,
    slotDrive,
    slotOutput,
    slotStrikeAzimuth,
    slotPerformer,
    slotVelocityCurve,
    slotOutputHighPass,
    slotEnsembleSize,
    slotEnsembleVariation,
    slotReverbRoom,
    slotReverbMix,
    slotCount
};

static_assert (static_cast<int> (slotCount) == ids::parameterCount,
               "the slot table and the declared parameter count must agree");

constexpr std::array<const char*, slotCount> parameterIds {
    ids::headDiameter, ids::bodyDepth, ids::tension, ids::headMaterial,
    ids::shellMaterial, ids::resonantTension, ids::cavityCoupling,
    ids::headDamping, ids::shellResonance, ids::pitch, ids::bachiHardness,
    ids::strikePosition, ids::velocityDepth, ids::tensionModulation,
    ids::strikeNoise, ids::humanise, ids::octaveBody, ids::micDistance,
    ids::micSpread, ids::stereoWidth, ids::drive, ids::output,
    ids::strikeAzimuth, ids::performer, ids::velocityCurve, ids::outputHighPass,
    ids::ensembleSize, ids::ensembleVariation, ids::reverbRoom, ids::reverbMix
};

float bipolarControllerValue (int rawValue) noexcept
{
    // MIDI's conventional centre is the exact value 64: 64 steps below it,
    // 63 above it.
    return rawValue < 64
        ? static_cast<float> (rawValue - 64) / 64.0f
        : static_cast<float> (rawValue - 64) / 63.0f;
}

std::unique_ptr<juce::RangedAudioParameter> makePercentParameter (
    const juce::String& id, const juce::String& name, float defaultValue,
    int versionHint = 1)
{
    return std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { id, versionHint }, name,
        juce::NormalisableRange<float> { 0.0f, 1.0f, 0.001f }, defaultValue,
        juce::AudioParameterFloatAttributes()
            .withLabel ("%")
            .withStringFromValueFunction ([] (float value, int)
            {
                return juce::String (juce::roundToInt (value * 100.0f));
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                return text.retainCharacters ("0123456789.-").getFloatValue() / 100.0f;
            }));
}

std::unique_ptr<juce::RangedAudioParameter> makeCentimetreParameter (
    const juce::String& id, const juce::String& name, float minimum, float maximum,
    float defaultValue, float step)
{
    return std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { id, 1 }, name,
        juce::NormalisableRange<float> { minimum, maximum, step }, defaultValue,
        juce::AudioParameterFloatAttributes()
            .withLabel ("cm")
            .withStringFromValueFunction ([] (float value, int)
            {
                return juce::String (value, 1);
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                return text.retainCharacters ("0123456789.-").getFloatValue();
            }));
}

// Thirty scalar parameters need only a few KiB. Bound bytes, recursion and
// node count before JUCE's recursive XML parser sees untrusted host state.
constexpr int maximumStateBytes = 1024 * 1024;

bool hasSafeXmlStructure (std::string_view text) noexcept
{
    int depth = 0, nodes = 0;
    for (std::size_t index = 0; index < text.size(); ++index)
    {
        if (text[index] == '\0')
            return false;
        if (text[index] != '<')
            continue;
        if (index + 1 == text.size() || text[index + 1] == '!')
            return false; // No DTDs, entities, CDATA or declarations in our format.
        const bool instruction = text[index + 1] == '?';
        const bool closing = text[index + 1] == '/';
        char quote = 0;
        auto end = index + 1;
        for (; end < text.size(); ++end)
        {
            const char c = text[end];
            if (c == '\0')
                return false;
            if (quote != 0)
            {
                if (c == quote)
                    quote = 0;
            }
            else if (c == '\'' || c == '"')
                quote = c;
            else if (c == '>')
                break;
            else if (c == '<')
                return false;
        }
        if (end == text.size())
            return false;
        if (instruction)
        {
            if (text[end - 1] != '?')
                return false;
        }
        else if (closing)
        {
            if (--depth < 0)
                return false;
        }
        else
        {
            if (++nodes > 512)
                return false;
            if (text[end - 1] != '/' && ++depth > 8)
                return false;
        }
        index = end;
    }
    return depth == 0 && nodes != 0;
}

bool readFiniteDecimal (const juce::String& input, double& value)
{
    const auto text = input.trim();
    if (text.isEmpty() || text.length() > 128)
        return false;
    auto cursor = text.getCharPointer();
    if (*cursor == '+' || *cursor == '-')
        ++cursor;
    const auto digits = [&cursor]
    {
        int count = 0;
        while (*cursor >= '0' && *cursor <= '9')
        {
            ++cursor;
            ++count;
        }
        return count;
    };
    int count = digits();
    if (*cursor == '.')
    {
        ++cursor;
        count += digits();
    }
    if (count == 0)
        return false;
    if (*cursor == 'e' || *cursor == 'E')
    {
        ++cursor;
        if (*cursor == '+' || *cursor == '-')
            ++cursor;
        if (digits() == 0)
            return false;
    }
    if (*cursor != 0)
        return false;
    // JUCE's text conversion accumulates the exponent in a signed int. State
    // can contain arbitrarily large exponent digits within our byte limit.
    // Classic-locale extraction handles overflow without integer wraparound
    // and correctly rounds representable values at the double boundary.
    std::istringstream stream (text.toStdString());
    stream.imbue (std::locale::classic());
    stream >> value;
    return ! stream.fail() && std::isfinite (value);
}

constexpr bool isValidArticulation (taikor::Articulation articulation) noexcept
{
    return static_cast<std::size_t> (articulation) < taikor::articulationCount;
}
} // namespace

TaikorAudioProcessor::TaikorAudioProcessor()
    : AudioProcessor (BusesProperties().withOutput ("Output",
                                                    juce::AudioChannelSet::stereo(),
                                                    true)),
      parameters (*this, nullptr, "TAIKOR_STATE", createParameterLayout())
{
    for (int slot = 0; slot < slotCount; ++slot)
    {
        parameterPointers[static_cast<std::size_t> (slot)] =
            parameters.getRawParameterValue (parameterIds[static_cast<std::size_t> (slot)]);
        jassert (parameterPointers[static_cast<std::size_t> (slot)] != nullptr);
        const auto* parameter = parameters.getParameter (parameterIds[static_cast<std::size_t> (slot)]);
        const auto& range = parameter->getNormalisableRange();
        parameterBounds[static_cast<std::size_t> (slot)] = {
            range.start, range.end,
            parameter->convertFrom0to1 (parameter->getDefaultValue())
        };
    }
}

juce::AudioProcessorValueTreeState::ParameterLayout
TaikorAudioProcessor::createParameterLayout()
{
    std::vector<std::unique_ptr<juce::RangedAudioParameter>> result;
    result.reserve (static_cast<std::size_t> (ids::parameterCount));

    // --- The drum -------------------------------------------------------
    result.push_back (makeCentimetreParameter (
        ids::headDiameter, "Head Diameter", ids::minimumDiameterCentimetres,
        ids::maximumDiameterCentimetres, 150.0f, 0.5f));
    result.push_back (makePercentParameter (ids::bodyDepth, "Body Depth", 0.5f));
    result.push_back (makePercentParameter (ids::tension, "Head Tension", 0.62f));
    result.push_back (makePercentParameter (ids::headMaterial, "Head Material", 0.75f));
    result.push_back (makePercentParameter (ids::shellMaterial, "Shell Material", 0.8f));
    result.push_back (makePercentParameter (
        ids::resonantTension, "Resonant Head", 0.5f));
    result.push_back (makePercentParameter (ids::cavityCoupling, "Air Coupling", 0.85f));
    result.push_back (makePercentParameter (ids::headDamping, "Head Damping", 0.50f));
    result.push_back (makePercentParameter (
        ids::shellResonance, "Shell Resonance", 0.4f));

    result.push_back (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ids::pitch, 1 }, "Pitch",
        juce::NormalisableRange<float> { -24.0f, 24.0f, 0.1f }, 0.0f,
        juce::AudioParameterFloatAttributes()
            .withLabel ("st")
            .withStringFromValueFunction ([] (float value, int)
            {
                return juce::String (value, 1);
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                return text.retainCharacters ("0123456789.-").getFloatValue();
            })));

    // --- The stroke -----------------------------------------------------
    result.push_back (makePercentParameter (
        ids::bachiHardness, "Bachi Hardness", 0.7f));

    result.push_back (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ids::strikePosition, 1 }, "Strike Position",
        juce::NormalisableRange<float> { -1.0f, 1.0f, 0.01f }, 0.0f,
        juce::AudioParameterFloatAttributes()
            .withStringFromValueFunction ([] (float value, int)
            {
                const auto amount = juce::roundToInt (std::abs (value) * 100.0f);
                if (amount == 0)
                    return juce::String ("As written");
                return (value < 0.0f ? juce::String ("Centre ") : juce::String ("Rim "))
                     + juce::String (amount);
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                const auto trimmed = text.trim().toUpperCase();
                const auto amount =
                    trimmed.retainCharacters ("0123456789.-").getFloatValue() / 100.0f;
                return trimmed.startsWith ("CENTRE") || trimmed.startsWith ("CENTER")
                    ? -std::abs (amount) : amount;
            })));

    result.push_back (makePercentParameter (
        ids::velocityDepth, "Velocity Depth", 0.75f));
    result.push_back (makePercentParameter (
        ids::tensionModulation, "Tension Mod", 0.4f));
    result.push_back (makePercentParameter (ids::strikeNoise, "Strike Noise", 0.35f));
    result.push_back (makePercentParameter (ids::humanise, "Humanise", 0.4f));

    result.push_back (std::make_unique<juce::AudioParameterChoice> (
        // Keep the established ID, order and normalised endpoints so old
        // projects restore, but advertise a real indexed choice to every host.
        // A legacy value below/above the midpoint snaps to the nearest endpoint.
        juce::ParameterID { ids::octaveBody, 1 }, "Drum Layout",
        juce::StringArray { "1 Drum", "4 Drums" }, 1,
        juce::AudioParameterChoiceAttributes()
            .withStringFromValueFunction ([] (int index, int)
            {
                return index == 0 ? juce::String ("1 Drum")
                                  : juce::String ("4 Drums");
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                const auto trimmed = text.trim().toUpperCase();
                if (trimmed.startsWith ("SINGLE")
                    || trimmed.startsWith ("TUNED")
                    || trimmed.startsWith ("1 DRUM"))
                    return 0;
                if (trimmed.startsWith ("FOUR")
                    || trimmed.startsWith ("FAMILY")
                    || trimmed.startsWith ("4 DRUM"))
                    return 1;
                const float number =
                    text.retainCharacters ("0123456789.-").getFloatValue();
                const float normalised = std::abs (number) > 1.0f
                    ? number * 0.01f : number;
                return normalised < 0.5f ? 0 : 1;
            })));

    // --- The close pair and the output ----------------------------------
    result.push_back (makeCentimetreParameter (
        ids::micDistance, "Mic Distance", ids::minimumMicDistanceCentimetres,
        ids::maximumMicDistanceCentimetres, 16.0f, 0.5f));
    result.push_back (makePercentParameter (ids::micSpread, "Mic Spread", 0.55f));
    result.push_back (makePercentParameter (ids::stereoWidth, "Stereo Width", 0.5f));
    result.push_back (makePercentParameter (ids::drive, "Drive", 0.0f));

    result.push_back (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ids::output, 1 }, "Output",
        // The model's fixed post-Drive reference calibration owns single-hit
        // headroom. Keep this public range and default stable so existing host
        // automation retains its exact normalised curve.
        juce::NormalisableRange<float> { -24.0f, 6.0f, 0.1f }, -22.5f,
        juce::AudioParameterFloatAttributes()
            .withLabel ("dB")
            .withStringFromValueFunction ([] (float value, int)
            {
                return juce::String (value, 1);
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                return text.retainCharacters ("0123456789.-").getFloatValue();
            })));

    // Appended after the established host slots so every earlier automation
    // index and saved parameter ID keeps its meaning.
    result.push_back (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ids::strikeAzimuth, 1 }, "Strike Azimuth",
        juce::NormalisableRange<float> { -180.0f, 180.0f, 0.1f }, 0.0f,
        juce::AudioParameterFloatAttributes()
            .withLabel (juce::String::fromUTF8 ("\xc2\xb0"))
            .withStringFromValueFunction ([] (float value, int)
            {
                return juce::String (value, 1);
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                return text.retainCharacters ("0123456789.-").getFloatValue();
            })));

    result.push_back (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { ids::performer, 1 }, "Performer",
        juce::StringArray { "P1", "P2", "P3", "P4" }, 0,
        juce::AudioParameterChoiceAttributes().withAutomatable (false)));

    result.push_back (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ids::velocityCurve, 1 }, "Velocity Curve",
        juce::NormalisableRange<float> { -1.0f, 1.0f, 0.01f }, 0.0f,
        juce::AudioParameterFloatAttributes()
            .withStringFromValueFunction ([] (float value, int)
            {
                const auto amount = juce::roundToInt (std::abs (value) * 100.0f);
                if (amount == 0)
                    return juce::String ("Linear");
                return (value < 0.0f ? juce::String ("Soft ")
                                     : juce::String ("Hard "))
                     + juce::String (amount);
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                const auto trimmed = text.trim().toUpperCase();
                if (trimmed.startsWith ("LINEAR"))
                    return 0.0f;
                const auto amount = std::abs (
                    trimmed.retainCharacters ("0123456789.-").getFloatValue())
                                  * 0.01f;
                return trimmed.startsWith ("SOFT") ? -amount : amount;
            })));

    // Append the new ID and use a higher version hint: AU hosts sort by this
    // before the ID hash, preserving every existing automation index.
    juce::NormalisableRange<float> highPassRange { 0.0f, 500.0f, 1.0f };
    highPassRange.setSkewForCentre (100.0f);
    result.push_back (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { ids::outputHighPass, 2 }, "Output High Pass",
        highPassRange, 0.0f,
        juce::AudioParameterFloatAttributes()
            .withStringFromValueFunction ([] (float value, int)
            {
                const auto frequency = juce::roundToInt (value);
                return frequency == 0 ? juce::String ("Off")
                                      : juce::String (frequency) + " Hz";
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                return text.retainCharacters ("0123456789.-").getFloatValue();
            })));

    result.push_back (std::make_unique<juce::AudioParameterInt> (
        juce::ParameterID { ids::ensembleSize, 3 }, "Ensemble Size", 1, 8, 1,
        juce::AudioParameterIntAttributes()
            .withStringFromValueFunction ([] (int value, int)
            {
                return value == 1 ? juce::String ("1 player")
                                  : juce::String (value) + " players";
            })
            .withValueFromStringFunction ([] (const juce::String& text)
            {
                return text.getIntValue();
            })));
    result.push_back (makePercentParameter (
        ids::ensembleVariation, "Ensemble Variation", 0.4f, 3));

    // New controls stay after the existing slots, including AU's version sort.
    result.push_back (std::make_unique<juce::AudioParameterChoice> (
        juce::ParameterID { ids::reverbRoom, 4 }, "Reverb Room",
        juce::StringArray { "Off", "Hall", "Theater", "Opera" }, 0));
    result.push_back (makePercentParameter (ids::reverbMix, "Reverb Dry/Wet", 0.2f, 4));

    jassert (static_cast<int> (result.size()) == ids::parameterCount);
    return { result.begin(), result.end() };
}

float TaikorAudioProcessor::readParameter (int slot) const noexcept
{
    const auto index = static_cast<std::size_t> (slot);
    const auto& bounds = parameterBounds[index];
    const auto* pointer = parameterPointers[index];
    const float value = pointer != nullptr
        ? pointer->load (std::memory_order_relaxed) : bounds.defaultValue;
    return std::isfinite (value)
        ? std::clamp (value, bounds.minimum, bounds.maximum) : bounds.defaultValue;
}

taikor::EngineParameters TaikorAudioProcessor::snapshotEngineParameters() const noexcept
{
    const auto read = [this] (int slot) noexcept
    {
        return readParameter (slot);
    };

    taikor::EngineParameters next;
    next.headDiameter = read (slotHeadDiameter) * 0.01f; // centimetres to metres
    next.bodyDepth = read (slotBodyDepth);
    next.tension = read (slotTension);
    next.headMaterial = read (slotHeadMaterial);
    next.shellMaterial = read (slotShellMaterial);
    next.resonantTension = read (slotResonantTension);
    next.cavityCoupling = read (slotCavityCoupling);
    next.headDamping = read (slotHeadDamping);
    next.shellResonance = read (slotShellResonance);
    next.pitch = read (slotPitch);
    next.bachiHardness = read (slotBachiHardness);
    next.strikePosition = read (slotStrikePosition);
    next.velocityDepth = read (slotVelocityDepth);
    next.tensionModulation = read (slotTensionModulation);
    next.strikeNoise = read (slotStrikeNoise);
    next.humanise = read (slotHumanise);
    next.octaveBody = read (slotOctaveBody);
    next.micDistance = juce::jmap (read (slotMicDistance),
                                   ids::minimumMicDistanceCentimetres,
                                   ids::maximumMicDistanceCentimetres, 0.0f, 1.0f);
    next.micSpread = read (slotMicSpread);
    next.stereoWidth = read (slotStereoWidth);
    next.drive = read (slotDrive);
    next.outputGain = juce::Decibels::decibelsToGain (read (slotOutput));
    next.strikeAzimuth = juce::degreesToRadians (read (slotStrikeAzimuth));
    next.performer = juce::roundToInt (read (slotPerformer));
    next.velocityCurve = read (slotVelocityCurve);
    next.outputHighPassHz = read (slotOutputHighPass);
    next.ensembleSize = juce::roundToInt (read (slotEnsembleSize));
    next.ensembleVariation = read (slotEnsembleVariation);
    return next;
}

taikor::TaikoEngine::DrumMeasurements TaikorAudioProcessor::measureDrum (
    int octaveOffset) const noexcept
{
    // The host's rate, because the pitch this reports has to be a partial the
    // engine will actually instantiate at it - see TaikoEngine::soundingMode.
    // Zero means the host has not prepared us yet, and then the measurement's
    // own default stands.
    auto current = snapshotEngineParameters();
    const int azimuth = strikeAzimuthController.load (std::memory_order_relaxed);
    const int position = strikePositionController.load (std::memory_order_relaxed);
    if (azimuth >= 0)
        current.strikeAzimuth = bipolarControllerValue (azimuth)
                              * juce::MathConstants<float>::pi;
    if (position >= 0)
        current.strikePosition = bipolarControllerValue (position);

    const auto rate = displaySampleRate.load (std::memory_order_relaxed);
    const int octave = std::clamp (octaveOffset, taikor::lowestOctaveOffset,
                                    taikor::highestOctaveOffset);
    auto& cached = drumReadouts[static_cast<std::size_t> (octave - taikor::lowestOctaveOffset)];
    const auto features = taikor::TaikoEngine::realismFeatures();
    if (cached.valid && cached.parameters == current && cached.sampleRate == rate
        && cached.features == features)
        return cached.measurements;
    cached.measurements = rate > 0.0
             ? taikor::TaikoEngine::measure (current, octave, 0.0f, rate)
             : taikor::TaikoEngine::measure (current, octave);
    cached.parameters = current;
    cached.sampleRate = rate;
    cached.features = features;
    cached.valid = true;
    return cached.measurements;
}

void TaikorAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    engineReady.store (false, std::memory_order_release);
    uiQueueGeneration.fetch_add (1, std::memory_order_acq_rel);
    discardUiTriggers();
    const auto rate = std::isfinite (sampleRate) && sampleRate > 0.0
        ? std::clamp (sampleRate, 8000.0, 384000.0) : 48000.0;
    const int blockHint = std::clamp (samplesPerBlock, 1, 65536);
    // Publish the parameters before prepare() resets the smoothers, so the
    // first stroke starts at the restored gain rather than gliding to it.
    updateEngineParameters();
    engine.prepare (rate, blockHint);
    engine.prepareOfflineRendering (isNonRealtime() ? 3 : 0);
    lastReverbRoom = juce::roundToInt (readParameter (slotReverbRoom));
    lastReverbMix = readParameter (slotReverbMix);
    reverb.prepare (rate, blockHint, lastReverbRoom, lastReverbMix);
    setLatencySamples (reverb.getLatency());
    meterRelease = static_cast<float> (std::exp (-1.0 / (0.22 * rate)));
    cpuMeter.reset (rate);
    meterLevels.fill (0.0f);
    for (auto& level : outputLevels)
        level.store (0.0f, std::memory_order_relaxed);
    strikeAzimuthController.store (-1, std::memory_order_relaxed);
    strikePositionController.store (-1, std::memory_order_relaxed);
    displaySampleRate.store (rate, std::memory_order_relaxed);
    activeVoiceCount.store (0, std::memory_order_relaxed);
    engineReady.store (true, std::memory_order_release);
}

void TaikorAudioProcessor::releaseResources()
{
    engineReady.store (false, std::memory_order_release);
    cpuMeter.reset();
    uiQueueGeneration.fetch_add (1, std::memory_order_acq_rel);
    discardUiTriggers();
    engine.releaseOfflineRendering();
    engine.allSoundsOff();
    engine.reset();
    reverb.reset();
    meterLevels.fill (0.0f);
    for (auto& level : outputLevels)
        level.store (0.0f, std::memory_order_relaxed);
    strikeAzimuthController.store (-1, std::memory_order_relaxed);
    strikePositionController.store (-1, std::memory_order_relaxed);
    activeVoiceCount.store (0, std::memory_order_relaxed);
    displaySampleRate.store (0.0, std::memory_order_relaxed);
}

bool TaikorAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    // The instrument is a stereo close pair by construction: the two
    // microphones are different points on the drum, and summing them to mono
    // would throw away the model's own image rather than fold it down.
    return layouts.getMainInputChannelSet().isDisabled()
        && layouts.getMainOutputChannelSet() == juce::AudioChannelSet::stereo();
}

void TaikorAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer,
                                         juce::MidiBuffer& midiMessages)
{
    const auto cpuStart = std::chrono::steady_clock::now();
    juce::ScopedNoDenormals noDenormals;
    buffer.clear();

    if (! engineReady.load (std::memory_order_acquire))
    {
        midiMessages.clear();
        return;
    }

    // This only selects an already-prepared pool; realtime always stays serial.
    engine.setOfflineRendering (isNonRealtime());
    if (panicRequested.exchange (false, std::memory_order_acq_rel))
    {
        engine.allSoundsOff();
        reverb.reset();
        engine.setHandDamping (0.0f);
        engine.setPitchBend (0.0f);
        engine.clearStrikeOverrides();
        engine.setRearHeadStrike (false);
        strikeAzimuthController.store (-1, std::memory_order_relaxed);
        strikePositionController.store (-1, std::memory_order_relaxed);
        meterLevels.fill (0.0f);
        activeVoiceCount.store (0, std::memory_order_relaxed);
    }

    const auto numSamples = buffer.getNumSamples();
    if (numSamples == 0 || buffer.getNumChannels() < 2)
    {
        // buffer.clear() above already made every available output finite.
        // No valid stereo time span exists. Keep UI auditions for a real block.
        midiMessages.clear();
        for (auto& level : outputLevels)
            level.store (0.0f, std::memory_order_relaxed);
        return;
    }

    updateEngineParameters();
    const int room = juce::roundToInt (readParameter (slotReverbRoom));
    const float mix = readParameter (slotReverbMix);
    if (room != lastReverbRoom || mix != lastReverbMix)
    {
        reverb.setParameters (room, mix);
        lastReverbRoom = room;
        lastReverbMix = mix;
    }

    // Editor pad strokes are deliberately quantised to the next block boundary.
    // UI strokes precede host MIDI at sample zero. Host events at the same
    // sample retain MidiBuffer insertion order, including notes and resets.
    dispatchUiTriggers();

    int renderedTo = 0;
    unsigned events = 0;

    for (const auto metadata : midiMessages)
    {
        if (events++ == maximumMidiEventsPerBlock)
        {
            // Fail silent when the work budget is exhausted; a stop message
            // later in the discarded suffix must never leave a ringing voice.
            engine.allSoundsOff();
            reverb.reset();
            break;
        }
        const auto eventSample = metadata.samplePosition;
        // The host owns scheduling. Never pull a future note into this block
        // or turn an invalid negative offset into an early hit.
        if (eventSample < 0 || eventSample >= numSamples)
            continue;

        if (eventSample > renderedTo)
        {
            renderAudio (buffer, renderedTo, eventSample - renderedTo);
            renderedTo = eventSample;
        }

        dispatchMidiData (metadata.data, metadata.numBytes);
    }
    midiMessages.clear(); // This instrument does not produce MIDI.

    if (renderedTo < numSamples)
        renderAudio (buffer, renderedTo, numSamples - renderedTo);

    activeVoiceCount.store (engine.getActiveVoiceCount(), std::memory_order_relaxed);
    for (std::size_t channel = 0; channel < meterLevels.size(); ++channel)
        outputLevels[channel].store (meterLevels[channel], std::memory_order_relaxed);
    cpuMeter.record (std::chrono::duration<double> (
        std::chrono::steady_clock::now() - cpuStart).count(), numSamples);
}

void TaikorAudioProcessor::renderAudio (juce::AudioBuffer<float>& buffer,
                                      int start, int samples) noexcept
{
    auto* left = buffer.getWritePointer (0, start);
    auto* right = buffer.getWritePointer (1, start);
    // Establish silence before rendering: a block that retires its last drum
    // can still contain peaks. Pending ensemble hits also prevent this path.
    const bool silent = reverb.isBypassed() && engine.isOutputFrozen();
    engine.process (left, right, samples);
    reverb.process (left, right, samples);
    // The room's exact dry bypass deliberately leaves its input untouched.
    // Keep the host boundary finite before the silent fast path or metering.
    // Finite samples are unchanged; this adds no second limiter or gain stage.
    for (int sample = 0; sample < samples; ++sample)
    {
        if (! std::isfinite (left[sample]))
            left[sample] = 0.0f;
        if (! std::isfinite (right[sample]))
            right[sample] = 0.0f;
    }
    if (silent)
    {
        if (meterLevels[0] != 0.0f || meterLevels[1] != 0.0f)
            for (int sample = 0; sample < samples; ++sample)
            {
                meterLevels[0] *= meterRelease;
                meterLevels[1] *= meterRelease;
            }
        return;
    }
    // Meter the finished wet/dry output, including a room ringing after drums
    // retire. Keeping processing inside MIDI slices makes panic sample-accurate.
    for (int sample = 0; sample < samples; ++sample)
    {
        meterLevels[0] = std::max (std::abs (left[sample]), meterLevels[0] * meterRelease);
        meterLevels[1] = std::max (std::abs (right[sample]), meterLevels[1] * meterRelease);
    }
}

void TaikorAudioProcessor::dispatchMidiData (const juce::uint8* data,
                                             int numBytes) noexcept
{
    // Only complete supported channel messages are meaningful to this omni
    // one-shot instrument. Channels share the same live gestures; note-offs
    // and zero-velocity note-ons do not truncate a struck drum's natural tail.
    if (data == nullptr || numBytes != 3 || data[1] >= 0x80u || data[2] >= 0x80u)
        return;

    const auto status = static_cast<unsigned> (data[0]);
    const auto kind = status & 0xf0u;

    if (kind == 0x90u && data[2] != 0)
    {
        const auto midiNote = static_cast<int> (data[1]);
        const auto velocity = static_cast<float> (data[2]) / 127.0f;

        if (engine.triggerMidi (midiNote, velocity))
            if (const auto articulation = taikor::articulationForMidiNote (midiNote))
                registerTrigger (*articulation);
    }
    else if (kind == 0xe0u)
    {
        // Pressing the head raises its tension, so the wheel bends the drum up.
        const auto raw = static_cast<int> (data[1] & 0x7fu)
                       | (static_cast<int> (data[2] & 0x7fu) << 7);
        engine.setPitchBend (static_cast<float> (raw - 8192) / 8192.0f);
    }
    else if (kind == 0xb0u)
    {
        const auto controller = data[1] & 0x7fu;
        const auto rawValue = static_cast<int> (data[2] & 0x7fu);
        const auto value = static_cast<float> (rawValue) / 127.0f;
        // Standard bipolar MIDI controls have a real centre at 64. Mapping the
        // whole byte through raw/127 would put that value slightly positive,
        // so use the 64 steps below centre and 63 above it explicitly.
        const float bipolarValue = bipolarControllerValue (rawValue);

        if (controller == 1u)
        {
            // A hand laid on the head. It damps what is still ringing without
            // changing anything about the strokes that follow.
            engine.setHandDamping (value);
        }
        else if (controller == 16u)
        {
            strikeAzimuthController.store (rawValue, std::memory_order_relaxed);
            engine.setStrikeAzimuthOverride (
                bipolarValue * juce::MathConstants<float>::pi);
        }
        else if (controller == 17u)
        {
            strikePositionController.store (rawValue, std::memory_order_relaxed);
            engine.setStrikePositionOverride (bipolarValue);
        }
        else if (controller == 18u)
        {
            engine.setRearHeadStrike (rawValue >= 64);
        }
        else if (controller == 121u)
        {
            // Reset All Controllers returns every live gesture to its host
            // parameter: hand off, wheel centred, and strike overrides cleared.
            engine.setHandDamping (0.0f);
            engine.setPitchBend (0.0f);
            engine.clearStrikeOverrides();
            engine.setRearHeadStrike (false);
            strikeAzimuthController.store (-1, std::memory_order_relaxed);
            strikePositionController.store (-1, std::memory_order_relaxed);
        }
        else if (controller == 120u || controller == 123u)
        {
            engine.allSoundsOff();
            reverb.reset();
        }
    }
}

void TaikorAudioProcessor::triggerFromUi (taikor::Articulation articulation,
                                          int octaveOffset, float velocity) noexcept
{
    // Capture the generation before checking readiness so prepare/release or
    // panic cannot relabel an in-flight old audition as a new one.
    const auto generation = uiQueueGeneration.load (std::memory_order_acquire);
    if (! engineReady.load (std::memory_order_acquire))
        return;
    enqueueUiTrigger (articulation, octaveOffset, velocity, generation);
}

void TaikorAudioProcessor::requestPanic() noexcept
{
    // Events published before this generation change are stale even if their
    // producer races the audio-thread queue flush.
    uiQueueGeneration.fetch_add (1, std::memory_order_acq_rel);
    strikeAzimuthController.store (-1, std::memory_order_relaxed);
    strikePositionController.store (-1, std::memory_order_relaxed);
    panicRequested.store (true, std::memory_order_release);
}

void TaikorAudioProcessor::enqueueUiTrigger (taikor::Articulation articulation,
                                             int octaveOffset, float velocity,
                                             std::uint32_t generation) noexcept
{
    if (! isValidArticulation (articulation) || ! std::isfinite (velocity)
        || velocity <= 0.0f)
        return;

    if (uiProducerBusy.test_and_set (std::memory_order_acquire))
        return;
    const auto write = uiWriteIndex.load (std::memory_order_relaxed);
    const auto next = (write + 1u) % uiQueueCapacity;

    if (next == uiReadIndex.load (std::memory_order_acquire))
    {
        // A one-shot audition can be dropped; the message thread must not wait.
        uiProducerBusy.clear (std::memory_order_release);
        return;
    }

    uiTriggerQueue[write] = {
        articulation,
        juce::jlimit (taikor::lowestOctaveOffset, taikor::highestOctaveOffset,
                      octaveOffset),
        juce::jlimit (0.0f, 1.0f, velocity),
        generation
    };
    uiWriteIndex.store (next, std::memory_order_release);
    uiProducerBusy.clear (std::memory_order_release);
}

void TaikorAudioProcessor::dispatchUiTriggers() noexcept
{
    auto read = uiReadIndex.load (std::memory_order_relaxed);
    const auto write = uiWriteIndex.load (std::memory_order_acquire);

    while (read != write)
    {
        const auto event = uiTriggerQueue[read];
        if (event.generation == uiQueueGeneration.load (std::memory_order_acquire))
        {
            engine.trigger (event.articulation, event.octaveOffset, event.velocity);
            registerTrigger (event.articulation);
        }
        read = (read + 1u) % uiQueueCapacity;
    }

    uiReadIndex.store (read, std::memory_order_release);
}

void TaikorAudioProcessor::discardUiTriggers() noexcept
{
    uiReadIndex.store (uiWriteIndex.load (std::memory_order_acquire),
                       std::memory_order_release);
}

void TaikorAudioProcessor::registerTrigger (taikor::Articulation articulation) noexcept
{
    if (isValidArticulation (articulation))
        triggerCounters[static_cast<std::size_t> (articulation)]
            .fetch_add (1, std::memory_order_release);
}

std::uint32_t TaikorAudioProcessor::getTriggerCounter (
    taikor::Articulation articulation) const noexcept
{
    return isValidArticulation (articulation)
        ? triggerCounters[static_cast<std::size_t> (articulation)]
              .load (std::memory_order_acquire)
        : 0u;
}

void TaikorAudioProcessor::updateEngineParameters() noexcept
{
    engine.setParameters (snapshotEngineParameters());
}

void TaikorAudioProcessor::getStateInformation (juce::MemoryBlock& destinationData)
{
    if (const auto xml = parameters.copyState().createXml())
        copyXmlToBinary (*xml, destinationData);
}

void TaikorAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    if (data == nullptr || sizeInBytes <= 8 || sizeInBytes > maximumStateBytes)
        return;
    const auto* bytes = static_cast<const char*> (data);
    const auto length = juce::ByteOrder::littleEndianInt (bytes + 4);
    if (length == 0 || length > static_cast<juce::uint32> (sizeInBytes - 8)
        || ! juce::CharPointer_UTF8::isValidString (bytes + 8, static_cast<int> (length))
        || ! hasSafeXmlStructure ({ bytes + 8, static_cast<std::size_t> (length) }))
        return;
    const auto xml = getXmlFromBinary (data, sizeInBytes);
    if (xml != nullptr && xml->hasTagName ("TAIKOR_STATE"))
    {
        juce::ValueTree restoredState { "TAIKOR_STATE" };
        for (const auto* id : parameterIds)
        {
            const auto* parameter = parameters.getParameter (id);
            const auto& range = parameter->getNormalisableRange();
            float value = parameter->convertFrom0to1 (parameter->getDefaultValue());
            // First matching entry wins, including when its value is invalid.
            // Unknown nodes do not enter the canonical parameter state.
            for (const auto* child : xml->getChildIterator())
            {
                if (! child->hasTagName ("PARAM") || child->getStringAttribute ("id") != id)
                    continue;
                double stored = 0.0;
                if (child->getNumChildElements() == 0
                    && readFiniteDecimal (child->getStringAttribute ("value"), stored))
                    value = range.snapToLegalValue (static_cast<float> (
                        std::clamp (stored, static_cast<double> (range.start),
                                    static_cast<double> (range.end))));
                break;
            }
            juce::ValueTree child { "PARAM" };
            child.setProperty ("id", id, nullptr);
            child.setProperty ("value", value, nullptr);
            restoredState.appendChild (child, nullptr);
        }
        parameters.replaceState (restoredState);
        requestPanic();
    }
}

juce::AudioProcessorEditor* TaikorAudioProcessor::createEditor()
{
    return new TaikorAudioProcessorEditor (*this);
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new TaikorAudioProcessor();
}
