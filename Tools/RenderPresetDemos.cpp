// Taikor's description for the shared preset demo renderer (README G3): the
// factory bank in program order, played through the same engine and MIDI
// mapping as the plug-in. Presets are dry, so the plug-in's convolution room
// is not part of their signal path.
#include "PresetDemoRenderer.h"

#include "DSP/FactoryPresets.h"
#include "DSP/MidiInput.h"
#include "DSP/OutputLimiter.h"

#include <cmath>
#include <memory>

namespace
{
class TaikorVoice final : public protocodus::demos::Voice
{
public:
    TaikorVoice (std::size_t preset, double sampleRate, int maxBlockFrames)
    {
        engine->prepare (sampleRate, maxBlockFrames);
        engine->setParameters (taikor::factoryPresets()[preset].parameters);
    }
    void Midi (const std::uint8_t* bytes, int size) override { (void) taikor::applyMidi (*engine, bytes, size); }
    void Render (float* left, float* right, int frames) override { engine->process (left, right, frames); }

private:
    // The engine is large; keep it off the worker thread's stack.
    std::unique_ptr<taikor::EnsembleEngine> engine = std::make_unique<taikor::EnsembleEngine>();
};
} // namespace

int main (int argc, char** argv)
{
    protocodus::demos::Instrument instrument;
    instrument.name = "Taikor";
    for (const auto& preset : taikor::factoryPresets())
        instrument.presets.emplace_back (preset.name);
    // The output limiter holds every sample at or below -1 dBFS; anything
    // above it means the limiter was bypassed.
    instrument.ceiling = std::nextafter (taikor::OutputLimiter::ceiling, 1.0f);
    // The engine's ensemble and noise variation is hashed from fixed internal
    // state, so the renderer's seed has nothing to set.
    instrument.open = [] (std::size_t preset, double sampleRate, int maxBlockFrames, std::uint32_t)
    {
        return std::make_unique<TaikorVoice> (preset, sampleRate, maxBlockFrames);
    };
    return protocodus::demos::Main (argc, argv, instrument);
}
