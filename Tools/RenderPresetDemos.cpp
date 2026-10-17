// Taikor's description for the shared preset demo renderer (README G3): the
// factory bank in program order, played through the same engine and MIDI
// mapping as the plug-in. Presets are dry, so the plug-in's convolution room
// is not part of their signal path.
//
//   RenderTaikorPresetDemos [--more] <shared renderer arguments>
//
// --more renders the README's More demos (G3d) instead: playing techniques
// the presets do not feature, each from Docs/audio/midi/more/<NN>-<slug>.mid
// through the factory preset named in the table below.
#include "PresetDemoRenderer.h"

#include "DSP/FactoryPresets.h"
#include "DSP/MidiInput.h"
#include "DSP/OutputLimiter.h"

#include <cmath>
#include <cstring>
#include <memory>
#include <vector>

namespace
{
struct Showcase
{
    const char* name;
    std::size_t preset;
};

// The order is the file numbering; renaming one renames its files.
constexpr Showcase showcases[] {
    { "Rolls and Flams", 2 },  // Close Studio: the pair hears each grace note
    { "Dynamics", 0 },
    { "Hand Damping", 0 },
    { "Pitch Wheel", 0 },
    { "Strike Placement", 2 },
};

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
    const bool more = argc > 1 && std::strcmp (argv[1], "--more") == 0;
    if (more)
    {
        argv[1] = argv[0];
        ++argv;
        --argc;
    }

    protocodus::demos::Instrument instrument;
    instrument.name = "Taikor";
    std::vector<std::size_t> presetFor;
    if (more)
        for (const auto& showcase : showcases)
        {
            instrument.presets.emplace_back (showcase.name);
            presetFor.push_back (showcase.preset);
        }
    else
        for (std::size_t index = 0; index < taikor::factoryPresets().size(); ++index)
        {
            instrument.presets.emplace_back (taikor::factoryPresets()[index].name);
            presetFor.push_back (index);
        }
    // The output limiter holds every sample at or below -1 dBFS; anything
    // above it means the limiter was bypassed.
    instrument.ceiling = std::nextafter (taikor::OutputLimiter::ceiling, 1.0f);
    // The engine's ensemble and noise variation is hashed from fixed internal
    // state, so the renderer's seed has nothing to set.
    instrument.open = [presetFor] (std::size_t demo, double sampleRate, int maxBlockFrames, std::uint32_t)
    {
        return std::make_unique<TaikorVoice> (presetFor[demo], sampleRate, maxBlockFrames);
    };
    return protocodus::demos::Main (argc, argv, instrument);
}
