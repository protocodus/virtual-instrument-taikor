// The standard's U4 parameter-extremes matrix over the plug-in processor, its
// parameter list, its factory programs and its state chunk.
//
//   TaikorParameterMatrixTests [--rotate] [--jobs N]
//
// --rotate runs each case at one rate and block size in turn instead of all
// nine, for the sanitizer pass on every push; the main build runs the full
// matrix.
#include "PluginProcessor.h"

#include "ParameterMatrix.h"

#include <cstring>
#include <memory>
#include <string>

namespace
{
class TaikorDevice final : public protocodus::matrix::JuceDevice
{
public:
    TaikorDevice() : JuceDevice (std::make_unique<TaikorAudioProcessor>()) {}

    protocodus::matrix::Limits LimitsFor (const protocodus::matrix::Case&) const override
    {
        protocodus::matrix::Limits limits;
        // A drum ignores note length: the matrix's chord is one stroke at
        // note-on (60 is the nagado-daiko's Don; 64 and 67 are not mapped), so
        // nothing is held, and a heavily damped head is rightly silent by the
        // time the chord's note-off arrives.
        limits.soundsWhileHeld = false;
        // Without these two limits 471 of 615 rotation runs failed on silent
        // held notes, and two all-pairs rows (moderate damping, no room) still
        // rang at -45 and -43 dBFS two seconds after the stroke. A stroke's
        // tail decays freely, so silence is due later.
        limits.releaseSeconds = 6.0;
        return limits;
    }
};
} // namespace

int main (int argc, char** argv)
{
    // JUCE's message manager must exist before a processor is created.
    juce::ScopedJuceInitialiser_GUI juce;
    protocodus::matrix::Options options;
    for (int i = 1; i < argc; ++i)
    {
        if (std::strcmp (argv[i], "--rotate") == 0)
            options.everyConfiguration = false;
        else if (std::strcmp (argv[i], "--jobs") == 0 && i + 1 < argc)
            options.jobs = std::stoi (argv[++i]);
    }
    return protocodus::matrix::Run ([] { return std::make_unique<TaikorDevice>(); }, options) == 0 ? 0 : 1;
}
