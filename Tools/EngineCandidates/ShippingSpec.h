#pragma once

// Everything the candidate engines take from the shipping engine: the drum,
// the stroke and the microphones it resolves for one set of controls, a
// render of the shipping engine itself, and its output stage as a filter that
// any raw microphone pair can be passed through.

#include "CandidateEngine.h"
#include "DSP/TaikoEngine.h"

#include <memory>
#include <string>
#include <vector>

namespace taikor::candidates
{

struct Take
{
    std::vector<float> left, right;
};

struct ShippingStroke
{
    EngineParameters parameters;   // the controls every letter shares
    int midiNote = 60;
    float velocity = 1.0f;         // 0..1, MIDI velocity / 127
    double sampleRate = 48000.0;
    int blockSize = 256;
    int strikeSample = 4800;       // the note-on lands on this sample
    int totalSamples = 216000;
};

// The factory controls of the plug-in (Mic Distance 16 cm, Output -22.5 dB)
// with Humanise 0, so the authored strike position is exact.
[[nodiscard]] EngineParameters listeningTestParameters();

struct ResolvedShipping
{
    DrumSpec drum;
    StrokeSpec stroke;
    MicSpec mics;
    std::string report;            // the resolved numbers, for key.md
};

// Reads the shipping engine's resolved drum, its first stroke after reset
// (including the seeded +/-0.5 % speed and +/-0.4 % contact-time variation
// that remains at Humanise 0) and its capsule positions.
[[nodiscard]] ResolvedShipping resolveShipping (const ShippingStroke&);

// Fits sigma0 + sigma1 xi(omega) to the head's shipping loss law in the
// relative least-squares sense over 60 Hz - 6 kHz (log-spaced), with both
// parameters kept non-negative. Writes head.sigma0 and head.sigma1.
void fitTwoParameterLoss (HeadSpec& head);

// Squared wavenumber of the stiff membrane at omega:
// T k^2 + D k^4 = sigma omega^2.
[[nodiscard]] double stiffMembraneXi (const HeadSpec& head, double omega) noexcept;

// The shipping engine rendered exactly as the plug-in renders it at
// Ensemble Size 1: EnsembleEngine::process in blockSize blocks.
[[nodiscard]] Take renderShipping (const ShippingStroke&);
// The same stroke's raw microphone pair (TaikoEngine::processRaw), which is
// what the output stage receives.
[[nodiscard]] Take renderShippingRaw (const ShippingStroke&);

// The shipping output stage (width, 12 Hz DC blocker, Drive, Low Cut, Output
// gain and the -1 dBFS limiter) applied to an external raw pair, by running an
// idle TaikoEngine's processWithEnsemble in blockSize blocks. Passing the
// shipping engine's own raw pair reproduces renderShipping bit for bit.
[[nodiscard]] Take applyShippingOutputStage (const ShippingStroke&, const Take& raw);

} // namespace taikor::candidates
