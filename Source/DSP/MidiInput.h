#pragma once

#include "EnsembleEngine.h"

#include <cstdint>

namespace taikor
{
// What one MIDI channel message did to the engine, so a front end can do its
// own bookkeeping (UI flashes, controller readouts, its reverb) on top.
enum class MidiAction
{
    none,
    struck,           // a note-on that triggered a drum
    strikeAzimuth,    // CC16
    strikePosition,   // CC17
    resetControllers, // CC121
    allSoundsOff      // CC120 or CC123
};

// MIDI's conventional centre is the exact value 64: 64 steps below it, 63
// above it. Mapping the whole byte through raw/127 would put the centre
// slightly positive.
[[nodiscard]] constexpr float bipolarControllerValue (int rawValue) noexcept
{
    return rawValue < 64
        ? static_cast<float> (rawValue - 64) / 64.0f
        : static_cast<float> (rawValue - 64) / 63.0f;
}

// Applies one channel message the way every front end plays the engine: the
// plug-in, the demo renderer and a future Rack adapter share this mapping.
MidiAction applyMidi (EnsembleEngine&, const std::uint8_t* data, int size) noexcept;
} // namespace taikor
