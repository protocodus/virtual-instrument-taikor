#include "MidiInput.h"

namespace taikor
{
MidiAction applyMidi (EnsembleEngine& engine, const std::uint8_t* data, int size) noexcept
{
    if (data == nullptr || size < 3)
        return MidiAction::none;

    const auto kind = data[0] & 0xf0u;
    const auto first = static_cast<int> (data[1] & 0x7fu);
    const auto second = static_cast<int> (data[2] & 0x7fu);

    if (kind == 0x90u && second != 0)
        return engine.triggerMidi (first, static_cast<float> (second) / 127.0f)
            ? MidiAction::struck : MidiAction::none;

    if (kind == 0xe0u)
    {
        // Pressing the head raises its tension, so the wheel bends the drum up.
        engine.setPitchBend (static_cast<float> ((first | (second << 7)) - 8192) / 8192.0f);
        return MidiAction::none;
    }

    if (kind != 0xb0u)
        return MidiAction::none;

    switch (first)
    {
        case 1:
            // A hand laid on the head. It damps what is still ringing without
            // changing anything about the strokes that follow.
            engine.setHandDamping (static_cast<float> (second) / 127.0f);
            return MidiAction::none;
        case 16:
            engine.setStrikeAzimuthOverride (bipolarControllerValue (second) * 3.14159265358979323846f);
            return MidiAction::strikeAzimuth;
        case 17:
            engine.setStrikePositionOverride (bipolarControllerValue (second));
            return MidiAction::strikePosition;
        case 18:
            engine.setRearHeadStrike (second >= 64);
            return MidiAction::none;
        case 121:
            // Reset All Controllers returns every live gesture to its host
            // parameter: hand off, wheel centred, and strike overrides cleared.
            engine.setHandDamping (0.0f);
            engine.setPitchBend (0.0f);
            engine.clearStrikeOverrides();
            engine.setRearHeadStrike (false);
            return MidiAction::resetControllers;
        case 120:
        case 123:
            engine.allSoundsOff();
            return MidiAction::allSoundsOff;
        default:
            return MidiAction::none;
    }
}
} // namespace taikor
