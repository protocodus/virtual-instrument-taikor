#pragma once

#include "TaikoEngine.h"

#include <array>
#include <cstddef>

namespace taikor
{
// The factory bank: the plug-in's program list and the source of one README
// demo per preset. A preset is a complete set of engine controls; every field
// not named in FactoryPresets.cpp keeps its EngineParameters default. Presets
// are dry: the plug-in's convolution room is switched off when one is loaded.
struct FactoryPreset
{
    const char* name;
    // Stable file stem for the preset's MIDI phrase and demo,
    // Docs/audio/midi/<NN>-<slug>.mid. Never rename a shipped slug.
    const char* slug;
    EngineParameters parameters;
};

inline constexpr std::size_t factoryPresetCount = 8;

[[nodiscard]] const std::array<FactoryPreset, factoryPresetCount>& factoryPresets() noexcept;
} // namespace taikor
