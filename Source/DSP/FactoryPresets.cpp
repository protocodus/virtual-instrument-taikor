#include "FactoryPresets.h"

namespace taikor
{
namespace
{
// Mic Distance is normalised over the engine's 3-40 cm range.
constexpr float micDistanceFor (float centimetres) noexcept
{
    return (centimetres - 3.0f) / (40.0f - 3.0f);
}

template <typename Edit>
FactoryPreset preset (const char* name, const char* slug, Edit edit) noexcept
{
    FactoryPreset result { name, slug, {} };
    edit (result.parameters);
    return result;
}
} // namespace

// Each preset is a placement or a playing context for the same modelled
// family, not a different instrument: they move the controls a player of a
// real kumi-daiko group would change (who plays, how hard the bachi are, how
// tight the heads are pulled, where the microphones stand).
const std::array<FactoryPreset, factoryPresetCount>& factoryPresets() noexcept
{
    static const std::array<FactoryPreset, factoryPresetCount> bank {{
        // The parameter defaults: the four drums of the family as documented.
        preset ("Kumi-daiko Stage", "kumi-daiko-stage", [] (EngineParameters&) {}),
        // Four players on each drum, as a festival group plays in unison.
        preset ("Festival Ensemble", "festival-ensemble", [] (EngineParameters& p)
        {
            p.ensembleSize = 4;
            p.ensembleVariation = 0.5f;
            p.humanise = 0.55f;
            p.micDistance = micDistanceFor (25.0f);
            p.micSpread = 0.7f;
        }),
        // The pair just off the head, a tight player and more stick on hide.
        preset ("Close Studio", "close-studio", [] (EngineParameters& p)
        {
            p.micDistance = micDistanceFor (6.0f);
            p.micSpread = 0.4f;
            p.humanise = 0.25f;
            p.strikeNoise = 0.45f;
        }),
        // Heads pulled up and partly damped, hard bachi: the dry, articulate
        // sound of a ji-uchi base pattern that has to cut through.
        preset ("Tight Ji-uchi", "tight-ji-uchi", [] (EngineParameters& p)
        {
            p.tension = 0.78f;
            p.headDamping = 0.68f;
            p.bachiHardness = 0.85f;
            p.tensionModulation = 0.3f;
        }),
        // Soft wood and a gentle controller response for quiet passages.
        preset ("Soft Bachi", "soft-bachi", [] (EngineParameters& p)
        {
            p.bachiHardness = 0.25f;
            p.strikeNoise = 0.2f;
            p.velocityCurve = -0.5f;
        }),
        // One o-daiko retuned across the keyboard instead of four drums.
        preset ("Single O-daiko", "single-odaiko", [] (EngineParameters& p)
        {
            p.octaveBody = 0.0f;
            p.tension = 0.55f;
            p.cavityCoupling = 0.9f;
        }),
        // Synthetic film heads on light laminated shells, as many modern
        // groups use for touring and practice drums.
        preset ("Synthetic Heads", "synthetic-heads", [] (EngineParameters& p)
        {
            p.headMaterial = 0.0f;
            p.shellMaterial = 0.3f;
            p.headDamping = 0.4f;
        }),
        // The pair at its furthest and widest: the drums as heard from the
        // front of a small stage.
        preset ("Distant Pair", "distant-pair", [] (EngineParameters& p)
        {
            p.micDistance = micDistanceFor (40.0f);
            p.micSpread = 0.85f;
            p.stereoWidth = 0.6f;
        }),
    }};
    return bank;
}
} // namespace taikor
