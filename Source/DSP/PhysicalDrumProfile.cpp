#include "DSP/PhysicalDrumProfile.h"

#include <algorithm>
#include <cmath>

namespace taikor::physical
{
namespace
{
// Construction evidence: natural cowhide on tacked nagado; horse/cowhide on
// lightweight rope-laced okedo; thick tensioned cowhide is a valid tsuke-shime
// design, rather than the previous thin-film-like shime endpoint.
// https://www.miyamoto-unosuke.co.jp/products/13030402014004
// https://www.asano.jp/files/kanade/ (different head thicknesses, independent ropes)
// https://www.asano.jp/product/145
// No source above supplies these constitutive constants. They are explicit,
// independently calibratable priors; 3.5 GPa is the existing treated-hide
// modulus assumption, not a measured modulus for horsehide or every shime.
constexpr std::array<Profile, 4> profiles {{
    { "Tacked heavy natural hide",
      { { 1.10f, 0.00110f, 3.50e9f, 0.30f, 5.2e-4f, 1.85e-7f },
        { 1.02f, 0.00102f, 3.35e9f, 0.30f, 5.6e-4f, 1.95e-7f }, 0.98f },
      0.0020f, 0.13f },
    { "Tacked medium natural hide",
      { { 0.94f, 0.00094f, 3.50e9f, 0.30f, 4.8e-4f, 1.60e-7f },
        { 0.89f, 0.00089f, 3.40e9f, 0.30f, 5.1e-4f, 1.75e-7f }, 1.02f },
      // Firmer supported-barrel mounting keeps the initial coupled boom
      // shorter than the freely ringing body. This is model voicing, not a
      // measured universal damping factor for every nagado stand.
      0.0028f, 0.39f, 1.75f },
    { "Rope-laced thin natural hide",
      { { 0.62f, 0.00062f, 3.20e9f, 0.30f, 3.8e-4f, 1.10e-7f },
        { 0.73f, 0.00073f, 3.10e9f, 0.30f, 4.5e-4f, 1.35e-7f }, 0.93f },
      0.0060f, 0.71f },
    { "Tensioned thick natural hide",
      { { 1.20f, 0.00120f, 3.50e9f, 0.30f, 5.7e-4f, 2.05e-7f },
        { 1.12f, 0.00112f, 3.40e9f, 0.30f, 6.0e-4f, 2.20e-7f }, 1.025f },
      0.0042f, 1.02f },
}};

float unit (float value, float fallback = 0.0f) noexcept
{
    return std::isfinite (value) ? std::clamp (value, 0.0f, 1.0f) : fallback;
}

float mix (float a, float b, float proportion) noexcept
{
    return a + (b - a) * proportion;
}

Skin mixSkin (const Skin& a, const Skin& b, float proportion) noexcept
{
    return { mix (a.arealDensity, b.arealDensity, proportion),
             mix (a.thickness, b.thickness, proportion),
             mix (a.youngsModulus, b.youngsModulus, proportion),
             mix (a.poissonRatio, b.poissonRatio, proportion),
             mix (a.lossAngle, b.lossAngle, proportion),
             mix (a.viscousLossSeconds, b.viscousLossSeconds, proportion) };
}

Skin trimSkin (Skin skin, float control) noexcept
{
    constexpr float identityBoundary = 0.20f;
    const float naturalControl = std::max (control, identityBoundary);
    const float thicknessScale = std::exp2 (1.35f * (naturalControl - 0.75f));
    skin.arealDensity *= thicknessScale;
    skin.thickness *= thicknessScale;
    // Loss is a material property with its own weak trim, not the cube of h.
    const float lossScale = std::exp2 (0.40f * (naturalControl - 0.75f));
    skin.lossAngle *= lossScale;
    skin.viscousLossSeconds *= lossScale;
    if (control >= identityBoundary)
        return skin;

    constexpr Skin film { 0.30f, 0.30f / 1390.0f, 4.0e9f, 0.30f,
                          1.6e-4f, 8.25e-8f };
    const float fraction = control / identityBoundary;
    const float smoothFraction = fraction * fraction * (3.0f - 2.0f * fraction);
    return mixSkin (film, skin, smoothFraction);
}
} // namespace

const Profile& profileForFamily (int family) noexcept
{
    return profiles[static_cast<std::size_t> (std::clamp (family, 0, 3))];
}

HeadPair resolveHeads (int family, float familyMix, float headMaterialControl) noexcept
{
    const auto& reference = profiles.front().heads;
    const auto& selected = profileForFamily (family).heads;
    const float blend = unit (familyMix);
    const float control = unit (headMaterialControl, 0.75f);
    return { trimSkin (mixSkin (reference.batter, selected.batter, blend), control),
             trimSkin (mixSkin (reference.rear, selected.rear, blend), control),
             mix (reference.rearTensionRatio, selected.rearTensionRatio, blend) };
}

float flexuralRigidity (const Skin& skin) noexcept
{
    return skin.youngsModulus * skin.thickness * skin.thickness * skin.thickness
         / (12.0f * (1.0f - skin.poissonRatio * skin.poissonRatio));
}

Tuning tuningForPad (int pad, float familyMix) noexcept
{
    // Factory neutral-contact voicing for the stated dimensions and skins.
    // The one-drum mode is an intentionally extended retuning, not an assertion
    // that a real tacked hide could survive these tensions. The first retuning
    // loosens the rear relative to the batter to avoid a dominant-mode gap.
    constexpr std::array<float, 4> oneDrum { 1.0f, 14.5514f, 73.5787f, 226.375f };
    constexpr std::array<float, 4> family { 1.0f, 1.05350f, 1.77301f, 2.54110f };
    constexpr std::array<float, 4> rear { 1.0f, 0.70f, 1.0f, 1.0f };
    const auto index = static_cast<std::size_t> (std::clamp (pad, 0, 3));
    const float blend = unit (familyMix);
    return { std::exp (mix (std::log (oneDrum[index]), std::log (family[index]), blend)),
             std::exp (mix (std::log (rear[index]), 0.0f, blend)) };
}

float AngularBasis::project (float angleRadians, int branch) const noexcept
{
    if (circumferentialOrder == 0)
        return branch == 0 ? 1.0f : 0.0f;
    const float angle = static_cast<float> (circumferentialOrder) * angleRadians;
    const auto& row = rotation[static_cast<std::size_t> (std::clamp (branch, 0, 1))];
    return row[0] * std::cos (angle) + row[1] * std::sin (angle);
}

AngularBasis basisForMode (int family, float familyMix, int circumferentialOrder,
                           int radialOrder) noexcept
{
    AngularBasis result;
    result.circumferentialOrder = std::max (circumferentialOrder, 0);
    if (result.circumferentialOrder == 0)
        return result;

    const float blend = unit (familyMix);
    const auto& reference = profiles.front();
    const auto& selected = profileForFamily (family);
    const float order = static_cast<float> (result.circumferentialOrder);
    const float radial = static_cast<float> (std::max (radialOrder, 1));
    const float split = mix (reference.pairSplit, selected.pairSplit, blend)
                      / (std::sqrt (radial) * (1.0f + 0.18f * (order - 1.0f)));
    // The perturbation is trace-free in squared frequency, so splitting does
    // not also impose an arbitrary common-mode tuning shift.
    result.frequencyScale = { std::sqrt (1.0f + 2.0f * split),
                              std::sqrt (1.0f - 2.0f * split) };
    const float axis = mix (reference.axisRadians, selected.axisRadians, blend);
    const float rotation = order * axis + 0.11f * (radial - 1.0f);
    const float c = std::cos (rotation), s = std::sin (rotation);
    result.rotation = {{ { c, s }, { -s, c } }};
    return result;
}
} // namespace taikor::physical
