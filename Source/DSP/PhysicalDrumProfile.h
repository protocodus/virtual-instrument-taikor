#pragma once

#include <array>
#include <string_view>

namespace taikor::physical
{
// Specimen priors, not measurements of a particular maker's instrument. Keep
// thickness, areal mass and loss explicit: thin natural hide is not polyester.
struct Skin
{
    float arealDensity { 1.0f };       // kg/m^2
    float thickness { 0.001f };       // metres
    float youngsModulus { 3.5e9f };    // Pa
    float poissonRatio { 0.30f };
    float lossAngle { 0.0005f };
    float viscousLossSeconds { 1.5e-7f };
};

struct HeadPair
{
    Skin batter;
    Skin rear;
    float rearTensionRatio { 1.0f };
};

struct Profile
{
    std::string_view name;
    HeadPair heads;
    // A small, fixed boundary/material perturbation. These are design priors,
    // not measured tolerances of rope or tack spacing.
    float pairSplit { 0.002f };
    float axisRadians { 0.0f };
    // Independent support-loss prior. This scales the existing smooth
    // low-frequency mounting loss, not a selected acoustic mode or EQ band.
    float supportLossScale { 1.0f };
};

[[nodiscard]] const Profile& profileForFamily (int family) noexcept;

// HEAD retains its 0..1 host range. 0.75 is each family's natural-hide specimen;
// above 0.20 it trims thickness without changing the specimen into synthetic
// film. Only the lowest fifth transitions to the explicit thin-film endpoint.
// Family mix zero always means the reference o-daiko, including its rear head.
[[nodiscard]] HeadPair resolveHeads (int family, float familyMix,
                                      float headMaterialControl) noexcept;
[[nodiscard]] float flexuralRigidity (const Skin&) noexcept;

struct Tuning
{
    float batterTensionScale { 1.0f };
    float rearTensionScale { 1.0f };
};
// Explicit model voicing at factory controls, not measured traditional tuning.
// Blends positive tensions geometrically; air/microphone settings never feed
// back into the instrument's dimensions or these construction priors.
[[nodiscard]] Tuning tuningForPad (int pad, float familyMix) noexcept;

struct AngularBasis
{
    int circumferentialOrder { 0 };
    std::array<float, 2> frequencyScale { 1.0f, 1.0f };
    // Rows map [cos(m theta), sin(m theta)] to the two principal shapes.
    std::array<std::array<float, 2>, 2> rotation {{ { 1.0f, 0.0f }, { 0.0f, 1.0f } }};

    [[nodiscard]] float project (float angleRadians, int branch) const noexcept;
};

// Use the SAME returned basis for frequency, strike, microphone and contact
// projections. Changing only the frequencies would rotate no nodal diameter.
// No RNG, history, velocity or Humanise input: the imperfection belongs to the
// physical head and cannot move its authored strike coordinates.
[[nodiscard]] AngularBasis basisForMode (int family, float familyMix,
                                         int circumferentialOrder,
                                         int radialOrder) noexcept;
} // namespace taikor::physical
