#pragma once

#include <array>
#include <cstddef>

namespace taikor::shellboundary
{
inline constexpr std::size_t modeCount = 6;

struct Input
{
    double headRadiusMetres = 0.275;
    double depthMetres = 0.5;
    double shellMaterial = 0.8;
    // O-daiko, nagado, okedo, tsuke-shime, respectively. Family mix zero
    // retains the first construction; one selects the requested construction.
    int family = 0;
    double familyMix = 1.0;
};

struct Modes
{
    std::array<double, modeCount> frequencyHz {};
    std::array<double, modeCount> amplitudeDecay {};
    std::array<double, modeCount> modalMassKg {};
    std::array<double, modeCount> boundaryDampingNsPerMetre {};
    // One cosine-oriented, radially moving ring coordinate of each order.
    // This reduced basis omits axial shell modes and the orthogonal doublet.
    std::array<int, modeCount> circumferentialOrder { 2, 3, 4, 5, 6, 7 };
    double axialYoungsModulusPa = 0.0;
    double transverseYoungsModulusPa = 0.0;
    double woodDensityKgPerM3 = 0.0;
    double wallThicknessMetres = 0.0;
    // Microphone voicing only, not mechanical mass or a measured radiation fit.
    // The thin stave shell needs less gain for realistic rim-accent headroom.
    double observationCalibration = 1100.0;
};

// Orthotropic ring reduction inspired by the significance of cross-grain
// stiffness measured by Hwang/Suzuki (2016), doi:10.1250/ast.37.115. The family
// constants are engineering priors, NOT fits to that paper's measured barrels.
// Stave joints and tensioned hoops are represented by effective properties.
[[nodiscard]] Modes makeModes (const Input&) noexcept;

struct Dashpot
{
    // J = impulsePerVelocity * u, u=sum_i b_i v_i.
    double impulsePerVelocity = 0.0;
    // Positive energy removed: energyLossPerVelocitySquared * u^2.
    double energyLossPerVelocitySquared = 0.0;
};

// Exact flow of M dv/dt = -c b b^T v, with g=b^T M^-1 b.
// q remains fixed. Apply v_i += (b_i/m_i)*J simultaneously to every port.
// The opposite sign of the shell port makes this a reciprocal relative-
// velocity boundary, not a second copy of the membrane's external force.
// Its kinetic-energy change is -(1-exp(-2*c*g*h))*u^2/(2*g), nonpositive
// for every time step; elastic energy is unchanged. Cache at configuration.
[[nodiscard]] Dashpot makeDashpot (double inverseBoundaryMass,
                                    double dampingNsPerMetre,
                                    double stepSeconds) noexcept;

// The pole stores q and its previous sample. Its quadrature coefficients
// recover the exact damped sinusoid's velocity; these helpers also make an
// instantaneous velocity impulse without changing current displacement.
[[nodiscard]] double modalVelocity (double omega, double amplitudeDecay,
                                     double quadratureCurrent,
                                     double quadraturePrevious,
                                     double current, double previous) noexcept;
[[nodiscard]] double velocityToPreviousScale (double omega,
                                              double quadraturePrevious) noexcept;
} // namespace taikor::shellboundary
