#include "DSP/ShellBoundary.h"

#include <algorithm>
#include <cmath>

namespace taikor::shellboundary
{
namespace
{
constexpr double pi = 3.1415926535897932384626433832795;

struct Construction
{
    double axialModulus, transverseModulus, density, wallRatio;
    double lossAngle, boundaryDamping, barrelRadiusRatio, hoopMass;
    double referenceControl, observationCalibration;
};

// Effective dry-wood construction priors, not measured replicas. Separating
// the cross-grain modulus avoids treating dense carved wood as isotropic.
// The light tub includes joint losses; the shime includes the ring hardware's
// inertial load. No result here establishes the sound of a particular maker.
constexpr std::array<Construction, 4> constructions {{
    { 11.0e9, 1.10e9, 780.0, 0.070, 0.040, 4.0, 1.16, 0.0, 0.80, 1100.0 },
    { 11.0e9, 1.10e9, 780.0, 0.067, 0.044, 3.0, 1.14, 0.0, 0.74, 1000.0 },
    {  7.0e9, 0.42e9, 390.0, 0.015, 0.095, 1.5, 1.00, 0.12, 0.20,  150.0 },
    { 10.0e9, 0.95e9, 730.0, 0.065, 0.060, 5.0, 1.02, 0.80, 0.92, 1500.0 }
}};

double finiteOr (double value, double fallback) noexcept
{
    return std::isfinite (value) ? value : fallback;
}
} // namespace

Modes makeModes (const Input& input) noexcept
{
    Modes result;
    const auto& selected = constructions[static_cast<std::size_t> (
        std::clamp (input.family, 0, 3))];
    const auto& reference = constructions[0];
    const double mix = std::clamp (finiteOr (input.familyMix, 1.0), 0.0, 1.0);
    const auto blend = [mix] (double a, double b) noexcept
    { return a + mix * (b - a); };
    result.observationCalibration = blend (reference.observationCalibration,
                                           selected.observationCalibration);
    const double radius = std::clamp (finiteOr (input.headRadiusMetres, 0.275),
                                      0.008, 3.75);
    const double depth = std::clamp (finiteOr (input.depthMetres, 0.5), 0.04, 9.5);
    const double material = std::clamp (finiteOr (input.shellMaterial, 0.8), 0.0, 1.0);
    const double referenceControl = blend (reference.referenceControl, selected.referenceControl);
    const double trim = material - referenceControl;
    result.axialYoungsModulusPa = blend (reference.axialModulus, selected.axialModulus)
                                  * std::exp2 (1.4 * trim);
    result.transverseYoungsModulusPa = blend (reference.transverseModulus, selected.transverseModulus)
                                       * std::exp2 (1.8 * trim);
    result.woodDensityKgPerM3 = blend (reference.density, selected.density)
                                * std::exp2 (0.35 * trim);
    result.wallThicknessMetres = 2.0 * radius
        * blend (reference.wallRatio, selected.wallRatio) * std::exp2 (0.25 * trim);
    const double shellRadius = radius
        * blend (reference.barrelRadiusRatio, selected.barrelRadiusRatio);
    const double hoopMass = blend (reference.hoopMass, selected.hoopMass)
                           * radius / 0.15;

    // Orthotropic plane-stress reciprocity: nu_theta,z / E_theta
    // = nu_z,theta / E_z. A ring with no axial curvature bends against D_theta;
    // longitudinal stiffness still enters the reciprocal Poisson constraint.
    constexpr double nuAxialTransverse = 0.30;
    const double nuTransverseAxial = nuAxialTransverse
        * result.transverseYoungsModulusPa / result.axialYoungsModulusPa;
    const double reciprocalPoisson = 1.0 - nuAxialTransverse * nuTransverseAxial;
    const double thickness = result.wallThicknessMetres;
    const double flexuralRigidity = result.transverseYoungsModulusPa
        * thickness * thickness * thickness / (12.0 * reciprocalPoisson);
    const double wallMass = 2.0 * pi * shellRadius * depth * thickness
                            * result.woodDensityKgPerM3;
    const double totalMass = wallMass + hoopMass;
    const double loss = blend (reference.lossAngle, selected.lossAngle)
                        * std::exp2 (-0.75 * trim);
    const double boundary = blend (reference.boundaryDamping, selected.boundaryDamping)
                            * radius / 0.275;

    for (std::size_t i = 0; i < modeCount; ++i)
    {
        const double n = static_cast<double> (result.circumferentialOrder[i]);
        // w=q cos(n theta), u=-q sin(n theta)/n preserves circumference.
        // T=1/2 * (wall+hoops)/2 * (1+1/n^2) * qdot^2.
        // V=1/2 * pi*L*D_theta/R^3 * (n^2-1)^2 * q^2.
        const double mass = 0.5 * totalMass * (1.0 + 1.0 / (n * n));
        const double stiffness = pi * depth * flexuralRigidity
            / (shellRadius * shellRadius * shellRadius)
            * (n * n - 1.0) * (n * n - 1.0);
        const double omega = std::sqrt (stiffness / mass);
        result.frequencyHz[i] = omega / (2.0 * pi);
        result.modalMassKg[i] = mass;
        result.amplitudeDecay[i] = 0.5 * omega * loss;
        result.boundaryDampingNsPerMetre[i] = boundary / (1.0 + 0.2 * static_cast<double> (i));
    }
    return result;
}

Dashpot makeDashpot (double inverseBoundaryMass, double dampingNsPerMetre,
                      double stepSeconds) noexcept
{
    if (! std::isfinite (inverseBoundaryMass) || inverseBoundaryMass <= 0.0
        || ! std::isfinite (dampingNsPerMetre) || dampingNsPerMetre <= 0.0
        || ! std::isfinite (stepSeconds) || stepSeconds <= 0.0)
        return {};
    const double exponent = dampingNsPerMetre * inverseBoundaryMass * stepSeconds;
    const double changedFraction = -std::expm1 (-exponent);
    const double lossFraction = -std::expm1 (-2.0 * exponent);
    return { -changedFraction / inverseBoundaryMass,
             0.5 * lossFraction / inverseBoundaryMass };
}

double modalVelocity (double omega, double amplitudeDecay,
                       double quadratureCurrent, double quadraturePrevious,
                       double current, double previous) noexcept
{
    return omega * (quadratureCurrent * current + quadraturePrevious * previous)
           - amplitudeDecay * current;
}

double velocityToPreviousScale (double omega, double quadraturePrevious) noexcept
{
    const double denominator = omega * quadraturePrevious;
    return std::isfinite (denominator) && std::abs (denominator) > 1.0e-12
        ? 1.0 / denominator : 0.0;
}
} // namespace taikor::shellboundary
