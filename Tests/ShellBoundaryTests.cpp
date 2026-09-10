#include "DSP/ShellBoundary.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>

namespace
{
using namespace taikor::shellboundary;
constexpr double pi = 3.1415926535897932384626433832795;
int failures = 0;

void expect (bool condition, const char* message)
{
    if (! condition)
    {
        if (failures < 12)
            std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

bool near (double actual, double expected, double tolerance = 1.0e-11)
{
    return std::abs (actual - expected)
        <= tolerance * std::max ({ 1.0, std::abs (actual), std::abs (expected) });
}

void passiveReciprocalImpulse()
{
    // Unequal physical masses exchange opposite impulses. Starting either
    // side proves this is a reciprocal coupling, rather than a duplicated
    // head excitation feeding an independent wooden resonator.
    const double headMass = 0.4, shellMass = 13.0;
    const double g = 1.0 / headMass + 1.0 / shellMass;
    for (double step : { 1.0 / 384000.0, 1.0 / 48000.0, 1.0 / 8000.0, 0.2, 100.0 })
    for (double damping : { 0.0, 0.1, 30.0, 10000.0, 1.0e200 })
    for (const auto initial : { std::array<double, 2> { 2.0, 0.0 },
                               std::array<double, 2> { 0.0, 2.0 } })
    {
        const auto damper = makeDashpot (g, damping, step);
        const double relative = initial[0] - initial[1];
        const double impulse = damper.impulsePerVelocity * relative;
        const double head = initial[0] + impulse / headMass;
        const double shell = initial[1] - impulse / shellMass;
        const double before = 0.5 * (headMass * initial[0] * initial[0]
                                      + shellMass * initial[1] * initial[1]);
        const double after = 0.5 * (headMass * head * head + shellMass * shell * shell);
        expect (after <= before + 1.0e-10, "passive for all positive time steps and strengths");
        expect (near (before - after, damper.energyLossPerVelocitySquared * relative * relative),
                "exact energy audit matches state update");
        expect (near (headMass * head + shellMass * shell,
                      headMass * initial[0] + shellMass * initial[1]),
                "equal opposite impulses preserve total momentum");
        if (damping > 0.0)
            expect (initial[0] == 0.0 ? head > 0.0 : shell > 0.0,
                    "energy transfers in both directions");
    }
}

void distributedPorts()
{
    std::mt19937 generator (42);
    std::uniform_real_distribution<double> random (-1.0, 1.0);
    // A generalized rim port contains several signed radial projections.
    // The modal coordinates are not physical point masses, so this tests
    // generalized energy/reciprocity instead of an invalid momentum sum.
    for (int trial = 0; trial < 2000; ++trial)
    {
        std::array<double, 12> velocity {}, inverseMass {}, projection {};
        double before = 0.0, relative = 0.0, g = 0.0;
        for (std::size_t i = 0; i < velocity.size(); ++i)
        {
            velocity[i] = random (generator);
            projection[i] = random (generator);
            inverseMass[i] = std::exp (5.0 * random (generator));
            before += 0.5 * velocity[i] * velocity[i] / inverseMass[i];
            relative += projection[i] * velocity[i];
            g += projection[i] * projection[i] * inverseMass[i];
        }
        const auto damper = makeDashpot (g, 50.0, 1.0 / 48000.0);
        const double impulse = damper.impulsePerVelocity * relative;
        double after = 0.0;
        for (std::size_t i = 0; i < velocity.size(); ++i)
        {
            velocity[i] += projection[i] * inverseMass[i] * impulse;
            after += 0.5 * velocity[i] * velocity[i] / inverseMass[i];
        }
        expect (after <= before + 1.0e-11, "signed distributed ports cannot inject energy");
        expect (near (before - after, damper.energyLossPerVelocitySquared * relative * relative),
                "distributed energy loss equals exact dissipation");
    }
    const auto full = makeDashpot (2.7, 30.0, 0.01);
    const auto half = makeDashpot (2.7, 30.0, 0.005);
    expect (near (1.0 + 2.7 * full.impulsePerVelocity,
                  std::pow (1.0 + 2.7 * half.impulsePerVelocity, 2.0)),
            "exact damper has sample-rate-independent semigroup");
}

void dampedPoleConversion()
{
    for (double rate : { 8000.0, 44100.0, 48000.0, 192000.0, 384000.0 })
    for (double fraction : { 0.00001, 0.01, 0.25, 0.48 })
    for (double decay : { 0.0, 1.0, 500.0 })
    {
        const double omega = 2.0 * pi * fraction * rate;
        const double theta = omega / rate;
        const double radius = std::exp (-decay / rate);
        const double qc = std::cos (theta) / std::sin (theta);
        const double qp = -radius / std::sin (theta);
        const double current = 0.002, initialVelocity = 0.7;
        const double previous = ((initialVelocity + decay * current) / omega
                                 - qc * current) / qp;
        expect (near (modalVelocity (omega, decay, qc, qp, current, previous), initialVelocity),
                "damped pole recovers physical velocity");
        const double increment = -0.13;
        const double adjusted = previous + increment * velocityToPreviousScale (omega, qp);
        expect (near (modalVelocity (omega, decay, qc, qp, current, adjusted),
                      initialVelocity + increment, 1.0e-9),
                "velocity kick changes velocity without changing displacement");
    }
}

void constructionAndScaling()
{
    const auto carved = makeModes ({ 0.275, 0.5, 0.8, 0, 1.0 });
    const auto enlarged = makeModes ({ 0.55, 1.0, 0.8, 0, 1.0 });
    const auto tub = makeModes ({ 0.275, 0.5, 0.2, 2, 1.0 });
    const auto hoop = makeModes ({ 0.275, 0.5, 0.92, 3, 1.0 });
    for (std::size_t i = 0; i < modeCount; ++i)
    {
        expect (near (enlarged.modalMassKg[i] / carved.modalMassKg[i], 8.0),
                "geometrically scaled carved shell mass follows volume");
        expect (near (enlarged.frequencyHz[i] / carved.frequencyHz[i], 0.5),
                "geometrically scaled shell frequency follows inverse size");
        expect (tub.modalMassKg[i] < carved.modalMassKg[i] * 0.2,
                "light staves are mechanically distinct from a carved barrel");
        expect (tub.frequencyHz[i] < carved.frequencyHz[i],
                "jointed thin shell has different effective flexural modes");
        expect (hoop.modalMassKg[i] > tub.modalMassKg[i], "tensioned ring hardware has inertia");
    }
    for (int family = 0; family < 4; ++family)
    for (double material : { 0.0, 0.5, 1.0 })
    for (double radius : { 0.008, 0.2, 0.75, 3.75 })
    {
        const auto modes = makeModes ({ radius, 1.1 * radius, material, family, 1.0 });
        expect (modes.transverseYoungsModulusPa < modes.axialYoungsModulusPa,
                "wood retains separate longitudinal and cross-grain stiffness");
        double previousFrequency = 0.0;
        for (std::size_t i = 0; i < modeCount; ++i)
        {
            expect (std::isfinite (modes.frequencyHz[i]) && modes.frequencyHz[i] > previousFrequency,
                    "ordered positive shell frequencies over supported geometry");
            expect (modes.modalMassKg[i] > 0.0 && modes.amplitudeDecay[i] > 0.0,
                    "positive physical modal mass and loss");
            previousFrequency = modes.frequencyHz[i];
        }
    }
    const auto collapsedFamily = makeModes ({ 0.275, 0.5, 0.8, 2, 0.0 });
    expect (collapsedFamily.frequencyHz == carved.frequencyHz, "zero family mix selects the base construction");
    expect (makeDashpot (0.0, 1.0, 0.1).impulsePerVelocity == 0.0, "zero mobility is bypassed");
    expect (makeDashpot (1.0, -1.0, 0.1).impulsePerVelocity == 0.0, "negative damping is rejected");
    expect (makeDashpot (1.0, 1.0, std::numeric_limits<double>::quiet_NaN()).impulsePerVelocity == 0.0,
            "invalid step is rejected");
}
} // namespace

int main()
{
    passiveReciprocalImpulse();
    distributedPorts();
    dampedPoleConversion();
    constructionAndScaling();
    if (failures != 0)
        return 1;
    std::cout << "Shell boundary passivity, reciprocity, pole conversion and construction tests passed.\n";
}
