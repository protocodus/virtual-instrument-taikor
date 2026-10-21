#pragma once

// The contract every alternative drum engine in the engine-candidate listening
// test implements. A candidate replaces the shipping engine's drum model and
// nothing else: the drum it models, the stroke it plays and the two capsules it
// is heard through all arrive here, resolved from the shipping engine's own
// state by ShippingSpec.cpp, and its microphone pressure is then passed through
// the shipping output stage by RenderEngineCandidates.cpp. A candidate must not
// substitute its own values for anything these specs carry.

#include <array>
#include <memory>
#include <string>

namespace taikor::candidates
{

// One head as the shipping engine resolves it for the factory controls.
struct HeadSpec
{
    double tension = 0.0;          // N/m
    double arealDensity = 0.0;     // kg/m^2
    double thickness = 0.0;        // m
    double youngsModulus = 0.0;    // Pa
    double poissonRatio = 0.3;
    double bendingStiffness = 0.0; // D = E h^3 / 12(1 - nu^2), N m

    // The shipping engine's non-radiative hide and rim loss, as an amplitude
    // decay rate: gamma(omega) = edgeLoss + hysteretic * omega / 2 + viscous * omega^2.
    // These are the shipping engine's priors (PhysicalDrumProfile.cpp and
    // TaikoEngine::resolveDrumFor), not measurements. Radiation and mounting
    // are deliberately absent: every candidate computes radiation itself, and
    // the shipping mount loss is a voicing term standing in for it.
    double edgeLoss = 0.0;         // 1/s
    double hysteretic = 0.0;       // loss angle, damping shape included
    double viscous = 0.0;          // s

    // The same law on the standard two-parameter form gamma = sigma0 + sigma1 xi,
    // where xi(omega) is the squared wavenumber of the stiff membrane at omega
    // (Bilbao, Numerical Sound Synthesis, 2009, sec. 11.4). A finite-difference
    // head can represent only this form, so every candidate uses this pair and
    // the hide loss is identical across letters. See fitTwoParameterLoss.
    double sigma0 = 0.0;           // 1/s
    double sigma1 = 0.0;           // m^2/s
};

// The drum: two clamped heads on a rigid barrel, axis along z.
struct DrumSpec
{
    int octave = 1;
    double radius = 0.0;           // head radius a, m
    double depth = 0.0;            // barrel length between the heads, m
    double shellWallThickness = 0.0; // m; the barrel is rigid in every candidate
    HeadSpec batter, rear;
    double airDensity = 1.2041;    // kg/m^3 (the shipping engine's constant)
    double soundSpeed = 343.0;     // m/s (the shipping engine's constant)
};

// The stroke. Coordinates are in the batter plane with the head centre at the
// origin; +z points out of the drum, towards the microphones.
struct StrokeSpec
{
    double impactSpeed = 0.0;      // m/s into the head at first contact
    double strikeX = 0.0, strikeY = 0.0; // m
    double stickMass = 0.0;        // kg, the whole bachi (shipping BachiModel)
    double stickLength = 0.0;      // m
    // Effective mass of a free rigid rod struck at its tip, m / (1 + m (L/2)^2 / I_G)
    // with I_G = m L^2 / 12, i.e. m / 4. A grip pivot about 0.1 m from the butt
    // gives 0.26 m, so the free-rod value stands for a held stick too. The
    // shipping engine treats the whole mass as a point at the tip.
    double tipMass = 0.0;          // kg
    double footprintRadius = 0.0;  // m, the tip disk that senses and loads the head
    double hertzStiffness = 0.0;   // N/m^alpha, the shipping stroke's own K
    double hertzExponent = 1.5;
};

// The close pair. capsule[0] is the left channel, capsule[1] the right.
struct MicSpec
{
    std::array<std::array<double, 3>, 2> capsule {};
};

class Engine
{
public:
    virtual ~Engine() = default;

    // One line naming the model and its references. It is written to key.md
    // only, never shown to the listener before they decide.
    [[nodiscard]] virtual std::string description() const = 0;

    // Build the drum for this sample rate. May take minutes and use threads.
    // Returns false and sets error if the drum cannot be built.
    virtual bool prepare (const DrumSpec& drum, const MicSpec& mics,
                          double sampleRate, std::string& error) = 0;

    // Arm the stroke. The stick meets the batter during the next processed
    // sample: it starts on the head surface moving in at impactSpeed.
    virtual void strike (const StrokeSpec& stroke) = 0;

    // Microphone pressure in pascals at the two capsules. Any block length.
    virtual void process (float* left, float* right, int samples) = 0;

    // Free-form diagnostics for key.md: contact duration, peak force, impulse,
    // energy balance, render cost, anything a reviewer needs to trust the take.
    [[nodiscard]] virtual std::string diagnostics() const = 0;
};

class NonlinearPotential;

// The finite-difference drum (Fdtd/). The batter's geometric nonlinearity is
// injected so the Berger and von Karman letters differ in that object alone;
// the rear head always uses Berger. Returns nullptr until implemented.
[[nodiscard]] std::unique_ptr<Engine> makeFdtdEngine (
    std::unique_ptr<NonlinearPotential> batterNonlinearity);
// The complete coupled-mode drum (Modal/). Returns nullptr until implemented.
[[nodiscard]] std::unique_ptr<Engine> makeModalEngine();

} // namespace taikor::candidates
