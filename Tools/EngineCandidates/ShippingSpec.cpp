#include "ShippingSpec.h"

#include "DSP/BachiModel.h"
#include "DSP/EnsembleEngine.h"
#include "DSP/PhysicalDrumProfile.h"
#include "DSP/ShellBoundary.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <sstream>
#include <stdexcept>

namespace taikor
{
// The engine declares this friend for its tests. The candidate tools never
// link a Tests/ translation unit, so this is the only definition in their
// link and it may read the resolved drum and the armed stroke directly.
struct TaikoEngineTestAccess
{
    static candidates::ResolvedShipping resolve (const candidates::ShippingStroke& take)
    {
        using namespace candidates;
        const auto octave = octaveOffsetForMidiNote (take.midiNote);
        const auto articulation = articulationForMidiNote (take.midiNote);
        if (! octave || ! articulation)
            throw std::runtime_error ("the MIDI note does not map to a stroke");

        const EngineParameters sanitised = TaikoEngine::sanitise (take.parameters);
        const EngineParameters applied = TaikoEngine::parametersForOctave (sanitised, *octave);
        const auto heads = physical::resolveHeads (applied.physicalFamily,
                                                   applied.physicalFamilyMix,
                                                   applied.headMaterial);
        const auto drum = TaikoEngine::resolveDrumFor (take.parameters, 0.0f, *octave);
        const auto shell = shellboundary::makeModes ({ drum.radius, drum.depth,
                                                       applied.shellMaterial,
                                                       drum.physicalFamily,
                                                       drum.physicalFamilyMix });

        ResolvedShipping result;
        auto& d = result.drum;
        d.octave = *octave;
        d.radius = drum.radius;
        d.depth = drum.depth;
        d.shellWallThickness = shell.wallThicknessMetres;

        const double a2 = static_cast<double> (drum.radius) * drum.radius;
        const auto fillHead = [a2] (HeadSpec& head, const physical::Skin& skin,
                                    double tension, double density, double thickness,
                                    double stiffnessB, double edge, double hysteretic,
                                    double viscous)
        {
            head.tension = tension;
            head.arealDensity = density;
            head.thickness = thickness;
            head.youngsModulus = skin.youngsModulus;
            head.poissonRatio = skin.poissonRatio;
            // B = D / (T a^2) is what the shipping modes are built from.
            head.bendingStiffness = stiffnessB * tension * a2;
            head.edgeLoss = edge;
            head.hysteretic = hysteretic;
            head.viscous = viscous;
            fitTwoParameterLoss (head);
        };
        fillHead (d.batter, heads.batter, drum.tension, drum.batterDensity,
                  drum.batterThickness, drum.stiffnessBatter, drum.edgeLoss,
                  drum.headLossFactor, drum.headViscousFactor);
        fillHead (d.rear, heads.rear, drum.resonantTension, drum.resonantDensity,
                  drum.resonantThickness, drum.stiffnessResonant, drum.rearEdgeLoss,
                  drum.rearHeadLossFactor, drum.rearHeadViscousFactor);

        for (int side = 0; side < 2; ++side)
        {
            const double angle = side == 0 ? drum.micAngleLeft : drum.micAngleRight;
            result.mics.capsule[static_cast<std::size_t> (side)] = {
                drum.micRadius * std::cos (angle), drum.micRadius * std::sin (angle),
                drum.micDistanceMetres };
        }

        // Arm the same stroke on a fresh engine and read what it armed.
        auto engine = std::make_unique<TaikoEngine>();
        engine->prepare (take.sampleRate, take.blockSize);
        engine->setParameters (take.parameters);
        engine->reset();
        if (! engine->triggerMidi (take.midiNote, take.velocity))
            throw std::runtime_error ("the shipping engine did not arm the stroke");
        const TaikoEngine::Voice* stroke = nullptr;
        for (const auto& voice : engine->voices_)
            if (voice.active && ! voice.physicalBank
                && (stroke == nullptr || voice.startOrder > stroke->startOrder))
                stroke = &voice;
        if (stroke == nullptr)
            throw std::runtime_error ("no armed stick contact found");
        DrumVisualState visual;
        engine->getVisualState (visual);

        auto& s = result.stroke;
        s.impactSpeed = (stroke->stickPosition - stroke->stickPrevious) * engine->sampleRate_;
        s.strikeX = visual.strikeRadius * drum.radius * std::cos (visual.strikeAngle);
        s.strikeY = visual.strikeRadius * drum.radius * std::sin (visual.strikeAngle);
        s.stickMass = stroke->stickMass;
        s.stickLength = drum.bachi.lengthMetres;
        const double inertia = s.stickMass * s.stickLength * s.stickLength / 12.0;
        const double lever = 0.5 * s.stickLength;
        s.tipMass = s.stickMass / (1.0 + s.stickMass * lever * lever / inertia);
        s.footprintRadius = bachi::contactRadius (drum.bachi, drum.radius, visual.strikeRadius);
        s.hertzStiffness = stroke->contactStiffness;
        s.hertzExponent = stroke->contactExponent;

        std::ostringstream report;
        report.precision (6);
        const auto head = [&report] (const char* name, const HeadSpec& h)
        {
            report << "- " << name << ": T " << h.tension << " N/m, sigma " << h.arealDensity
                   << " kg/m^2, h " << h.thickness * 1000.0 << " mm, E " << h.youngsModulus
                   << " Pa, nu " << h.poissonRatio << ", D " << h.bendingStiffness
                   << " N m; loss edge " << h.edgeLoss << " 1/s + " << h.hysteretic
                   << " omega/2 + " << h.viscous << " s omega^2 -> sigma0 " << h.sigma0
                   << " 1/s, sigma1 " << h.sigma1 << " m^2/s\n";
        };
        report << "- drum: octave " << d.octave << ", radius " << d.radius << " m, depth "
               << d.depth << " m, rigid wall " << d.shellWallThickness * 1000.0 << " mm, air "
               << d.airDensity << " kg/m^3 at " << d.soundSpeed << " m/s\n";
        head ("batter", d.batter);
        head ("rear", d.rear);
        report << "- stroke: " << s.impactSpeed << " m/s at (" << s.strikeX << ", " << s.strikeY
               << ") m; bachi " << s.stickMass << " kg, " << s.stickLength
               << " m, rigid-rod tip mass " << s.tipMass << " kg; footprint "
               << s.footprintRadius * 1000.0 << " mm; Hertz K " << s.hertzStiffness
               << " N/m^" << s.hertzExponent << "\n";
        for (int side = 0; side < 2; ++side)
        {
            const auto& c = result.mics.capsule[static_cast<std::size_t> (side)];
            report << "- " << (side == 0 ? "left" : "right") << " capsule: (" << c[0] << ", "
                   << c[1] << ", " << c[2] << ") m\n";
        }
        result.report = report.str();
        return result;
    }
};
} // namespace taikor

namespace taikor::candidates
{

EngineParameters listeningTestParameters()
{
    EngineParameters parameters;
    parameters.humanise = 0.0f;
    parameters.octaveBody = 1.0f;
    // The plug-in's factory Mic Distance (16 cm on its 3..40 cm range) and
    // Output (-22.5 dB), which EngineParameters' own defaults round.
    parameters.micDistance = (16.0f - 3.0f) / (40.0f - 3.0f);
    parameters.outputGain = std::pow (10.0f, -22.5f / 20.0f);
    return parameters;
}

ResolvedShipping resolveShipping (const ShippingStroke& take)
{
    return TaikoEngineTestAccess::resolve (take);
}

double stiffMembraneXi (const HeadSpec& head, double omega) noexcept
{
    const double T = head.tension, D = head.bendingStiffness, s = head.arealDensity;
    if (D <= 0.0)
        return s * omega * omega / T;
    return (-T + std::sqrt (T * T + 4.0 * D * s * omega * omega)) / (2.0 * D);
}

void fitTwoParameterLoss (HeadSpec& head)
{
    constexpr int points = 64;
    constexpr double low = 60.0, high = 6000.0;
    double aa = 0.0, ab = 0.0, bb = 0.0, a1 = 0.0, b1 = 0.0;
    for (int i = 0; i < points; ++i)
    {
        const double f = low * std::pow (high / low, i / (points - 1.0));
        const double omega = 2.0 * 3.14159265358979323846 * f;
        const double gamma = head.edgeLoss + 0.5 * head.hysteretic * omega
                             + head.viscous * omega * omega;
        const double a = 1.0 / gamma, b = stiffMembraneXi (head, omega) / gamma;
        aa += a * a;
        ab += a * b;
        bb += b * b;
        a1 += a;
        b1 += b;
    }
    const double det = aa * bb - ab * ab;
    double s0 = (a1 * bb - b1 * ab) / det;
    double s1 = (aa * b1 - ab * a1) / det;
    if (s0 < 0.0)
    {
        s0 = 0.0;
        s1 = b1 / bb;
    }
    if (s1 < 0.0)
    {
        s1 = 0.0;
        s0 = a1 / aa;
    }
    head.sigma0 = s0;
    head.sigma1 = s1;
}

namespace
{
void renderBlocks (Take& take, int& position, int until, int blockSize,
                   const std::function<void (float*, float*, int)>& render)
{
    while (position < until)
    {
        const int count = std::min (blockSize, until - position);
        render (take.left.data() + position, take.right.data() + position, count);
        position += count;
    }
}

Take silence (const ShippingStroke& stroke)
{
    Take take;
    take.left.assign (static_cast<std::size_t> (stroke.totalSamples), 0.0f);
    take.right.assign (static_cast<std::size_t> (stroke.totalSamples), 0.0f);
    return take;
}
} // namespace

Take renderShipping (const ShippingStroke& stroke)
{
    auto engine = std::make_unique<EnsembleEngine>();
    engine->setParameters (stroke.parameters);
    engine->prepare (stroke.sampleRate, stroke.blockSize);
    Take take = silence (stroke);
    int position = 0;
    const auto render = [&engine] (float* l, float* r, int n) { engine->process (l, r, n); };
    renderBlocks (take, position, stroke.strikeSample, stroke.blockSize, render);
    if (! engine->triggerMidi (stroke.midiNote, stroke.velocity))
        throw std::runtime_error ("the shipping engine ignored the note");
    renderBlocks (take, position, stroke.totalSamples, stroke.blockSize, render);
    return take;
}

Take renderShippingRaw (const ShippingStroke& stroke)
{
    auto engine = std::make_unique<TaikoEngine>();
    engine->prepare (stroke.sampleRate, stroke.blockSize);
    engine->setParameters (stroke.parameters);
    engine->reset();
    Take take = silence (stroke);
    int position = 0;
    const auto render = [&engine] (float* l, float* r, int n) { engine->processRaw (l, r, n); };
    renderBlocks (take, position, stroke.strikeSample, stroke.blockSize, render);
    if (! engine->triggerMidi (stroke.midiNote, stroke.velocity))
        throw std::runtime_error ("the shipping engine ignored the note");
    renderBlocks (take, position, stroke.totalSamples, stroke.blockSize, render);
    return take;
}

Take applyShippingOutputStage (const ShippingStroke& stroke, const Take& raw)
{
    auto stage = std::make_unique<TaikoEngine>();
    stage->prepare (stroke.sampleRate, stroke.blockSize);
    stage->setParameters (stroke.parameters);
    stage->reset();
    Take take = silence (stroke);
    int position = 0;
    const auto render = [&stage, &raw, &position] (float* l, float* r, int n)
    {
        const auto offset = static_cast<std::size_t> (position);
        stage->processWithEnsemble (l, r, n, raw.left.data() + offset,
                                    raw.right.data() + offset, nullptr, true, nullptr);
    };
    renderBlocks (take, position, stroke.totalSamples, stroke.blockSize, render);
    return take;
}

} // namespace taikor::candidates
