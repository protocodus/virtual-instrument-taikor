#include "DSP/TaikoEngine.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
int failures = 0;
void expect (bool condition, const std::string& message)
{
    if (! condition)
    {
        if (failures < 24)
            std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

double relativeError (double actual, double expected, double floor = 1.0e-12)
{
    return std::abs (actual - expected)
         / std::max ({ floor, std::abs (actual), std::abs (expected) });
}

// Independent static solve for the integrated physical K matrix. Comparing
// this displacement with a modal force/readout sum catches missing radial
// coordinates, wrong modal mass and inconsistent front/rear projections.
taikor::cavity::Vector solveStatic (taikor::cavity::Matrix matrix,
                                   taikor::cavity::Vector force, std::size_t count)
{
    for (std::size_t pivot = 0; pivot < count; ++pivot)
    {
        std::size_t row = pivot;
        for (std::size_t other = pivot + 1; other < count; ++other)
            if (std::abs (matrix[other][pivot]) > std::abs (matrix[row][pivot]))
                row = other;
        std::swap (matrix[pivot], matrix[row]);
        std::swap (force[pivot], force[row]);
        expect (std::abs (matrix[pivot][pivot]) > 1.0e-16,
                "integrated shared stiffness must have a nonsingular static solve");
        if (std::abs (matrix[pivot][pivot]) <= 1.0e-16)
            return {};
        for (std::size_t other = pivot + 1; other < count; ++other)
        {
            const double factor = matrix[other][pivot] / matrix[pivot][pivot];
            for (std::size_t column = pivot; column < count; ++column)
                matrix[other][column] -= factor * matrix[pivot][column];
            force[other] -= factor * force[pivot];
        }
    }
    taikor::cavity::Vector result {};
    for (std::size_t reverse = count; reverse > 0; --reverse)
    {
        const auto row = reverse - 1;
        double value = force[row];
        for (std::size_t column = row + 1; column < count; ++column)
            value -= matrix[row][column] * result[column];
        result[row] = value / matrix[row][row];
    }
    return result;
}

double dot (const taikor::cavity::Vector& a, const taikor::cavity::Vector& b)
{
    double result = 0.0;
    for (std::size_t index = 0; index < a.size(); ++index)
        result += a[index] * b[index];
    return result;
}
} // namespace

namespace taikor
{
// This executable owns its own friend accessor; no production API is added.
struct TaikoEngineTestAccess
{
    static void checkIntegratedBasis (int family, float coupling)
    {
        EngineParameters parameters;
        parameters.humanise = 0.0f;
        parameters.cavityCoupling = coupling;
        const auto drum = TaikoEngine::resolveDrumFor (parameters, 0.0f, family);
        const std::string where = " (family " + std::to_string (family)
                                + ", coupling " + std::to_string (coupling) + ")";
        expect (drum.sharedAirValid, "factory engine must activate the shared cavity" + where);
        if (! drum.sharedAirValid)
            return;
        const auto& solution = drum.sharedAir;
        expect (solution.coordinateCount > cavity::headCoordinateCount,
                "integrated solve must retain explicit air coordinates" + where);

        TaikoEngine engine;
        engine.prepare (48000.0, 256);
        const auto& profile = TaikoEngine::strikeProfile (Articulation::Don);
        auto voice = std::make_unique<TaikoEngine::Voice>();
        voice->physicalBank = true;
        voice->strikeRadius = 0.41f;
        voice->strikeAngle = 0.37f;
        engine.buildVoiceModes (*voice, drum, profile, 0.0f, false);
        std::array<bool, TaikoEngine::resonatorCount> ids {};
        int sharedCount = 0;
        for (int index = 0; index < voice->modeCount; ++index)
        {
            const auto& mode = voice->modes[static_cast<std::size_t> (index)];
            const auto id = static_cast<std::size_t> (mode.physicalIndex);
            expect (id < ids.size() && ! ids[id],
                    "every integrated pole must have a distinct bounded physical ID" + where);
            if (id < ids.size())
                ids[id] = true;
            if (! mode.sharedCavityMode)
                continue;
            ++sharedCount;
            expect (mode.inverseModalMass == 1.0f,
                    "mixed cavity modes must retain unit modal mass in the renderer" + where);
            double modalMass = 0.0;
            for (std::size_t i = 0; i < solution.coordinateCount; ++i)
                for (std::size_t j = 0; j < solution.coordinateCount; ++j)
                    modalMass += mode.cavityBasis[i] * solution.mass[i][j] * mode.cavityBasis[j];
            expect (relativeError (modalMass, 1.0) < 2.0e-6,
                    "rendered full cavity basis must be normalized in physical mass" + where);

            double expectedContact = 0.0;
            for (std::size_t radial = 0; radial < cavity::radialModeCount; ++radial)
            {
                const double lambda = TaikoEngine::membraneModes()[radial].besselZero;
                const double footprint = bachi::contactRadius (drum.bachi, drum.radius, voice->strikeRadius);
                expectedContact += TaikoEngine::besselJ (0, lambda * voice->strikeRadius)
                    * bachi::diskProjection (lambda * footprint / drum.radius)
                    * mode.cavityBasis[radial];
            }
            expectedContact *= profile.membraneGain * profile.levelScale;
            expect (relativeError (mode.contactShape, expectedContact, 1.0e-5) < 3.0e-5,
                    "rendered contact must use every batter coordinate of its shared eigenvector" + where);
        }
        int eligibleShared = 0;
        for (std::size_t index = 0; index < solution.coordinateCount; ++index)
        {
            const auto data = TaikoEngine::sharedModeData (drum, static_cast<int> (index),
                                                          voice->strikeRadius, drum.micDistanceMetres);
            if (data.omega / (2.0f * 3.14159265358979f) < TaikoEngine::renderedModeCeilingHz (48000.0)
                && data.omega > 0.0f)
                ++eligibleShared;
        }
        expect (sharedCount > 0 && sharedCount == eligibleShared,
                "factory shared bank must retain every in-band head and air coordinate" + where);

        cavity::Vector forceB {}, forceR {};
        constexpr float radiusB = 0.23f, radiusR = 0.61f;
        for (std::size_t radial = 0; radial < cavity::radialModeCount; ++radial)
        {
            const double lambda = TaikoEngine::membraneModes()[radial].besselZero;
            forceB[radial] = TaikoEngine::besselJ (0, lambda * radiusB)
                * bachi::diskProjection (lambda * bachi::contactRadius (
                      drum.bachi, drum.radius, radiusB) / drum.radius);
            forceR[cavity::radialModeCount + radial] = TaikoEngine::besselJ (0, lambda * radiusR)
                * bachi::diskProjection (lambda * bachi::contactRadius (
                      drum.bachi, drum.radius, radiusR) / drum.radius);
        }
        const auto responseB = solveStatic (solution.stiffness, forceB, solution.coordinateCount);
        const auto responseR = solveStatic (solution.stiffness, forceR, solution.coordinateCount);
        const double directBB = dot (forceB, responseB);
        const double directBR = dot (forceR, responseB);
        const double directRB = dot (forceB, responseR);
        double modalBB = 0.0, modalBR = 0.0;
        for (std::size_t index = 0; index < solution.coordinateCount; ++index)
        {
            const auto b = TaikoEngine::sharedModeData (drum, static_cast<int> (index),
                                                       radiusB, drum.micDistanceMetres);
            const auto r = TaikoEngine::sharedModeData (drum, static_cast<int> (index),
                                                       radiusR, drum.micDistanceMetres, true);
            const double omegaSquared = solution.angularFrequencies[index]
                                      * solution.angularFrequencies[index];
            modalBB += static_cast<double> (b.contact) * b.contact / omegaSquared;
            modalBR += static_cast<double> (b.contact) * r.contact / omegaSquared;
            if (coupling == 0.0f)
                expect (std::abs (static_cast<double> (b.contact) * r.contact) < 1.0e-10,
                        "zero coupling must separate the two head contact ports" + where);
        }
        expect (directBB > 0.0 && relativeError (modalBB, directBB) < 2.0e-5,
                "unit-mass contact compliance must match the independent full stiffness solve" + where);
        expect (relativeError (directBR, directRB, 1.0e-16) < 1.0e-8,
                "the integrated front/rear ports must obey mechanical reciprocity" + where);
        expect (relativeError (modalBR, directBR, 1.0e-11) < 4.0e-5,
                "both head projections must match the full physical cross-compliance" + where);
        if (coupling == 0.0f)
            expect (std::abs (modalBR) < 1.0e-14 && std::abs (directBR) < 1.0e-14,
                    "an open shared cavity must transfer no contact force between heads" + where);
        else
            expect (std::abs (directBR) > 1.0e-9,
                    "reciprocity must be tested on a nontrivially coupled pair" + where);
    }

    static void checkPhysicalRetuning()
    {
        constexpr double pi = 3.14159265358979323846;
        for (int family = 0; family < 4; ++family)
        {
            EngineParameters parameters;
            parameters.humanise = parameters.tensionModulation = 0.0f;
            const auto base = TaikoEngine::resolveDrumFor (parameters, 0.0f, family);
            expect (relativeError (2.0 * base.radius,
                        getDrumDescription (family).headDiameterMetres) < 2.0e-6,
                    "factory family geometry must remain its explicitly declared size");
            for (float coupling : { 0.0f, 0.001f, 0.85f, 1.0f })
            {
                parameters.cavityCoupling = coupling;
                const auto changed = TaikoEngine::resolveDrumFor (parameters, 0.0f, family);
                expect (changed.radius == base.radius && changed.depth == base.depth
                            && changed.tension == base.tension,
                        "Air automation must never retune or resize the physical instrument");
            }
            parameters.cavityCoupling = 0.85f;
            TaikoEngine engine;
            engine.setParameters (parameters);
            engine.prepare (48000.0, 64);
            engine.ensurePhysicalDrum (family, base);
            auto& bank = engine.physicalDrums_[static_cast<std::size_t> (family)];
            constexpr float shift = 1.122462048f; // two semitones of head wave speed
            engine.applyTensionShift (bank, shift);
            for (int index = 0; index < bank.modeCount; ++index)
            {
                const auto& mode = bank.modes[static_cast<std::size_t> (index)];
                if (! mode.membrane || ! (mode.liveOmega > 0.0))
                    continue;
                const double expected = mode.omega * std::sqrt (1.0
                    + (static_cast<double> (shift) * shift - 1.0)
                       * (mode.batterTensionFraction + mode.rearTensionFraction));
                expect (relativeError (mode.liveOmega, expected) < 4.0e-7,
                        "live Pitch must project head tension while leaving air and bending stiffness fixed");
            }
            engine.pitchBend_ = engine.pitchBendTarget_ = 1.0f;
            const auto bent = TaikoEngine::resolveDrumFor (parameters, 2.0f, family);
            engine.ensurePhysicalDrum (family, bent);
            expect (bank.tuningAtStrike == 2.0f,
                    "a new bent strike must refresh the physical bank's tuning baseline");
            for (int index = 0; index < bank.modeCount; ++index)
            {
                const auto& mode = bank.modes[static_cast<std::size_t> (index)];
                if (! mode.membrane)
                    continue;
                const auto observed = TaikoEngine::observeIdentity (
                    bent, TaikoEngine::identityForMode (mode), 0.23f);
                expect (relativeError (mode.omega / (2.0 * pi), observed.frequencyHz) < 4.0e-6,
                        "a refreshed bent bank must contain the fully resolved physical poles");
            }
        }
        for (float layout : { 0.0f, 1.0f })
        {
            EngineParameters parameters;
            parameters.octaveBody = layout;
            float previous = 0.0f;
            for (int family = 0; family < 4; ++family)
            {
                const auto result = TaikoEngine::measure (parameters, family);
                if (family > 0)
                    expect (std::abs (1200.0f * std::log2 (result.soundingHz / previous)
                                      - 1200.0f) < 20.0f,
                            "explicit factory tension priors must preserve the neutral heard-octave ladder");
                previous = result.soundingHz;
                if (layout == 0.0f)
                    expect (result.radiusMetres == 0.75f,
                            "1 Drum retuning must leave the reference drum's dimensions fixed");
            }
        }
    }

    struct PhysicalState { cavity::Vector position {}, velocity {}, pending {}; int modes = 0; };

    static void checkReadoutAndPalmPersistence()
    {
        EngineParameters parameters;
        parameters.humanise = parameters.tensionModulation = 0.0f;
        TaikoEngine engine;
        engine.setParameters (parameters);
        engine.prepare (48000.0, 256);
        const auto drum = TaikoEngine::resolveDrumFor (parameters, 0.0f, 0);
        engine.ensurePhysicalDrum (0, drum);
        auto& bank = engine.physicalDrums_[0];
        bank.active = true;
        bool testedAdditionalAirPole = false;
        for (int index = 0; index < bank.modeCount; ++index)
        {
            const auto& mode = bank.modes[static_cast<std::size_t> (index)];
            if (! mode.membrane)
                continue;
            const auto identity = TaikoEngine::identityForMode (mode);
            expect (identity.branch < 2,
                    "public mode metadata must not decode an extra air/rear ID as a giant branch");
            const auto observed = TaikoEngine::observeIdentity (drum, identity, 0.23f, 0.37f);
            expect (relativeError (observed.frequencyHz, mode.omega / 6.283185307179586) < 3.0e-6,
                    "every explicit shared/rear readout identity must resolve to the rendered pole");
            testedAdditionalAirPole = testedAdditionalAirPole || identity.sharedIndex >= 8;
        }
        expect (testedAdditionalAirPole, "readout identity audit must include the additional air coordinates");
        const auto measured = TaikoEngine::measure (parameters, 0, 0.0f, 48000.0);
        bool foundLowest = false, foundBreathing = false, foundSounding = false;
        for (int index = 0; index < bank.modeCount; ++index)
        {
            const auto& mode = bank.modes[static_cast<std::size_t> (index)];
            const double hz = mode.omega / 6.283185307179586;
            if (! mode.membrane)
                continue;
            foundLowest = foundLowest || relativeError (measured.loadedFundamentalHz, hz) < 3.0e-6;
            foundBreathing = foundBreathing || relativeError (measured.breathingModeHz, hz) < 3.0e-6;
            foundSounding = foundSounding || relativeError (measured.soundingHz, hz) < 3.0e-6;
        }
        expect (measured.sharedCavityActive && foundLowest && foundBreathing && foundSounding,
                "public measurements must name poles present in the active shared bank");
        expect (measured.cavityEnergyFraction > 0.0f && measured.cavityEnergyFraction <= 1.0f
                    && measured.tailSeconds > 0.0f && measured.tailSeconds <= maximumTailSeconds,
                "shared energy participation and structural tail estimate must have physical bounds");

        const auto& palm = TaikoEngine::strikeProfile (Articulation::Tsu);
        engine.dampPhysicalDrum (bank, palm, 0.22f, drum, false);
        engine.dampPhysicalDrum (bank, palm, 0.63f, drum, true);
        const auto ticks = bank.localMuteTicksRemaining;
        parameters.bodyDepth = 0.72f;
        parameters.headMaterial = 0.86f;
        engine.setParameters (parameters);
        const auto changed = TaikoEngine::resolveDrumFor (parameters, 0.0f, 0);
        engine.ensurePhysicalDrum (0, changed);
        expect (bank.localMutePatchCount == 2 && bank.localMuteTicksRemaining == ticks,
                "structural rebuild must retain both physical palm patches and their remaining duration");
        int dampedShared = 0;
        for (int index = 0; index < bank.modeCount; ++index)
        {
            const auto& mode = bank.modes[static_cast<std::size_t> (index)];
            if (! mode.sharedCavityMode)
                continue;
            const float expected = palm.muteAmount * palm.muteAmount * std::max (
                TaikoEngine::sharedPalmDamping (changed, mode, 0.22f, false),
                TaikoEngine::sharedPalmDamping (changed, mode, 0.63f, true));
            expect (relativeError (mode.localMuteDampingRate, expected, 1.0e-6) < 2.0e-6,
                    "held palms must be reprojected through the new complete shared basis");
            dampedShared += mode.localMuteDampingRate > 0.0f ? 1 : 0;
        }
        expect (dampedShared > 0, "the preserved local mute must still damp shared modes");
        engine.handDamping_ = engine.handDampingTarget_ = 0.4f;
        engine.updateVoiceControl (bank);
        for (int index = 0; index < bank.activeModeCount; ++index)
        {
            const auto& mode = bank.modes[static_cast<std::size_t> (index)];
            if (! mode.membrane)
                continue;
            const double expected = 0.5 * (mode.handDampingRate * 0.16f + mode.localMuteDampingRate);
            expect (relativeError (mode.appliedPalmDecay, expected, 1.0e-6) < 2.0e-6,
                    "recomputed local palms and held CC1 must reach live pole damping");
        }
        for (int tick = 0; tick < ticks; ++tick)
            engine.updateVoiceControl (bank);
        expect (bank.localMutePatchCount == 0 && bank.localMuteTicksRemaining == 0,
                "expired local palms must release their physical patch metadata");
    }

    static double velocityOf (const TaikoEngine::Mode& mode)
    {
        const double cosine = -mode.resonator.a1 / (2.0 * mode.poleRadius);
        const double quadrature = (mode.resonator.y1 * cosine
            - mode.poleRadius * mode.resonator.y2) / mode.resonator.b0;
        return mode.liveOmega * quadrature
             - (mode.decayRate + mode.appliedPalmDecay) * mode.resonator.y1;
    }

    static PhysicalState snapshot (const TaikoEngine::Voice& voice)
    {
        PhysicalState state;
        for (int index = 0; index < voice.modeCount; ++index)
        {
            const auto& mode = voice.modes[static_cast<std::size_t> (index)];
            if (! mode.sharedCavityMode)
                continue;
            ++state.modes;
            const double velocity = velocityOf (mode);
            const double pending = mode.resonator.b0 * voice.modalInput[mode.physicalIndex];
            for (std::size_t coordinate = 0; coordinate < state.position.size(); ++coordinate)
            {
                state.position[coordinate] += mode.cavityBasis[coordinate] * mode.resonator.y1;
                state.velocity[coordinate] += mode.cavityBasis[coordinate] * velocity;
                state.pending[coordinate] += mode.cavityBasis[coordinate] * pending;
            }
        }
        return state;
    }

    static void checkStructuralRemap()
    {
        EngineParameters parameters;
        parameters.humanise = parameters.tensionModulation = 0.0f;
        TaikoEngine engine;
        engine.setParameters (parameters);
        engine.prepare (192000.0, 256);
        constexpr int family = 2;
        const auto oldDrum = TaikoEngine::resolveDrumFor (parameters, 0.0f, family);
        engine.ensurePhysicalDrum (family, oldDrum);
        auto& physical = engine.physicalDrums_[family];
        physical.active = true;
        for (int index = 0; index < physical.modeCount; ++index)
        {
            auto& mode = physical.modes[static_cast<std::size_t> (index)];
            if (! mode.sharedCavityMode)
                continue;
            const double displacement = 1.0e-5 * std::cos (0.71 * mode.physicalIndex + 0.4);
            const double velocity = 0.025 * std::sin (0.43 * mode.physicalIndex + 0.3);
            mode.resonator.y1 = displacement;
            const double cosine = -mode.resonator.a1 / (2.0 * mode.poleRadius);
            const double quadrature = (velocity
                + (mode.decayRate + mode.appliedPalmDecay) * displacement) / mode.liveOmega;
            mode.resonator.y2 = (displacement * cosine
                - quadrature * mode.resonator.b0) / mode.poleRadius;
            physical.modalInput[mode.physicalIndex] = static_cast<float> (
                1.0e-7 * std::cos (0.27 * mode.physicalIndex) / mode.resonator.b0);
        }
        const auto before = snapshot (physical);
        const auto revision = physical.configurationRevision;
        parameters.headMaterial = 0.81f;
        parameters.resonantTension = 0.69f;
        parameters.bodyDepth = 0.64f;
        parameters.cavityCoupling = 0.57f;
        engine.setParameters (parameters);
        const auto newDrum = TaikoEngine::resolveDrumFor (parameters, 0.0f, family);
        engine.ensurePhysicalDrum (family, newDrum);
        const auto after = snapshot (physical);
        expect (physical.configurationRevision != revision && physical.active,
                "structural automation must rebuild the active canonical drum");
        expect (before.modes == static_cast<int> (oldDrum.sharedAir.coordinateCount)
                    && after.modes == static_cast<int> (newDrum.sharedAir.coordinateCount),
                "remap audit must retain the complete old and new head/air bases");
        expect (std::abs (before.position[cavity::headCoordinateCount]) > 1.0e-8,
                "remap audit must include moving air, not just the head coordinates");
        const auto checkVector = [] (const auto& old, const auto& next, const char* name)
        {
            double error = 0.0, norm = 0.0;
            for (std::size_t coordinate = 0; coordinate < old.size(); ++coordinate)
            {
                error += (old[coordinate] - next[coordinate]) * (old[coordinate] - next[coordinate]);
                norm += old[coordinate] * old[coordinate];
            }
            expect (norm > 0.0 && std::sqrt (error / norm) < 4.0e-6,
                    std::string ("structural automation must preserve full physical ") + name);
        };
        checkVector (before.position, after.position, "displacement");
        checkVector (before.velocity, after.velocity, "velocity");
        checkVector (before.pending, after.pending, "pending force displacement");
    }
};
} // namespace taikor

namespace
{
std::vector<float> renderFamily (taikor::TaikoEngine& engine, int family,
                                 double rate, int blockSize)
{
    engine.reset();
    const int length = static_cast<int> (rate * 0.06);
    const int secondHit = length / 3;
    std::vector<float> left (static_cast<std::size_t> (length));
    std::vector<float> right (static_cast<std::size_t> (length));
    engine.trigger (taikor::Articulation::Don, family, 0.83f);
    for (int sample = 0; sample < length;)
    {
        if (sample == secondHit)
            engine.trigger (taikor::Articulation::DonRim, family, 0.61f);
        const int nextEvent = sample < secondHit ? secondHit : length;
        const int count = std::min (blockSize, nextEvent - sample);
        engine.process (left.data() + sample, right.data() + sample, count);
        sample += count;
    }
    left.insert (left.end(), right.begin(), right.end());
    return left;
}

void checkRenderAndPartition()
{
    for (double rate : { 8000.0, 48000.0, 192000.0 })
    {
        auto engine = std::make_unique<taikor::TaikoEngine>();
        taikor::EngineParameters parameters;
        parameters.humanise = 0.37f;
        engine->setParameters (parameters);
        engine->prepare (rate, 257);
        for (int family = 0; family < 4; ++family)
        {
            const std::string where = " (family " + std::to_string (family)
                                    + ", rate " + std::to_string (rate) + ")";
            const auto reference = renderFamily (*engine, family, rate, 64);
            const auto partitioned = renderFamily (*engine, family, rate, 257);
            expect (reference == partitioned,
                    "shared physical drums must replay bit-for-bit across buffer partitions" + where);
            double energy = 0.0;
            bool finite = true;
            for (float sample : reference)
            {
                finite = finite && std::isfinite (sample) && std::abs (sample) <= 0.891252f;
                energy += static_cast<double> (sample) * sample;
            }
            expect (finite && energy > 1.0e-9,
                    "every shared-cavity family/rate must render finite, protected, nonzero audio" + where);
        }
    }
}
} // namespace

int main()
{
    expect (taikor::TaikoEngine::hasRealismFeature (taikor::TaikoEngine::sharedCavity),
            "shipping defaults must enable the integrated shared cavity");
    for (int family = 0; family < 4; ++family)
    {
        taikor::TaikoEngineTestAccess::checkIntegratedBasis (family, 0.85f);
        taikor::TaikoEngineTestAccess::checkIntegratedBasis (family, 0.0f);
    }
    taikor::TaikoEngineTestAccess::checkPhysicalRetuning();
    taikor::TaikoEngineTestAccess::checkStructuralRemap();
    taikor::TaikoEngineTestAccess::checkReadoutAndPalmPersistence();
    checkRenderAndPartition();
    if (failures == 0)
        std::cout << "Shared cavity integration tests passed\n";
    return failures == 0 ? 0 : 1;
}
