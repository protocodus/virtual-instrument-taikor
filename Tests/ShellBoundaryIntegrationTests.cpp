#include "DSP/ShellBoundary.h"
#include "DSP/TaikoEngine.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>

namespace taikor
{
struct ShellBoundaryIntegrationTestAccess
{
    static double modalEnergy (const TaikoEngine::Voice& voice)
    {
        double result = 0.0;
        for (int index = 0; index < voice.activeModeCount; ++index)
        {
            const auto& mode = voice.modes[static_cast<std::size_t> (index)];
            if (! (mode.inverseModalMass > 0.0f))
                continue;
            const double decay = mode.decayRate + mode.appliedPalmDecay;
            const double velocity = shellboundary::modalVelocity (
                mode.liveOmega, decay, mode.quadratureFromCurrent, mode.quadratureFromPrevious,
                mode.resonator.y1, mode.resonator.y2);
            const double restoring = mode.liveOmega * mode.liveOmega + decay * decay;
            result += (velocity * velocity + restoring * mode.resonator.y1 * mode.resonator.y1)
                      / (2.0 * mode.inverseModalMass);
        }
        return result;
    }

    static bool checkRuntime()
    {
        // Use the real builder and process hook. The synthetic-port check below
        // cannot detect a missing builder projection or an uncalled render step.
        for (int octave = 0; octave < drumCount; ++octave)
        {
            TaikoEngine engine;
            EngineParameters parameters;
            parameters.tensionModulation = 0.0f;
            parameters.strikeNoise = 0.0f;
            parameters.humanise = 0.0f;
            engine.setParameters (parameters);
            engine.prepare (48000.0, 64);
            engine.trigger (Articulation::Ka, octave, 0.7f);
            auto& physical = engine.physicalDrums_[static_cast<std::size_t> (octave)];
            int headPorts = 0, shellPorts = 0;
            for (int index = 0; index < physical.modeCount; ++index)
            {
                const auto& mode = physical.modes[static_cast<std::size_t> (index)];
                if (mode.shellBoundaryImpulseScale != 0.0)
                    (mode.membrane ? headPorts : shellPorts)++;
                if (! mode.membrane)
                    for (const auto& strike : engine.voices_)
                        if (strike.active && strike.modeProjection[mode.physicalIndex] != 0.0f)
                            return false; // an ordinary head stroke cannot duplicate its force
            }
            if (headPorts == 0 || shellPorts == 0)
                return false;
            std::array<float, 64> left {}, right {};
            for (int block = 0; block < 16; ++block)
                engine.process (left.data(), right.data(), 64);
            double woodState = 0.0;
            for (const auto& mode : physical.modes)
                if (! mode.membrane)
                    woodState += mode.resonator.y1 * mode.resonator.y1
                                 + mode.resonator.y2 * mode.resonator.y2;
            if (! std::isfinite (woodState) || woodState <= 1.0e-30)
                return false; // the shell must actually receive boundary energy

            // Reverse experiment: start the wall moving, stop all external
            // contacts, and observe energy enter the head. The exact free-pole
            // and boundary composition must stay passive without forcing.
            for (auto& strike : engine.voices_)
                strike.active = false;
            bool started = false;
            for (int index = 0; index < physical.modeCount; ++index)
            {
                auto& mode = physical.modes[static_cast<std::size_t> (index)];
                mode.resonator.clear();
                if (! started && ! mode.membrane && mode.shellBoundaryImpulseScale != 0.0)
                {
                    mode.resonator.y2 = 0.01 * mode.shellBoundaryVelocityToPrevious;
                    started = true;
                }
            }
            physical.modalInput.fill (0.0f);
            double previousEnergy = modalEnergy (physical);
            bool headMoved = false;
            for (int sample = 0; sample < 1024; ++sample)
            {
                engine.process (left.data(), right.data(), 1);
                const double energy = modalEnergy (physical);
                if (! std::isfinite (energy) || energy > previousEnergy + 1.0e-10)
                    return false;
                previousEnergy = energy;
                if (sample == 16)
                    for (const auto& mode : physical.modes)
                        if (mode.membrane && std::abs (mode.resonator.y1) > 1.0e-20)
                            headMoved = true;
            }
            if (! headMoved)
                return false;
        }
        return true;
    }

    static bool checkDuringContact()
    {
        // Audit the boundary step on states reached during actual front/rear
        // contact, including direct rim excitation. Restore the probe's pole
        // history afterward so the normal solver still advances exactly once.
        // This proves the boundary's passivity; the contact solver has its own
        // stricter stick/spring/reciprocal-mode energy audit in TaikoEngineTests.
        for (double rate : { 8000.0, 48000.0, 192000.0 })
        for (int octave = 0; octave < drumCount; ++octave)
        for (bool rear : { false, true })
        for (auto articulation : { Articulation::Ka, Articulation::DonRim })
        {
            TaikoEngine engine;
            EngineParameters parameters;
            parameters.tensionModulation = parameters.strikeNoise = parameters.humanise = 0.0f;
            engine.setParameters (parameters);
            engine.prepare (rate, 1);
            engine.setRearHeadStrike (rear);
            engine.trigger (articulation, octave, 0.9f);
            bool audited = false;
            for (int sample = 0; sample < static_cast<int> (rate * 0.03); ++sample)
            {
                float left {}, right {};
                engine.process (&left, &right, 1);
                bool hasForce = false;
                for (const auto& strike : engine.voices_)
                    hasForce = hasForce || strike.solvedContactForce > 0.0;
                if (! hasForce)
                    continue;
                auto& physical = engine.physicalDrums_[static_cast<std::size_t> (octave)];
                std::array<double, TaikoEngine::resonatorCount> previous {};
                for (int index = 0; index < physical.activeModeCount; ++index)
                    previous[static_cast<std::size_t> (index)] = physical.modes[static_cast<std::size_t> (index)].resonator.y2;
                const double before = modalEnergy (physical);
                TaikoEngine::applyShellBoundary (physical);
                const double after = modalEnergy (physical);
                if (! std::isfinite (after) || after > before + 1.0e-11 * std::max (before, 1.0))
                    return false;
                for (int index = 0; index < physical.activeModeCount; ++index)
                    physical.modes[static_cast<std::size_t> (index)].resonator.y2 = previous[static_cast<std::size_t> (index)];
                audited = true;
            }
            if (! audited)
                return false;
        }
        return true;
    }

    static bool checkLiveObservation()
    {
        auto changed = std::make_unique<TaikoEngine>();
        auto reference = std::make_unique<TaikoEngine>();
        EngineParameters parameters;
        parameters.tensionModulation = parameters.humanise = 0.0f;
        parameters.shellResonance = 0.0f;
        std::array<float, 64> left {}, right {}, referenceLeft {}, referenceRight {};
        for (auto* engine : { changed.get(), reference.get() })
        {
            engine->setParameters (parameters);
            engine->prepare (48000.0, 64);
            engine->trigger (Articulation::DonRim, 2, 0.9f);
            engine->process (left.data(), right.data(), 64);
        }
        const auto revision = changed->physicalConfigurationRevision_;
        parameters.shellResonance = 1.0f;
        changed->setParameters (parameters);
        if (changed->physicalConfigurationRevision_ != revision)
            return false;
        double audibleDifference = 0.0;
        for (int block = 0; block < 8; ++block)
        {
            changed->process (left.data(), right.data(), 64);
            reference->process (referenceLeft.data(), referenceRight.data(), 64);
            for (std::size_t sample = 0; sample < left.size(); ++sample)
                audibleDifference = std::max (audibleDifference,
                    std::abs (static_cast<double> (left[sample] - referenceLeft[sample])));
        }
        const auto& a = changed->physicalDrums_[2];
        const auto& b = reference->physicalDrums_[2];
        for (int index = 0; index < a.modeCount; ++index)
        {
            const auto& actual = a.modes[static_cast<std::size_t> (index)];
            const auto& neutral = b.modes[static_cast<std::size_t> (index)];
            if (actual.resonator.y1 != neutral.resonator.y1
                || actual.resonator.y2 != neutral.resonator.y2
                || actual.resonator.a1 != neutral.resonator.a1
                || actual.resonator.a2 != neutral.resonator.a2)
                return false;
        }
        return audibleDifference > 1.0e-4;
    }

    static bool check()
    {
        TaikoEngine engine;
        bool valid = true;
        const auto near = [] (double a, double b)
        { return std::abs (a - b) <= 2.0e-8 * std::max ({ 1.0, std::abs (a), std::abs (b) }); };
        for (double rate : { 8000.0, 48000.0, 192000.0 })
        {
            engine.prepare (rate, 64);
            auto storage = std::make_unique<TaikoEngine::Voice>();
            auto& voice = *storage;
            voice.physicalBank = true;
            voice.modeCount = voice.activeModeCount = 2;
            voice.shellBoundaryDamping[0] = 30.0;
            for (std::size_t i = 0; i < 2; ++i)
            {
                auto& mode = voice.modes[i];
                const float frequency = i == 0 ? 190.0f : 74.0f;
                TaikoEngine::beginMode (mode, 6.2831853f * frequency, frequency, 2.5f);
                mode.inverseModalMass = i == 0 ? 2.0f : 0.125f;
                mode.membrane = i == 0;
                mode.shellBoundaryGroup = 0;
                mode.shellBoundaryProjection = i == 0 ? 0.3f : -1.0f;
                engine.configureResonator (mode.resonator, frequency, mode.decayRate,
                                            1.0f, &mode.poleRadius);
                TaikoEngine::updateQuadratureScales (mode);
            }
            engine.configureShellBoundary (voice);
            const auto velocity = [] (const TaikoEngine::Mode& mode)
            {
                return shellboundary::modalVelocity (
                    mode.liveOmega, static_cast<double> (mode.decayRate + mode.appliedPalmDecay),
                    mode.quadratureFromCurrent, mode.quadratureFromPrevious,
                    mode.resonator.y1, mode.resonator.y2);
            };
            for (std::size_t movingSide = 0; movingSide < 2; ++movingSide)
            {
                for (std::size_t i = 0; i < 2; ++i)
                {
                    auto& mode = voice.modes[i];
                    mode.resonator.y1 = 0.001 * static_cast<double> (i + 1);
                    const double requestedVelocity = i == movingSide ? 0.9 : 0.0;
                    mode.resonator.y2 = ((requestedVelocity + mode.decayRate * mode.resonator.y1)
                                         / mode.liveOmega - mode.quadratureFromCurrent * mode.resonator.y1)
                                        / mode.quadratureFromPrevious;
                }
                const std::array<double, 2> before { velocity (voice.modes[0]), velocity (voice.modes[1]) };
                const double relative = static_cast<double> (voice.modes[0].shellBoundaryProjection) * before[0]
                                        - before[1];
                const double g = std::pow (static_cast<double> (voice.modes[0].shellBoundaryProjection), 2.0)
                                 * voice.modes[0].inverseModalMass + voice.modes[1].inverseModalMass;
                const auto damper = shellboundary::makeDashpot (g, 30.0, 1.0 / rate);
                const double impulse = damper.impulsePerVelocity * relative;
                TaikoEngine::applyShellBoundary (voice);
                double energyBefore = 0.0, energyAfter = 0.0;
                for (std::size_t i = 0; i < 2; ++i)
                {
                    const auto& mode = voice.modes[i];
                    const double actual = velocity (mode);
                    const double expected = before[i] + mode.inverseModalMass * mode.shellBoundaryProjection * impulse;
                    valid = valid && near (actual, expected);
                    valid = valid && mode.resonator.y1 == 0.001 * static_cast<double> (i + 1);
                    energyBefore += before[i] * before[i] / (2.0 * mode.inverseModalMass);
                    energyAfter += actual * actual / (2.0 * mode.inverseModalMass);
                }
                valid = valid && energyAfter <= energyBefore + 1.0e-10;
                valid = valid && near (energyBefore - energyAfter,
                                        damper.energyLossPerVelocitySquared * relative * relative);
                valid = valid && std::abs (velocity (voice.modes[1 - movingSide])) > 1.0e-8;
            }
            // A filtered-out ring must not turn the same port into a damper
            // attached to an immovable ground that is absent from the model.
            voice.modes[1].shellBoundaryGroup = -1;
            engine.configureShellBoundary (voice);
            valid = valid && voice.modes[0].shellBoundaryImpulseScale == 0.0;
            valid = valid && voice.modes[1].shellBoundaryImpulseScale == 0.0;
        }
        return valid;
    }
};
} // namespace taikor

int main()
{
    if (! taikor::ShellBoundaryIntegrationTestAccess::check()
        || ! taikor::ShellBoundaryIntegrationTestAccess::checkRuntime()
        || ! taikor::ShellBoundaryIntegrationTestAccess::checkDuringContact()
        || ! taikor::ShellBoundaryIntegrationTestAccess::checkLiveObservation())
    {
        std::cerr << "FAIL: canonical shell/head modal coupling is not passive and reciprocal\n";
        return 1;
    }
    std::cout << "Canonical resonator shell/head coupling tests passed.\n";
}
