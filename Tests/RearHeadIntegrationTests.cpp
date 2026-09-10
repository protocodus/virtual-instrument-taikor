#include "DSP/EnsembleEngine.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace
{
int failures = 0;
void expect (bool condition, const char* message)
{
    if (! condition)
    {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}
}

namespace taikor
{
struct TaikoEngineTestAccess
{
    static void verifyBank()
    {
        auto engine = std::make_unique<TaikoEngine>();
        engine->prepare (48000.0, 256);
        auto front = std::make_unique<TaikoEngine::Voice>();
        auto rear = std::make_unique<TaikoEngine::Voice>();
        const auto& stroke = TaikoEngine::strikeProfile (Articulation::Don);
        for (int family = 0; family < drumCount; ++family)
        {
            const auto drum = TaikoEngine::resolveDrumFor (EngineParameters {}, 0.0f, family);
            front->strikeRadius = rear->strikeRadius = 0.4f;
            front->strikeAngle = rear->strikeAngle = 0.3f;
            rear->rearStrike = true;
            engine->buildVoiceModes (*front, drum, stroke, 0.0f, true);
            engine->buildVoiceModes (*rear, drum, stroke, 0.0f, true);
            expect (front->modeCount == rear->modeCount, "head selection must preserve the physical bank");
            const auto find = [] (const TaikoEngine::Voice& bank, int id) -> const TaikoEngine::Mode*
            {
                for (int i = 0; i < bank.modeCount; ++i)
                    if (bank.modes[static_cast<std::size_t> (i)].physicalIndex == id)
                        return &bank.modes[static_cast<std::size_t> (i)];
                return nullptr;
            };
            int rearCount = 0;
            for (int i = 0; i < rear->modeCount; ++i)
            {
                auto& mode = rear->modes[static_cast<std::size_t> (i)];
                const auto* other = find (*front, mode.physicalIndex);
                expect (other != nullptr, "both strike sides need the same unique physical coordinates");
                if (other == nullptr)
                    continue;
                expect (mode.omega == other->omega && mode.micLeft == other->micLeft
                        && mode.micLeftQuadrature == other->micLeftQuadrature,
                        "switching the struck head must not move a pole or microphone");
                if (! mode.membrane || mode.circumferentialOrder == 0)
                    continue;
                if (mode.rearHeadMode)
                {
                    ++rearCount;
                    expect (other->drive == 0.0f && other->contactShape == 0.0f,
                            "a front stroke must not directly force the rear angular bank");
                    if (mode.modeEntry < TaikoEngine::legacyModeEntryCount)
                        expect (std::abs (mode.contactShape) > 1.0e-9f,
                                "a rear stroke needs reciprocal non-axisymmetric contact");
                    expect (mode.batterTensionFraction == 0.0f && mode.rearTensionFraction > 0.0f,
                            "rear restoring energy must belong to the rear skin");
                }
                else
                    expect (mode.drive == 0.0f && mode.contactShape == 0.0f,
                            "a rear stroke must not directly force the batter angular bank");
            }
            expect (rearCount > 24, "rear-head strikes need an extended angular bank");

            // One displaced rear angular mode stretches only that skin. A
            // global same-head label would incorrectly harden the batter too.
            rear->activeModeCount = rear->modeCount;
            for (int i = 0; i < rear->modeCount; ++i)
            {
                auto& mode = rear->modes[static_cast<std::size_t> (i)];
                mode.resonator.clear();
                if (mode.rearHeadMode && mode.modeEntry == 4)
                    mode.resonator.y1 = 0.002;
            }
            expect (TaikoEngine::membraneSquaredSlope (*rear) == 0.0,
                    "rear angular displacement must not appear as batter strain");
            expect (TaikoEngine::membraneSquaredSlope (*rear, true) > 0.0,
                    "rear angular displacement must produce rear strain");
        }
    }
    static bool allContactsRear (const TaikoEngine& engine)
    {
        bool found = false;
        for (const auto& voice : engine.voices_)
            if (voice.active)
            {
                found = true;
                if (! voice.rearStrike)
                    return false;
            }
        return found;
    }
    static void verifyContactAutomation()
    {
        auto engine = std::make_unique<TaikoEngine>();
        EngineParameters params;
        params.humanise = 0.0f;
        engine->setParameters (params);
        engine->prepare (48000.0, 64);
        engine->setRearHeadStrike (true);
        engine->trigger (Articulation::Don, 2, 0.8f);
        float left {}, right {};
        engine->process (&left, &right, 1);
        auto& contact = engine->voices_[0];
        expect (contact.nonlinearContactActive && contact.rearStrike,
                "automation probe must begin during a rear-head collision");
        params.cavityCoupling = 0.4f;
        params.headMaterial = 0.65f;
        engine->setParameters (params);
        engine->refreshDrumIfNeeded();
        const auto& drum = engine->drumCache_[2];
        engine->ensurePhysicalDrum (2, drum);
        auto expected = std::make_unique<TaikoEngine::Voice>();
        expected->strikeRadius = contact.strikeRadius;
        expected->strikeAngle = contact.strikeAngle;
        expected->rearStrike = true;
        engine->buildVoiceModes (*expected, drum,
            TaikoEngine::strikeProfile (Articulation::Don), 0.0f, false);
        for (int i = 0; i < expected->modeCount; ++i)
        {
            const auto& mode = expected->modes[static_cast<std::size_t> (i)];
            expect (contact.contactProjection[mode.physicalIndex] == mode.contactShape
                    && contact.modeProjection[mode.physicalIndex] == mode.drive,
                    "structural automation must retain the latched head on both contact paths");
        }
    }
};
struct EnsembleEngineTestAccess
{
    static void verifyScheduledSide()
    {
        auto ensemble = std::make_unique<EnsembleEngine>();
        EngineParameters params;
        params.ensembleSize = 4;
        params.ensembleVariation = 1.0f;
        ensemble->setParameters (params);
        ensemble->prepare (48000.0, 256);
        ensemble->setRearHeadStrike (true);
        ensemble->trigger (Articulation::Don, 2, 0.6f);
        expect (ensemble->pendingCount > 0, "test must schedule delayed companions");
        ensemble->setRearHeadStrike (false);
        std::array<float, 256> left {}, right {};
        for (int i = 0; i < 7; ++i)
            ensemble->process (left.data(), right.data(), 256);
        for (int member = 0; member < 4; ++member)
            expect (TaikoEngineTestAccess::allContactsRear (*ensemble->players[member]),
                    "later CC18 changes must not change an already scheduled stroke's side");
        ensemble->reset();
        expect (! ensemble->rearHeadStrike, "reset must restore front-head strokes");
    }
};
}

namespace
{
std::vector<float> performance (double rate, int block, bool rear)
{
    auto engine = std::make_unique<taikor::TaikoEngine>();
    taikor::EngineParameters params;
    params.humanise = 0.0f;
    engine->setParameters (params);
    engine->prepare (rate, block);
    std::vector<float> result;
    for (int family = 0; family < taikor::drumCount; ++family)
    {
        engine->setRearHeadStrike (rear);
        engine->trigger (taikor::Articulation::Don, family, 0.75f);
        const int count = static_cast<int> (0.12 * rate);
        std::vector<float> left (count), right (count);
        for (int offset = 0; offset < count; offset += block)
            engine->process (left.data() + offset, right.data() + offset,
                             std::min (block, count - offset));
        for (int sample = 0; sample < count; ++sample)
        {
            expect (std::isfinite (left[sample]) && std::isfinite (right[sample]),
                    "both head paths must render finite samples");
            result.push_back (left[sample]);
            result.push_back (right[sample]);
        }
    }
    return result;
}
}

int main()
{
    taikor::TaikoEngineTestAccess::verifyBank();
    taikor::TaikoEngineTestAccess::verifyContactAutomation();
    taikor::EnsembleEngineTestAccess::verifyScheduledSide();
    for (const double rate : { 8000.0, 48000.0, 192000.0 })
    {
        const auto rear = performance (rate, 64, true);
        expect (rear == performance (rate, 257, true), "rear performance must be block-partition deterministic");
        const auto front = performance (rate, 64, false);
        expect (front != rear, "independently constructed skins must produce a front/back contrast");
        double energy = 0.0;
        for (const auto sample : rear)
            energy += static_cast<double> (sample) * sample;
        expect (energy > 1.0e-6, "rear strokes must reach the microphones");
    }
    if (failures != 0)
        return 1;
    std::cout << "Rear-head contact, independent strain and scheduled-side replay passed.\n";
}
