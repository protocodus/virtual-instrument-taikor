#include "DSP/ShellBoundary.h"
#include "DSP/TaikoEngine.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>

namespace taikor
{
struct TaikoEngineTestAccess
{
    static bool verify (bool rear)
    {
        EngineParameters before;
        before.humanise = before.tensionModulation = 0.0f;
        before.pitch = 12.0f;
        constexpr int family = 3;
        const auto oldDrum = TaikoEngine::resolveDrumFor (before, 0.0f, family);
        auto probe = std::make_unique<TaikoEngine>();
        probe->setParameters (before);
        probe->prepare (384000.0, 64);
        probe->ensurePhysicalDrum (family, oldDrum);
        const auto modeAt = [] (TaikoEngine::Voice& bank, int id) -> TaikoEngine::Mode*
        {
            for (int index = 0; index < bank.modeCount; ++index)
                if (bank.modes[static_cast<std::size_t> (index)].physicalIndex == id)
                    return &bank.modes[static_cast<std::size_t> (index)];
            return nullptr;
        };
        auto& full = probe->physicalDrums_[family];
        int firstId = -1;
        double sampleRate = 0.0;
        for (int entry = 0; entry < TaikoEngine::modeEntryCount; ++entry)
        {
            if (TaikoEngine::membraneModes()[static_cast<std::size_t> (entry)].circumferentialOrder == 0)
                continue;
            const int id = (rear ? TaikoEngine::rearMembraneOffset : 0) + 2 * entry;
            const auto* high = modeAt (full, id);
            const auto* low = modeAt (full, id + 1);
            if (high == nullptr || low == nullptr)
                continue;
            const double rate = (high->liveOmega + low->liveOmega)
                              / (2.0 * 6.283185307179586 * 0.49);
            if (rate >= 8000.0 && rate <= 192000.0)
            {
                firstId = id;
                sampleRate = rate;
                break;
            }
        }
        if (firstId < 0)
        {
            std::cerr << "No eligible split pair\n";
            return false;
        }

        auto engine = std::make_unique<TaikoEngine>();
        engine->setParameters (before);
        // Put the real Nyquist cutoff strictly between the two split poles.
        engine->prepare (sampleRate, 64);
        engine->ensurePhysicalDrum (family, oldDrum);
        auto& physical = engine->physicalDrums_[family];
        auto* surviving = modeAt (physical, firstId + 1);
        if (modeAt (physical, firstId) != nullptr || surviving == nullptr)
        {
            std::cerr << "Initial pair does not straddle cutoff at " << sampleRate << " Hz\n";
            return false;
        }
        const auto oldBasis = surviving->angularBasis;
        constexpr double q = 1.0e-4, v = 0.2, pending = 2.0e-7;
        surviving->resonator.y1 = q;
        surviving->resonator.y2 = ((v + surviving->decayRate * q) / surviving->liveOmega
            - surviving->quadratureFromCurrent * q) / surviving->quadratureFromPrevious;
        physical.modalInput[static_cast<std::size_t> (firstId + 1)] =
            static_cast<float> (pending / surviving->resonator.b0);
        surviving->localMuteDampingRate = 0.7f;

        auto after = before;
        after.pitch = 0.0f;
        after.octaveBody = 0.0f; // switches to the reference head's nodal axes
        engine->setParameters (after);
        const auto newDrum = TaikoEngine::resolveDrumFor (after, 0.0f, family);
        engine->ensurePhysicalDrum (family, newDrum);
        const auto close = [] (double a, double b)
        { return std::abs (a - b) <= 4.0e-6 * std::max ({ std::abs (a), std::abs (b), 1.0e-12 }); };
        std::array<double, 2> rebuiltQ {}, rebuiltV {}, rebuiltPending {};
        for (int branch = 0; branch < 2; ++branch)
        {
            const auto* mode = modeAt (physical, firstId + branch);
            if (mode == nullptr)
            {
                std::cerr << "Partner still above cutoff after structural change\n";
                return false;
            }
            const double projection = mode->angularBasis[0] * oldBasis[0]
                                    + mode->angularBasis[1] * oldBasis[1];
            if (branch == 0 && std::abs (projection) < 0.01)
            {
                std::cerr << "Insufficient basis rotation: " << projection << '\n';
                return false; // the newly admitted mode must carry real motion
            }
            const double actualV = shellboundary::modalVelocity (
                mode->liveOmega, mode->decayRate, mode->quadratureFromCurrent,
                mode->quadratureFromPrevious, mode->resonator.y1, mode->resonator.y2);
            const double actualPending = mode->resonator.b0
                * physical.modalInput[static_cast<std::size_t> (firstId + branch)];
            if (! close (mode->resonator.y1, q * projection)
                || ! close (actualV, v * projection)
                || ! close (actualPending, pending * projection)
                || ! close (mode->localMuteDampingRate, 0.7))
            {
                std::cerr << "Projection mismatch rear=" << rear << " branch=" << branch
                          << " q=" << mode->resonator.y1 << "/" << q * projection
                          << " v=" << actualV << "/" << v * projection
                          << " pending=" << actualPending << "/" << pending * projection
                          << " mute=" << mode->localMuteDampingRate << '\n';
                return false;
            }
            for (std::size_t axis = 0; axis < 2; ++axis)
            {
                rebuiltQ[axis] += mode->angularBasis[axis] * mode->resonator.y1;
                rebuiltV[axis] += mode->angularBasis[axis] * actualV;
                rebuiltPending[axis] += mode->angularBasis[axis] * actualPending;
            }
        }
        for (std::size_t axis = 0; axis < 2; ++axis)
            if (! close (rebuiltQ[axis], oldBasis[axis] * q)
                || ! close (rebuiltV[axis], oldBasis[axis] * v)
                || ! close (rebuiltPending[axis], oldBasis[axis] * pending))
            {
                std::cerr << "Physical reconstruction mismatch at axis " << axis << '\n';
                return false;
            }
        return true;
    }
};
}

int main()
{
    if (! taikor::TaikoEngineTestAccess::verify (false)
        || ! taikor::TaikoEngineTestAccess::verify (true))
    {
        std::cerr << "FAIL: a newly in-band angular partner lost its physical projection\n";
        return 1;
    }
    std::cout << "Front/rear angular Nyquist-crossing remaps preserve physical motion.\n";
}
