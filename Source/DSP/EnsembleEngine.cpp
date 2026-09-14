#include "EnsembleEngine.h"

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

namespace taikor
{
// Offline only. The calling thread participates in the work and waits for all
// helpers before ordered mixing or any parameter/MIDI/lifecycle mutation.
// There are no detached jobs, and no worker ever touches another player's DSP.
struct EnsembleEngine::OfflineRenderPool
{
    using Job = void (*) (void*, int) noexcept;

    explicit OfflineRenderPool (int count)
    {
        try
        {
            workers.reserve (static_cast<std::size_t> (count));
            for (int index = 0; index < count; ++index)
                workers.emplace_back ([this] { work(); });
        }
        catch (...)
        {
            stop();
            throw;
        }
    }

    ~OfflineRenderPool() { stop(); }

    void stop() noexcept
    {
        {
            const std::lock_guard lock (mutex);
            stopping = true;
        }
        available.notify_all();
        for (auto& worker : workers)
            if (worker.joinable())
                worker.join();
    }

    void run (void* nextContext, Job nextJob, int count) noexcept
    {
        {
            const std::lock_guard lock (mutex);
            context = nextContext;
            job = nextJob;
            jobCount = count;
            next.store (0, std::memory_order_relaxed);
            remaining = static_cast<int> (workers.size());
            // Inherit the caller's rounding/denormal mode, including a host's
            // scoped floating-point controls, rather than the prepare thread's.
            environmentValid = std::fegetenv (&environment) == 0;
            ++generation;
        }
        available.notify_all();
        execute();
        std::unique_lock lock (mutex);
        completed.wait (lock, [this] { return remaining == 0; });
    }

    void execute() noexcept
    {
        for (int index = next.fetch_add (1, std::memory_order_relaxed);
             index < jobCount; index = next.fetch_add (1, std::memory_order_relaxed))
            job (context, index);
    }

    void work() noexcept
    {
        std::uint64_t observed = 0;
        for (;;)
        {
            std::unique_lock lock (mutex);
            available.wait (lock, [this, observed]
                { return stopping || generation != observed; });
            if (stopping)
                return;
            observed = generation;
            lock.unlock();
            if (environmentValid)
                (void) std::fesetenv (&environment);
            execute();
            lock.lock();
            if (--remaining == 0)
                completed.notify_one();
        }
    }

    std::vector<std::thread> workers;
    std::mutex mutex;
    std::condition_variable available, completed;
    std::atomic<int> next { 0 };
    int jobCount = 0, remaining = 0;
    std::uint64_t generation = 0;
    void* context = nullptr;
    Job job = nullptr;
    std::fenv_t environment {};
    bool environmentValid = false, stopping = false;
};

EnsembleEngine::EnsembleEngine()
{
    for (int member = 0; member < maximumEnsembleSize; ++member)
    {
        players[member] = std::make_unique<TaikoEngine>();
        players[member]->setEnsembleMember (member);
    }
}

EnsembleEngine::~EnsembleEngine() = default;

void EnsembleEngine::prepareOfflineRendering (int workerCount) noexcept
{
    releaseOfflineRendering();
    const auto cores = std::thread::hardware_concurrency();
    const int availableWorkers = cores > 0 ? static_cast<int> (std::min (cores - 1u, 3u)) : 3;
    workerCount = std::clamp (workerCount, 0, availableWorkers);
    if (workerCount == 0)
        return;
    try
    {
        offlinePool = std::make_unique<OfflineRenderPool> (workerCount);
    }
    catch (...)
    {
        // Resource exhaustion must not prevent the instrument from rendering.
        offlinePool.reset();
    }
}

void EnsembleEngine::releaseOfflineRendering() noexcept
{
    offlineRendering = false;
    offlinePool.reset();
}

int EnsembleEngine::getOfflineWorkerCount() const noexcept
{
    return offlinePool ? static_cast<int> (offlinePool->workers.size()) : 0;
}

void EnsembleEngine::prepare (double sampleRate, int maxBlockSize) noexcept
{
    // Match the physical engine's supported clock before scheduling delays.
    rate = std::isfinite (sampleRate) ? std::clamp (sampleRate, 8000.0, 384000.0)
                                    : 48000.0;
    gainSmoothing = -std::expm1 (-1.0 / (0.015 * rate));
    for (auto& player : players)
        player->prepare (rate, maxBlockSize);
    prepared = true;
    reset();
}

void EnsembleEngine::reset() noexcept
{
    for (auto& player : players)
        player->reset();
    pendingCount = 0;
    sampleClock = strokeSequence = eventSequence = 0;
    gain = 1.0f / std::sqrt (static_cast<float> (parameters.ensembleSize));
    pan.fill ({});
    panTarget.fill ({});
    updateStage (true);
    clearStrikeOverrides();
    rearHeadStrike = false;
    publishVoices();
}

void EnsembleEngine::allSoundsOff() noexcept
{
    for (auto& player : players)
        player->allSoundsOff();
    pendingCount = 0;
    // Panic clears scheduled hits but preserves each player's gesture sequence.
    gain = 1.0f / std::sqrt (static_cast<float> (parameters.ensembleSize));
    updateStage (true);
    publishVoices();
}

EngineParameters EnsembleEngine::memberParameters (int member) const noexcept
{
    // The lead plays the drum the controls describe. Every companion plays a
    // physically distinct instrument: no two drums of a real ensemble share
    // one hide, one tension or one shell, and eight copies of one solve read
    // as one drum sampled eight times rather than as eight drums. Each
    // companion therefore carries a stable, per-member offset on the four
    // construction controls, scaled by Ensemble Variation so that zero still
    // means identical instruments. The ranges are expressive control ranges
    // in control units - a few per cent of hide density and shell, a tension
    // spread worth about a quarter-tone at full variation - not measured
    // ensemble statistics, and every member still keeps the family's own
    // geometry, so the four rows stay the four drums.
    auto result = parameters;
    if (member <= 0 || parameters.ensembleVariation <= 0.0f
        || ! TaikoEngine::hasRealismFeature (TaikoEngine::distinctEnsembleDrums))
        return result;
    const float spread = parameters.ensembleVariation;
    const auto unit = [member] (std::uint32_t salt) noexcept
    {
        const auto seed = hash (static_cast<std::uint32_t> (member) * 0x9e3779b9u
                                ^ salt);
        return 2.0f * static_cast<float> (seed) / 4294967295.0f - 1.0f;
    };
    result.tension = std::clamp (
        parameters.tension + 0.018f * spread * unit (0x1b873593u), 0.0f, 1.0f);
    result.headMaterial = std::clamp (
        parameters.headMaterial + 0.03f * spread * unit (0x85ebca6bu), 0.0f, 1.0f);
    result.shellMaterial = std::clamp (
        parameters.shellMaterial + 0.03f * spread * unit (0xc2b2ae35u), 0.0f, 1.0f);
    result.headDamping = std::clamp (
        parameters.headDamping + 0.04f * spread * unit (0x27d4eb2fu), 0.0f, 1.0f);
    return result;
}

void EnsembleEngine::setParameters (const EngineParameters& next) noexcept
{
    const int previousSize = parameters.ensembleSize;
    parameters = TaikoEngine::sanitiseParameters (next);
    for (int member = 0; member < maximumEnsembleSize; ++member)
        players[member]->setParameters (memberParameters (member));
    const bool silent = pendingCount == 0 && getActiveVoiceCount() == 0;
    if (silent)
        gain = 1.0f / std::sqrt (static_cast<float> (parameters.ensembleSize));
    if (parameters.ensembleSize != previousSize)
    {
        // A tail removed during an earlier move stays where it is now.
        for (int member = parameters.ensembleSize; member < previousSize; ++member)
            panTarget[member] = pan[member];
        updateStage (silent);
    }
}

float EnsembleEngine::stagePosition (int member, int size) noexcept
{
    if (size <= 1)
        return 0.0f;
    const int half = size / 2;
    // Odd ensembles include center; even ensembles fill both halves equally.
    const int slot = member - half + (size % 2 == 0 && member >= half ? 1 : 0);
    return static_cast<float> (slot) / static_cast<float> (half);
}

void EnsembleEngine::updateStage (bool snap) noexcept
{
    for (int member = 0; member < parameters.ensembleSize; ++member)
    {
        panTarget[member] = StereoPan::atPosition (
            stagePosition (member, parameters.ensembleSize));
        // A silent player can take its seat before the next attack. Only an
        // already-ringing microphone pair needs a gradual move.
        if (snap || players[member]->getActiveVoiceCount() == 0)
            pan[member] = panTarget[member];
    }
    // Removed members keep their stage positions as existing tails finish.
}

std::uint32_t EnsembleEngine::hash (std::uint32_t value) noexcept
{
    value ^= value >> 16;
    value *= 0x7feb352du;
    value ^= value >> 15;
    value *= 0x846ca68bu;
    return value ^ (value >> 16);
}

bool EnsembleEngine::later (const Hit& a, const Hit& b) noexcept
{
    return a.due > b.due || (a.due == b.due && a.order > b.order);
}

void EnsembleEngine::fire (const Hit& hit) noexcept
{
    auto& player = *players[hit.member];
    // Capture placement when MIDI arrives. Later controller changes cannot
    // move an already-scheduled player's hand to a different point.
    player.setStrikePositionOverride (hit.position);
    player.setStrikeAzimuthOverride (hit.azimuth);
    player.setRearHeadStrike (hit.rear);
    player.trigger (hit.articulation, hit.octave, hit.velocity,
                    hit.radial, hit.tangential);
    player.clearStrikeOverrides();
}

void EnsembleEngine::dispatchDue() noexcept
{
    while (pendingCount > 0 && pending[0].due <= sampleClock)
    {
        std::pop_heap (pending.begin(), pending.begin() + pendingCount, later);
        fire (pending[--pendingCount]);
    }
}

void EnsembleEngine::trigger (Articulation articulation, int octave, float velocity) noexcept
{
    if (! prepared || static_cast<std::size_t> (articulation) >= articulationCount
        || ! std::isfinite (velocity) || velocity <= 0.0f)
        return;
    dispatchDue();
    const auto sequence = ++strokeSequence;
    const float spread = parameters.ensembleVariation;
    for (int member = 0; member < parameters.ensembleSize; ++member)
    {
        Hit hit;
        hit.order = ++eventSequence;
        hit.due = sampleClock;
        hit.member = member;
        hit.octave = std::clamp (octave, lowestOctaveOffset, highestOctaveOffset);
        hit.articulation = articulation;
        hit.velocity = std::clamp (velocity, 0.0f, 1.0f);
        hit.position = positionOverridden ? positionOverride : parameters.strikePosition;
        hit.azimuth = azimuthOverridden ? azimuthOverride : parameters.strikeAzimuth;
        hit.rear = rearHeadStrike;
        if (member > 0)
        {
            const auto seed = hash (static_cast<std::uint32_t> (sequence))
                            ^ hash (static_cast<std::uint32_t> (sequence >> 32))
                            ^ hash (static_cast<std::uint32_t> (member) * 0x9e3779b9u)
                            ^ hash (static_cast<std::uint32_t> (parameters.performer));
            const auto unit = [] (std::uint32_t value)
            { return static_cast<double> (hash (value)) / 4294967295.0; };
            // The leader anchors the MIDI timestamp. Companions are late by
            // 0..30 ms, so live playing needs no added latency/lookahead.
            hit.due += static_cast<std::uint64_t> (std::llround (
                maximumEnsembleDelaySeconds * rate * spread * unit (seed)));
            hit.radial = static_cast<float> (0.055 * spread * (2.0 * unit (seed + 1u) - 1.0));
            hit.tangential = static_cast<float> (0.055 * spread * (2.0 * unit (seed + 2u) - 1.0));
        }
        if (hit.due == sampleClock)
            fire (hit);
        else if (pendingCount < queueCapacity)
        {
            pending[pendingCount++] = hit;
            std::push_heap (pending.begin(), pending.begin() + pendingCount, later);
        }
        // A pathological MIDI flood can fill the bounded delay queue. Drop
        // only excess delayed companions; the on-time leader always plays.
    }
    publishVoices();
}

bool EnsembleEngine::triggerMidi (int note, float velocity) noexcept
{
    const auto articulation = articulationForMidiNote (note);
    const auto octave = octaveOffsetForMidiNote (note);
    if (! prepared || ! articulation || ! octave || ! std::isfinite (velocity)
        || velocity <= 0.0f)
        return false;
    trigger (*articulation, *octave, velocity);
    return true;
}

void EnsembleEngine::process (float* left, float* right, int samples) noexcept
{
    if (left == nullptr || right == nullptr || samples <= 0)
        return;
    if (! prepared)
    {
        std::fill_n (left, samples, 0.0f);
        std::fill_n (right, samples, 0.0f);
        return;
    }
    int rendered = 0;
    while (rendered < samples)
    {
        dispatchDue();
        int count = std::min (blockCapacity, samples - rendered);
        if (pendingCount > 0)
            count = static_cast<int> (std::min<std::uint64_t> (
                static_cast<std::uint64_t> (count), pending[0].due - sampleClock));
        std::fill_n (extraLeft.data(), count, 0.0f);
        std::fill_n (extraRight.data(), count, 0.0f);
        bool extraActive = false;
        bool extraSilent = true;
        const auto silent = [] (float value) { return value == 0.0f; };
        std::array<bool, maximumEnsembleSize> wasActive {};
        bool parallel = false;
        if (offlineRendering && offlinePool && count >= 64)
        {
            int activeCompanions = 0;
            for (int member = 1; member < maximumEnsembleSize; ++member)
            {
                wasActive[member] = players[member]->getActiveVoiceCount() > 0;
                activeCompanions += wasActive[member] ? 1 : 0;
            }
            parallel = activeCompanions >= 2;
            if (parallel)
            {
                struct Batch { EnsembleEngine* engine; int samples; } batch { this, count };
                offlinePool->run (&batch, [] (void* context, int index) noexcept
                {
                    const auto& work = *static_cast<Batch*> (context);
                    auto& engine = *work.engine;
                    const auto member = static_cast<std::size_t> (index + 1);
                    engine.players[member]->processRaw (engine.memberLeft[member].data(),
                        engine.memberRight[member].data(), work.samples);
                }, maximumEnsembleSize - 1);
            }
        }
        for (int member = 1; member < maximumEnsembleSize; ++member)
        {
            auto& player = *players[member];
            const bool memberActive = parallel ? wasActive[member] : player.getActiveVoiceCount() > 0;
            extraActive = extraActive || memberActive;
            // Even silent players advance held-palm/pitch smoothers. Removing
            // a player from the size control lets their existing tail finish.
            auto* rawLeft = parallel ? memberLeft[member].data() : scratchLeft.data();
            auto* rawRight = parallel ? memberRight[member].data() : scratchRight.data();
            if (! parallel)
                player.processRaw (rawLeft, rawRight, count);
            // The published count covers ringing drums, so also check the raw
            // samples before skipping a possible remaining contact transient.
            if (! memberActive
                && std::all_of (rawLeft, rawLeft + count, silent)
                && std::all_of (rawRight, rawRight + count, silent))
            {
                // A silent member may still be moving after a size change.
                // Keep that recurrence exact; a settled matrix needs no work.
                const auto& current = pan[member];
                const auto& target = panTarget[member];
                if (current.ll != target.ll || current.lr != target.lr
                    || current.rl != target.rl || current.rr != target.rr)
                    for (int sample = 0; sample < count; ++sample)
                        pan[member].approach (target, gainSmoothing);
                continue;
            }
            extraSilent = false;
            for (int sample = 0; sample < count; ++sample)
            {
                pan[member].approach (panTarget[member], gainSmoothing);
                pan[member].apply (rawLeft[sample], rawRight[sample]);
                extraLeft[sample] += rawLeft[sample];
                extraRight[sample] += rawRight[sample];
            }
        }
        const float target = 1.0f / std::sqrt (static_cast<float> (parameters.ensembleSize));
        const bool unityGain = gain == 1.0f && target == 1.0f;
        const bool centeredLead = pan[0].isCentered() && panTarget[0].isCentered();
        if (! extraActive && unityGain && centeredLead
            && (extraSilent
                || (std::all_of (extraLeft.begin(), extraLeft.begin() + count, silent)
                    && std::all_of (extraRight.begin(), extraRight.begin() + count, silent))))
            players[0]->process (left + rendered, right + rendered, count);
        else
        {
            for (int sample = 0; sample < count; ++sample)
            {
                pan[0].approach (panTarget[0], gainSmoothing);
                leadPan[sample] = pan[0];
                gain += gainSmoothing * (target - gain);
                if (std::abs (gain - target) < 1.0e-10)
                    gain = target;
                mixGain[sample] = static_cast<float> (gain);
            }
            players[0]->processWithEnsemble (left + rendered, right + rendered, count,
                extraSilent ? nullptr : extraLeft.data(),
                extraSilent ? nullptr : extraRight.data(), mixGain.data(), extraActive,
                centeredLead ? nullptr : leadPan.data());
        }
        rendered += count;
        sampleClock += static_cast<std::uint64_t> (count);
    }
    publishVoices();
}

bool EnsembleEngine::isOutputFrozen() const noexcept
{
    return pendingCount == 0 && std::all_of (players.begin(), players.end(),
        [] (const auto& player) { return player->isOutputFrozen(); });
}

void EnsembleEngine::publishVoices() noexcept
{
    int count = 0;
    for (const auto& player : players)
        count += player->getActiveVoiceCount();
    activeVoices.store (count, std::memory_order_relaxed);
}

void EnsembleEngine::setHandDamping (float amount) noexcept
{ for (auto& player : players) player->setHandDamping (amount); }
void EnsembleEngine::setPitchBend (float amount) noexcept
{ for (auto& player : players) player->setPitchBend (amount); }
void EnsembleEngine::setStrikePositionOverride (float amount) noexcept
{ positionOverride = amount; positionOverridden = true; }
void EnsembleEngine::setStrikeAzimuthOverride (float amount) noexcept
{ azimuthOverride = amount; azimuthOverridden = true; }
void EnsembleEngine::clearStrikeOverrides() noexcept
{ positionOverridden = azimuthOverridden = false; }
} // namespace taikor
