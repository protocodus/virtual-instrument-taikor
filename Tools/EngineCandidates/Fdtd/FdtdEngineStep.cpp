// The time step of the finite-difference drum (see FdtdEngine.h for the
// equations and the phase structure).

#include "FdtdEngine.h"

#include <chrono>
#include <cmath>

namespace taikor::candidates::fdtd
{
namespace
{
inline int mod3 (long long n) noexcept { return static_cast<int> (((n % 3) + 3) % 3); }
} // namespace

template <typename Real>
void FdtdEngine<Real>::phase1Item (int thread, int item, long long n)
{
    const auto& it = items[static_cast<std::size_t> (item)];
    switch (it[0])
    {
        case 0: // nonlinear potential of head it[1] at w^n (member it[2] of its team), and the SAV vector
        {
            auto& head = heads[it[1]];
            const double* w = head.w[mod3 (n)].data();
            const double* wPrev = head.w[mod3 (n - 1)].data();
            if (head.teamSize > 1)
            {
                // Every member returns once the whole gradient is written
                // (VonKarmanPotential::evaluateTeam ends at the team's
                // barrier); member 0 goes on to the SAV vector, the others
                // back to the queue.
                const double V = head.team->evaluateMember (w, head.gradient.data(), it[2]);
                if (it[2] != 0)
                    break;
                head.V = V;
            }
            else
                head.V = head.potential->evaluate (w, head.gradient.data());
            const double scale = 1.0 / std::sqrt (2.0 * head.V + numerics.savEpsilon);
            double gg = 0.0, gc = 0.0, gw = 0.0;
            const auto& grid = head.grid;
            for (int j = 0; j < grid.rows(); ++j)
            {
                const int base = j * grid.grid().nx;
                for (int i = grid.rowBegin (j); i < grid.rowEnd (j); ++i)
                {
                    const auto l = static_cast<std::size_t> (base + i);
                    const double value = head.gradient[l] * scale;
                    head.g[l] = value;
                    gg += value * value;
                    gc += value * head.footprint[l];
                    gw += value * wPrev[l];
                }
            }
            head.gg = gg;
            head.gc = gc;
            head.gwPrev = gw;
            // The tension the SAV step applies this step (StretchingProbe.h).
            if (head.probe)
            {
                head.stretchTension = head.probe->peakStretchingTension();
                head.savTension = head.stretchTension * head.psi[(n + 1) & 1] * scale;
            }
            // J g (bg) is spread over the threads in P2a: row by row the same
            // sums, off this item's path (an in-plane solve makes it the
            // longest item of P1).
            break;
        }
        case 1: // the air planes on both sides of head it[1]
        {
            auto& head = heads[it[1]];
            const int planes = head.coupling.planeSize();
            for (int c = 0; c < planes; ++c)
            {
                if (head.hasAbove)
                    domain.stepMasked (n, head.cellAbove[static_cast<std::size_t> (c)]);
                if (head.hasBelow)
                    domain.stepMasked (n, head.cellBelow[static_cast<std::size_t> (c)]);
            }
            // The wall's ring of plane cells: its viscous corrections and the
            // thermal partial update, before the plane RHS reads them.
            if (! head.wallAbove.empty())
                for (int c = 0; c < planes; ++c)
                {
                    const auto u = static_cast<std::size_t> (c);
                    if (head.wallAbove[u] >= 0)
                        wall.slicePartial (n, head.wallAbove[u]);
                    if (head.wallBelow[u] >= 0)
                        wall.slicePartial (n, head.wallBelow[u]);
                }
            head.coupling.spread.multiply (head.w[mod3 (n - 1)].data(), head.jwPrev.data(), 0, planes);
            const Real* next = domain.level (n + 1);
            const Real* prev = domain.level (n - 1);
            for (int c = 0; c < planes; ++c)
            {
                const int above = head.cellAbove[static_cast<std::size_t> (c)], below = head.cellBelow[static_cast<std::size_t> (c)];
                const double qBelow = head.hasBelow ? static_cast<double> (next[below]) - prev[below] : 0.0;
                const double qAbove = head.hasAbove ? static_cast<double> (next[above]) - prev[above] : 0.0;
                if (head.sliceScale.empty())
                    head.y[static_cast<std::size_t> (c)] = head.loadScale * (qBelow - qAbove) + head.mu * head.jwPrev[static_cast<std::size_t> (c)];
                else // mu I S J w^{n-1} (HeadState::sliceScale)
                    head.y[static_cast<std::size_t> (c)] = head.loadScale * (qBelow - qAbove)
                                                           + head.mu * head.sliceScale[static_cast<std::size_t> (c)] * head.jwPrev[static_cast<std::size_t> (c)];
            }
            head.coupling.gram.multiply (head.y.data(), head.gy.data(), 0, planes);
            break;
        }
        case 2: // explicit part of head it[1], row chunk it[2]
        {
            auto& head = heads[it[1]];
            int j0, j1;
            splitRange (head.grid.rows(), head.chunks, it[2], j0, j1);
            auto& lap = head.scratch[static_cast<std::size_t> (thread)];
            const double* w = head.w[mod3 (n)].data();
            head.grid.laplacian (w, lap.data(), j0 - 1, j1 + 1);
            head.grid.explicitPart (w, head.w[mod3 (n - 1)].data(), lap.data(),
                                    head.lap[static_cast<std::size_t> ((n + 1) & 1)].data(), head.rExpl.data(), j0, j1);
            auto& keep = head.lap[static_cast<std::size_t> (n & 1)];
            const int nx = head.grid.grid().nx;
            for (int j = j0; j < j1; ++j)
                for (int i = head.grid.rowBegin (j); i < head.grid.rowEnd (j); ++i)
                    keep[static_cast<std::size_t> (j * nx + i)] = lap[static_cast<std::size_t> (j * nx + i)];
            break;
        }
        case 3: domain.stepRegular (n, it[1], it[2]); break;
        case 4: domain.stepIrregular (n, it[1], it[2]); break;
        case 5: domain.stepPml (n, it[1], it[2]); break;
        default: break;
    }
}

template <typename Real>
void FdtdEngine<Real>::phase2a (int thread, long long n)
{
    const int threads = pool->size();
    auto& slot = partials[static_cast<std::size_t> (thread)];
    slot.fill (0.0);
    for (int index = 0; index < 2; ++index)
    {
        auto& head = heads[index];
        int c0, c1;
        splitRange (head.coupling.planeSize(), threads, thread, c0, c1);
        head.coupling.spread.multiply (head.rExpl.data(), head.br.data(), c0, c1);
        for (int c = c0; c < c1; ++c)
            head.br[static_cast<std::size_t> (c)] += head.gy[static_cast<std::size_t> (c)];
        if (head.nonlinear)
        {
            head.coupling.spread.multiply (head.g.data(), head.bg.data(), c0, c1);
            int j0, j1;
            splitRange (head.grid.rows(), threads, thread, j0, j1);
            const int nx = head.grid.grid().nx;
            double sum = 0.0;
            for (int j = j0; j < j1; ++j)
                for (int i = head.grid.rowBegin (j); i < head.grid.rowEnd (j); ++i)
                    sum += head.g[static_cast<std::size_t> (j * nx + i)] * head.rExpl[static_cast<std::size_t> (j * nx + i)];
            slot[static_cast<std::size_t> (index * PartialCount + PartialGR)] = sum;
        }
    }
    if (thread == 0 && struckYet)
    {
        const auto& head = heads[0];
        const int nx = head.grid.grid().nx;
        double sum = 0.0;
        for (int j = head.footprintBegin; j < head.footprintEnd; ++j)
            for (int i = head.grid.rowBegin (j); i < head.grid.rowEnd (j); ++i)
                sum += head.footprint[static_cast<std::size_t> (j * nx + i)] * head.rExpl[static_cast<std::size_t> (j * nx + i)];
        slot[PartialCR] = sum;
    }
    // The barrel wall's boundary layer on every cell but the head planes'
    // (level n + 1 of the air is complete here, level n final).
    if (wallActive)
    {
        int e0, e1, c0, c1;
        splitRange (wall.frictionCount(), threads, thread, e0, e1);
        splitRange (wall.ordinaryCellCount(), threads, thread, c0, c1);
        wallLost[static_cast<std::size_t> (thread)] += wall.stepFriction (n, e0, e1) + wall.finishCells (n, c0, c1);
    }
}

template <typename Real>
void FdtdEngine<Real>::phase2b (int thread, long long)
{
    const int threads = pool->size();
    auto& slot = partials[static_cast<std::size_t> (thread)];
    for (int index = 0; index < 2; ++index)
    {
        auto& head = heads[index];
        int c0, c1;
        splitRange (head.coupling.planeSize(), threads, thread, c0, c1);
        head.coupling.inverse.multiply (head.br.data(), head.zr.data(), c0, c1);
        if (head.nonlinear)
            head.coupling.inverse.multiply (head.bg.data(), head.zg.data(), c0, c1);
        double s[PartialCount] = {};
        if (head.sliceScale.empty())
            for (int c = c0; c < c1; ++c)
            {
                const auto u = static_cast<std::size_t> (c);
                s[PartialGJY] += head.bg[u] * head.y[u];
                s[PartialCJY] += head.jc[u] * head.y[u];
                s[PartialGZR] += head.bg[u] * head.zr[u];
                s[PartialGZG] += head.bg[u] * head.zg[u];
                s[PartialCZR] += head.jc[u] * head.zr[u];
                s[PartialCZG] += head.jc[u] * head.zg[u];
                s[PartialGZC] += head.bg[u] * head.zc[u];
            }
        else // x . A^-1 v = (x . v - nu x . I S z) / alpha (Coupling.h setSliceScale)
            for (int c = c0; c < c1; ++c)
            {
                const auto u = static_cast<std::size_t> (c);
                const double sc = head.sliceScale[u];
                s[PartialGJY] += head.bg[u] * head.y[u];
                s[PartialCJY] += head.jc[u] * head.y[u];
                s[PartialGZR] += head.bg[u] * sc * head.zr[u];
                s[PartialGZG] += head.bg[u] * sc * head.zg[u];
                s[PartialCZR] += head.jc[u] * sc * head.zr[u];
                s[PartialCZG] += head.jc[u] * sc * head.zg[u];
                s[PartialGZC] += head.bg[u] * sc * head.zc[u];
            }
        for (int p = PartialGJY; p < PartialCount; ++p)
            slot[static_cast<std::size_t> (index * PartialCount + p)] = s[p];
    }
}

template <typename Real>
void FdtdEngine<Real>::phase3 (int thread, long long n)
{
    const int threads = pool->size();
    double sums[2 * PartialCount] = {};
    for (int t = 0; t < threads; ++t)
        for (int p = 0; p < 2 * PartialCount; ++p)
            sums[p] += partials[static_cast<std::size_t> (t)][static_cast<std::size_t> (p)];
    double aG[2] = {}, aC[2] = {}, gw1[2] = {};
    double f = 0.0;
    for (int index = 0; index < 2; ++index)
    {
        const auto& head = heads[index];
        const double* S = sums + index * PartialCount;
        const double q = head.coupling.transposeScale, nu = head.nu, alpha = head.alpha, beta = head.beta, gamma = head.gamma;
        const double gr0 = S[PartialGR] + q * S[PartialGJY];
        const double cr0 = S[PartialCR] + q * S[PartialCJY];
        const double gXr = (gr0 - nu * q * S[PartialGZR]) / alpha;
        const double gXg = (head.gg - nu * q * S[PartialGZG]) / alpha;
        const double gXc = (head.gc - nu * q * S[PartialGZC]) / alpha;
        const double cXr = (cr0 - nu * q * S[PartialCZR]) / alpha;
        const double cXg = (head.gc - nu * q * S[PartialCZG]) / alpha;
        const double psiPrev = head.psi[(n + 1) & 1];
        const double zeta = beta * head.gwPrev - gamma * psiPrev;
        const double denominator = 1.0 + beta * gXg;
        const double eR = zeta - beta * (gXr + zeta * gXg) / denominator;
        const double eC = -beta * gXc / denominator;
        double contactForce = 0.0;
        if (index == 0 && struckYet)
        {
            const double cx0 = cXr + eR * cXg, cxc = head.cXc + eC * cXg;
            const double yNow = stickY[mod3 (n)], yPrev = stickY[mod3 (n - 1)];
            const double etaFree = cx0 - 2.0 * yNow + yPrev;
            const double m = gamma * cxc + k * k / stickMass;
            const auto& wPrev = head.w[mod3 (n - 1)];
            const int nx = head.grid.grid().nx;
            double cw = 0.0;
            for (int j = head.footprintBegin; j < head.footprintEnd; ++j)
                for (int i = head.grid.rowBegin (j); i < head.grid.rowEnd (j); ++i)
                    cw += head.footprint[static_cast<std::size_t> (j * nx + i)] * wPrev[static_cast<std::size_t> (j * nx + i)];
            hertz.solve (cw - yPrev, etaFree, m, contactForce);
            f = contactForce;
        }
        aC[index] = -gamma * contactForce;
        aG[index] = eR - gamma * contactForce * eC;
        gw1[index] = gXr + aG[index] * gXg + aC[index] * gXc;
    }
    // w^{n+1} by rows and the air planes' coupling term by cells.
    for (int index = 0; index < 2; ++index)
    {
        auto& head = heads[index];
        double* next = head.w[mod3 (n + 1)].data();
        const auto& I = head.coupling.interpolate;
        const double invAlpha = 1.0 / head.alpha, nu = head.nu, ag = aG[index], ac = aC[index];
        int j0, j1;
        splitRange (head.grid.rows(), threads, thread, j0, j1);
        const int nx = head.grid.grid().nx;
        for (int j = j0; j < j1; ++j)
            for (int i = head.grid.rowBegin (j); i < head.grid.rowEnd (j); ++i)
            {
                const int l = j * nx + i;
                double air = 0.0;
                if (head.sliceScale.empty())
                    for (int e = I.start[static_cast<std::size_t> (l)]; e < I.start[static_cast<std::size_t> (l) + 1]; ++e)
                    {
                        const auto c = static_cast<std::size_t> (I.column[static_cast<std::size_t> (e)]);
                        air += I.value[static_cast<std::size_t> (e)] * (head.y[c] - nu * (head.zr[c] + ag * head.zg[c]));
                    }
                else // (v - nu I S z) / alpha (Coupling.h setSliceScale)
                    for (int e = I.start[static_cast<std::size_t> (l)]; e < I.start[static_cast<std::size_t> (l) + 1]; ++e)
                    {
                        const auto c = static_cast<std::size_t> (I.column[static_cast<std::size_t> (e)]);
                        air += I.value[static_cast<std::size_t> (e)]
                               * (head.y[c] - nu * (head.sliceScale[c] * (head.zr[c] + ag * head.zg[c])));
                    }
                next[l] = (head.rExpl[static_cast<std::size_t> (l)] + ag * head.g[static_cast<std::size_t> (l)] + air) * invAlpha
                          + ac * head.xc[static_cast<std::size_t> (l)];
            }
        if (options.withAir)
        {
            int c0, c1;
            splitRange (head.coupling.planeSize(), threads, thread, c0, c1);
            Real* level = domain.level (n + 1);
            const double flux = drumSpec.soundSpeed * drumSpec.soundSpeed * k / (2.0 * domain.delta);
            if (head.wallAbove.empty())
                for (int c = c0; c < c1; ++c)
                {
                    const auto u = static_cast<std::size_t> (c);
                    const double jw1 = (head.zr[u] + ag * head.zg[u] + ac * head.zc[u]) * invAlpha;
                    const double s = flux * (jw1 - head.jwPrev[u]);
                    if (head.hasAbove)
                        level[head.cellAbove[u]] = static_cast<Real> (level[head.cellAbove[u]] + s);
                    if (head.hasBelow)
                        level[head.cellBelow[u]] = static_cast<Real> (level[head.cellBelow[u]] - s);
                }
            else
            {
                // The wall's ring: the head's flux scaled by 1 / (1 + beta), then
                // the thermal branches from the final level n + 1.
                double lost = 0.0;
                for (int c = c0; c < c1; ++c)
                {
                    const auto u = static_cast<std::size_t> (c);
                    const double jw1 = (head.zr[u] + ag * head.zg[u] + ac * head.zc[u]) * invAlpha;
                    const double s = flux * (jw1 - head.jwPrev[u]);
                    if (head.hasAbove)
                        level[head.cellAbove[u]] = static_cast<Real> (level[head.cellAbove[u]] + s * head.factorAbove[u]);
                    if (head.hasBelow)
                        level[head.cellBelow[u]] = static_cast<Real> (level[head.cellBelow[u]] - s * head.factorBelow[u]);
                    if (head.wallAbove[u] >= 0)
                        lost += wall.sliceFinish (n, head.wallAbove[u]);
                    if (head.wallBelow[u] >= 0)
                        lost += wall.sliceFinish (n, head.wallBelow[u]);
                }
                wallLost[static_cast<std::size_t> (thread)] += lost;
            }
        }
    }
    if (thread != 0)
        return;
    queue.store (0, std::memory_order_relaxed);
    // Level n is final here (step n+1 has not started); the shift is applied
    // at the start of step n+1 to its two input levels n+1 and n.
    if (options.withAir && ((n + 1) & 63) == 0)
        interiorShift = -domain.interiorSample (n);
    force = f;
    for (int index = 0; index < 2; ++index)
    {
        auto& head = heads[index];
        // psi is advanced from the gradient, never recomputed from V, so it
        // keeps whatever offset from sqrt(2 V) the attack leaves; energy is
        // conserved with psi, not with V (Bilbao, Ducceschi & Zama 2023). For
        // von Karman the offset is measured below and discussed in
        // VonKarmanPotential.h ("SAV consistency"); the Berger letter is not
        // monitored here (tests/fdtd savTracking checks its tracking).
        const double psiBefore = head.psi[(n + 1) & 1];
        head.psi[n & 1] = psiBefore + 0.5 * (gw1[index] - head.gwPrev);
        if (head.monitorSav)
        {
            // psi^{n-1/2} against sqrt(V^n + V^{n-1} + eps) (HeadState::SavConsistency).
            auto& m = head.consistency;
            const double eps = numerics.savEpsilon;
            m.ratio = psiBefore / std::sqrt (head.V + head.vPrevious + eps);
            m.spurious = head.stretchTension * (m.ratio - 1.0);
            if (struckYet)
            {
                const double t = static_cast<double> (n) * k;
                double& worst = t < sizingWindow ? m.worstAttack : m.worstTail;
                double& when = t < sizingWindow ? m.worstAttackTime : m.worstTailTime;
                if (std::abs (m.spurious) > worst)
                {
                    worst = std::abs (m.spurious);
                    when = t;
                }
                m.worstOffset = std::max (m.worstOffset, std::abs (0.5 * (psiBefore * psiBefore - eps)
                                                                   - 0.5 * (head.V + head.vPrevious)));
            }
            head.vPrevious = head.V;
        }
        if (struckYet && head.probe)
        {
            if (head.savTension > head.peakSavTension)
            {
                head.peakSavTension = head.savTension;
                head.peakSavTensionTime = static_cast<double> (n) * k;
            }
            if (static_cast<double> (n) * k < sizingWindow)
                head.windowSavTension = std::max (head.windowSavTension, head.savTension);
        }
    }
    if (struckYet)
    {
        stickY[mod3 (n + 1)] = 2.0 * stickY[mod3 (n)] - stickY[mod3 (n - 1)] + k * k * f / stickMass;
        const double t = static_cast<double> (n + 1) * k;
        if (static_cast<double> (stats.force.size()) * k < sizingWindow)
            stats.force.push_back (f);
        if (f > 0.0)
        {
            if (! stats.inContact)
            {
                ++stats.episodes;
                if (stats.firstContact < 0.0)
                    stats.firstContact = t;
            }
            stats.inContact = true;
            stats.lastContact = t;
            stats.peakForce = std::max (stats.peakForce, f);
            stats.impulse += f * k;
        }
        else if (stats.inContact)
        {
            stats.inContact = false;
            stats.exitSpeed = (stickY[mod3 (n + 1)] - stickY[mod3 (n)]) / k;
        }
    }
    if (options.withAir && options.exteriorAir)
        for (int capsule = 0; capsule < 2; ++capsule)
            lastMic[capsule] = domain.pressure (n, micSpec.capsule[static_cast<std::size_t> (capsule)]);
}

// Float precision of the sealed interior. The stroke compresses the barrel's
// air, so its potential acquires an offset (1/V d/dt int Psi dV = -c^2 dV/V,
// about 1-2 m^2/s after a Don) that is 10^3-10^4 times the potential's
// change per step; in float32 that would leave only three or four digits of
// the interior pressure. Every 64 steps the interior is shifted back by its
// mean (AirDomain::shiftInterior: an exact symmetry of the scheme), so its
// values stay of the order of the fluctuation itself.
template <typename Real>
void FdtdEngine<Real>::recentre (int thread, long long n)
{
    if (! options.withAir || (n & 63) != 0 || interiorShift == 0.0)
        return;
    int c0, c1;
    splitRange (domain.interiorCount(), pool->size(), thread, c0, c1);
    domain.shiftInterior (n, interiorShift, c0, c1);
    barrier->wait();
}

template <typename Real>
void FdtdEngine<Real>::stepPhases (int thread, long long n)
{
    recentre (thread, n);
    if (profiling)
    {
        using clock = std::chrono::steady_clock;
        auto& p = profile[static_cast<std::size_t> (thread)];
        auto t = clock::now();
        const auto lap = [&t] { const auto now = clock::now(); const double s = std::chrono::duration<double> (now - t).count(); t = now; return s; };
        for (;;)
        {
            const int item = queue.fetch_add (1, std::memory_order_relaxed);
            if (item >= itemCount)
                break;
            phase1Item (thread, item, n);
            p[static_cast<std::size_t> (items[static_cast<std::size_t> (item)][0])] += lap();
        }
        barrier->wait();
        p[9] += lap();
        phase2a (thread, n);
        p[6] += lap();
        barrier->wait();
        p[9] += lap();
        phase2b (thread, n);
        p[7] += lap();
        barrier->wait();
        p[9] += lap();
        phase3 (thread, n);
        p[8] += lap();
        barrier->wait();
        p[9] += lap();
        return;
    }
    for (;;)
    {
        const int item = queue.fetch_add (1, std::memory_order_relaxed);
        if (item >= itemCount)
            break;
        phase1Item (thread, item, n);
    }
    barrier->wait();
    phase2a (thread, n);
    barrier->wait();
    phase2b (thread, n);
    barrier->wait();
    phase3 (thread, n);
    barrier->wait();
}

template <typename Real>
void FdtdEngine<Real>::runSteps (int steps, float* left, float* right)
{
    int written = 0;
    const long long first = step;
    pool->run ([&] (int thread)
    {
        for (int s = 0; s < steps; ++s)
        {
            stepPhases (thread, first + s);
            if (thread == 0 && left != nullptr)
            {
                double outL, outR;
                const bool dueL = decimator[0].push (lastMic[0], outL);
                decimator[1].push (lastMic[1], outR);
                if (dueL)
                {
                    left[written] = static_cast<float> (outL);
                    right[written] = static_cast<float> (outR);
                    ++written;
                }
            }
        }
    });
    step += steps;
}

template <typename Real>
void FdtdEngine<Real>::advance (int steps)
{
    runSteps (steps, nullptr, nullptr);
}

template <typename Real>
void FdtdEngine<Real>::process (float* left, float* right, int samples)
{
    if (! struckYet || ! failed.empty())
    {
        // Silence before the stroke, and after a failed rebuild (failure()),
        // which the render tool then reports as a silent take.
        for (int i = 0; i < samples; ++i)
            left[i] = right[i] = 0.0f;
        return;
    }
    const auto t0 = std::chrono::steady_clock::now();
    if (outputPreroll)
    {
        // The output FIR is linear-phase with a delay of D internal samples
        // (Decimator.h, D = 48, a multiple of the decimation factor). Running
        // the drum D steps ahead of its output and dropping the D / factor
        // outputs cancels that delay, so output sample j after the stroke is
        // the band-limited pressure at the stroke's time plus j output samples.
        outputPreroll = false;
        const int delay = decimator[0].delayInternalSamples();
        std::vector<float> droppedLeft (static_cast<std::size_t> (delay / numerics.oversample + 1), 0.0f);
        std::vector<float> droppedRight (droppedLeft);
        runSteps (delay, droppedLeft.data(), droppedRight.data());
        processedSteps += delay;
    }
    runSteps (samples * numerics.oversample, left, right);
    processSeconds += std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
    processedSteps += static_cast<long long> (samples) * numerics.oversample;
}

#define TAIKOR_ENGINE_INSTANTIATE(Real)                                               \
    template void FdtdEngine<Real>::phase1Item (int, int, long long);                 \
    template void FdtdEngine<Real>::phase2a (int, long long);                         \
    template void FdtdEngine<Real>::phase2b (int, long long);                         \
    template void FdtdEngine<Real>::phase3 (int, long long);                          \
    template void FdtdEngine<Real>::stepPhases (int, long long);                      \
    template void FdtdEngine<Real>::recentre (int, long long);                        \
    template void FdtdEngine<Real>::runSteps (int, float*, float*);                   \
    template void FdtdEngine<Real>::advance (int);                                    \
    template void FdtdEngine<Real>::process (float*, float*, int);
TAIKOR_ENGINE_INSTANTIATE (float)
TAIKOR_ENGINE_INSTANTIATE (double)
#undef TAIKOR_ENGINE_INSTANTIATE

} // namespace taikor::candidates::fdtd
