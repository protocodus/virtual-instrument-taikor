// Cost breakdown of the in-plane solver kernels on the 48 kHz grid (report
// only; no checks).

#include "InVacuoHead.h"

namespace vktest
{

void testProfile()
{
    const double fs = 48000.0;
    const double h = 1.001 * InVacuoHead::minimumSpacing (batter(), fs);
    Fixture f = makeFixture (headRadius(), h, 1.0e-8);
    if (! f.vk)
        return check (false, "prepare: " + f.error);
    auto& solver = const_cast<vk::InPlaneSolver&> (f.vk->solver());
    const int levels = solver.levelCount();
    std::string sizes;
    for (int l = 0; l < levels; ++l)
    {
        const auto& L = solver.level (l);
        long free = 0;
        for (int k = 0; k < L.N; ++k)
            free += L.freeU[static_cast<std::size_t> (k)] + L.freeV[static_cast<std::size_t> (k)];
        sizes += format ("L%d n=%d free=%ld taps=%zu rim rows=%zu lmax=%.3f; ", l, L.n, free, L.stencil[0].size(), L.irregularRow.size(), L.lambdaMax);
    }
    note (sizes);
    auto timeIt = [] (auto&& fn, int reps)
    {
        const double t0 = cpuSeconds();
        for (int r = 0; r < reps; ++r)
            fn();
        return (cpuSeconds() - t0) / reps * 1e6;
    };
    for (int l = 0; l < levels; ++l)
    {
        const auto& L = solver.level (l);
        std::vector<double> x (static_cast<std::size_t> (L.size()), 0.0), y (x.size(), 0.0);
        for (int k = 0; k < L.N; ++k)
        {
            if (L.freeU[static_cast<std::size_t> (k)]) x[static_cast<std::size_t> (L.base (0) + k)] = 1.0;
            if (L.freeV[static_cast<std::size_t> (k)]) x[static_cast<std::size_t> (L.base (1) + k)] = 1.0;
        }
        const double tApply = timeIt ([&] { solver.apply (l, x.data(), y.data()); }, 200);
        std::vector<float> xf (x.begin(), x.end()), yf (y.size());
        const double tApplyF = timeIt ([&] { solver.apply (l, xf.data(), yf.data()); }, 200);
        std::string line = format ("level %d: apply %.1f us (float %.1f us)", l, tApply, tApplyF);
        if (l > 0)
        {
            const auto& F = solver.level (l - 1);
            std::vector<float> ff (static_cast<std::size_t> (F.size()), 0.0f);
            line += format (", prolong %.1f us", timeIt ([&] { solver.prolong (l, xf.data(), ff.data()); }, 200));
            line += format (", restrict %.1f us", timeIt ([&] { solver.restrict_ (l, ff.data(), xf.data()); }, 200));
        }
        note (line);
    }
    const auto& L0 = solver.level (0);
    std::vector<double> r (static_cast<std::size_t> (L0.size()), 0.0), z (r.size(), 0.0);
    for (int k = 0; k < L0.N; ++k)
        if (L0.freeU[static_cast<std::size_t> (k)]) r[static_cast<std::size_t> (L0.base (0) + k)] = 1.0;
    note (format ("V-cycle (with double<->float conversion): %.1f us", timeIt ([&] { solver.precondition (r.data(), z.data()); }, 50)));

    // Warm-start quality over the first 20 ms of the test stroke.
    InVacuoHead head;
    head.init (f.grid, batter(), fs, f.vk.get());
    head.setDisplacement (std::vector<double> (static_cast<std::size_t> (f.grid.size()), 0.0));
    head.armStick (shipping().stroke);
    for (int order : { 0, 1, 2 })
    {
        f.vk->setWarmStartOrder (order);
        f.vk->resetStats();
        head.setDisplacement (std::vector<double> (static_cast<std::size_t> (f.grid.size()), 0.0));
        head.armStick (shipping().stroke);
        f.vk->resetWarmStart();
        const double t0 = cpuSeconds();
        for (int n = 0; n < 960; ++n)
            head.step();
        const auto& st = f.vk->stats();
        note (format ("warm order %d: initial residual mean %.2e max %.2e; PCG its mean %.2f max %d; %.0f us per step",
                      order, st.initialResidualSum / static_cast<double> (st.calls), st.initialResidualMax,
                      static_cast<double> (st.iterations) / static_cast<double> (st.calls), st.maxIterations,
                      (cpuSeconds() - t0) / 960 * 1e6));
    }
}

} // namespace vktest
