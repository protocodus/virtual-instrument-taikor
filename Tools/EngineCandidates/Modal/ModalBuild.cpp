#include "ModalInternal.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <sstream>
#include <thread>

namespace taikor::candidates::modal
{

void parallelFor (int count, int threads, const std::function<void (int)>& body)
{
    std::atomic<int> next { 0 };
    const auto worker = [&]
    {
        for (int i = next++; i < count; i = next++)
            body (i);
    };
    const int n = std::max (1, std::min (threads, count));
    std::vector<std::thread> pool;
    for (int t = 1; t < n; ++t)
        pool.emplace_back (worker);
    worker();
    for (auto& t : pool)
        t.join();
}

bool prepareContext (const DrumSpec& drum, const BuildOptions& o, DrumContext& ctx, std::string& error, const MicSpec* mics)
{
    const auto material = [&] (const HeadSpec& h)
    {
        HeadMaterial m;
        m.tension = h.tension;
        m.arealDensity = h.arealDensity;
        m.bending = h.bendingStiffness;
        m.radius = drum.radius;
        return m;
    };
    if (drum.radius <= 0.0 || drum.depth <= 0.0 || drum.airDensity <= 0.0 || drum.soundSpeed <= 0.0)
    {
        error = "drum geometry and air must be positive";
        return false;
    }
    for (const HeadSpec* h : { &drum.batter, &drum.rear })
        if (h->tension <= 0.0 || h->arealDensity <= 0.0 || h->bendingStiffness <= 0.0)
        {
            error = "head tension, areal density and bending stiffness must be positive";
            return false;
        }
    ctx.options = o;
    ctx.air = { drum.radius, drum.depth, drum.airDensity, drum.soundSpeed };
    ctx.batterMaterial = material (drum.batter);
    ctx.rearMaterial = material (drum.rear);
    ctx.batterLoss = { drum.batter.sigma0, drum.batter.sigma1, drum.batter.arealDensity };
    ctx.rearLoss = { drum.rear.sigma0, drum.rear.sigma1, drum.rear.arealDensity };
    ctx.omegaCoupled = 2.0 * pi * o.coupledMaxHz;
    ctx.omegaAir = 2.0 * pi * o.airMaxHz;
    ctx.batter = solveHeadModes (ctx.batterMaterial, 2.0 * pi * o.headMaxHz);
    // Rear modes above the coupled band would be lone, never excited modes.
    ctx.rear = solveHeadModes (ctx.rearMaterial, ctx.omegaCoupled);
    const int orders = static_cast<int> (std::max (ctx.batter.size(), ctx.rear.size()));
    ctx.batterGradient.resize (static_cast<std::size_t> (orders));
    ctx.rearGradient.resize (static_cast<std::size_t> (orders));
    ctx.cross.resize (static_cast<std::size_t> (orders));
    const double a = drum.radius;
    parallelFor (orders, o.threads, [&] (int m)
    {
        const auto um = static_cast<std::size_t> (m);
        double b1 = 0.0;
        if (um < ctx.batter.size())
        {
            ctx.batterGradient[um] = gradientForm (ctx.batter[um], ctx.batterMaterial, o.gradientTolerance);
            b1 = std::max (b1, ctx.batter[um].back().b1);
        }
        if (um < ctx.rear.size())
        {
            ctx.rearGradient[um] = gradientForm (ctx.rear[um], ctx.rearMaterial, o.gradientTolerance);
            b1 = std::max (b1, ctx.rear[um].back().b1);
        }
        // Cross-modes to 3 b1 + 40/a: the interior sums over q converge as
        // gamma^-4 beyond the head wavenumbers (tests/modal cavity_convergence).
        const double gammaMax = std::max (o.crossGammaFactor * b1 + 40.0 / a, ctx.omegaAir / drum.soundSpeed + 40.0 / a);
        if (o.interior)
            ctx.cross[um] = crossModes (m, a, gammaMax);
    });
    double b1Max = 0.0;
    for (const auto& set : ctx.batter)
        b1Max = std::max (b1Max, set.back().b1);
    for (const auto& set : ctx.rear)
        b1Max = std::max (b1Max, set.back().b1);
    ctx.exterior = std::make_unique<ExteriorTable> (a, 12.0 * b1Max + 60.0 / a, orders + 1);
    ctx.wall = std::max (drum.shellWallThickness, 0.0);
    if (o.exterior && o.bodyExterior)
    {
        // Columns: the coupled-band head modes below the table's top, in the
        // blocks' order (ascending frequency per order).
        const double top = 2.0 * pi * (o.bodyMaxHz > 0.0 ? o.bodyMaxHz : 1.3 * o.coupledMaxHz);
        std::vector<std::vector<const HeadMode*>> b (static_cast<std::size_t> (orders)), r (static_cast<std::size_t> (orders));
        for (int m = 0; m < orders; ++m)
        {
            const auto um = static_cast<std::size_t> (m);
            if (um < ctx.batter.size())
                for (const auto& mode : ctx.batter[um])
                    if (mode.omega < std::min (ctx.omegaCoupled, top))
                        b[um].push_back (&mode);
            if (um < ctx.rear.size())
                for (const auto& mode : ctx.rear[um])
                    if (mode.omega < std::min (ctx.omegaCoupled, top))
                        r[um].push_back (&mode);
        }
        while (! b.empty() && b.back().empty() && r.back().empty())
        {
            b.pop_back();
            r.pop_back();
        }
        BodyExteriorOptions bo;
        bo.maxHz = o.bodyMaxHz > 0.0 ? o.bodyMaxHz : 1.3 * o.coupledMaxHz;
        bo.gridPerOctave = o.gridPerOctave;
        bo.panelsPerWavelength = o.bodyPanelsPerWavelength;
        bo.chiefPoints = o.bodyChiefPoints;
        bo.threads = o.threads;
        if (mics != nullptr)
        {
            // Capsules as field points; the reference is the engine's own transfer of a
            // column (ModalMics.cpp): baffled Rayleigh for the batter, the compact
            // monopole for the rear (order 0; none for higher orders), per unit velocity.
            std::array<std::array<double, 2>, 2> points {};
            for (std::size_t j = 0; j < 2; ++j)
            {
                points[j] = { std::hypot (mics->capsule[j][0], mics->capsule[j][1]), mics->capsule[j][2] };
                bo.fieldPoints.push_back (points[j]);
            }
            const AirSpec air = ctx.air;
            const std::array<std::array<double, 3>, 2> capsules = mics->capsule;
            bo.reference = [b, r, air, points, capsules] (int m, int column, double omega, int point)
            {
                const auto um = static_cast<std::size_t> (m);
                const std::complex<double> jw (0.0, omega);
                const auto nbc = static_cast<int> (b[um].size());
                if (column < nbc)
                {
                    const HeadMode& h = *b[um][static_cast<std::size_t> (column)];
                    const auto shape = [&h, &air] (double rr) { return evaluateShape (h, air.radius, rr).value; };
                    return air.density * jw * rayleighIntegral (m, air.radius, shape, h.b1, jw, air.soundSpeed,
                                                                points[static_cast<std::size_t> (point)][0],
                                                                points[static_cast<std::size_t> (point)][1], 1);
                }
                if (m != 0)
                    return std::complex<double> (0.0);
                const HeadMode& h = *r[um][static_cast<std::size_t> (column - nbc)];
                const double path = rearDistance (air.length, capsules[static_cast<std::size_t> (point)]);
                return air.density * jw * h.volume * std::exp (-jw * (path / air.soundSpeed)) / (4.0 * pi * path);
            };
        }
        ctx.body = std::make_unique<BodyExterior> (ctx.air, ctx.wall, b, r, bo);
    }
    return true;
}

std::unique_ptr<ModalModel> buildModel (const DrumSpec& drum, const MicSpec& mics, const BuildOptions& o,
                                        std::string& error)
{
    const auto t0 = std::chrono::steady_clock::now();
    DrumContext ctx;
    if (! prepareContext (drum, o, ctx, error, &mics))
        return nullptr;
    const auto t1 = std::chrono::steady_clock::now();
    const int orders = static_cast<int> (ctx.batterGradient.size());
    std::vector<std::vector<CoupledMode>> coupled (static_cast<std::size_t> (orders)), high (static_cast<std::size_t> (orders));
    // Low orders carry the largest blocks: they are taken first.
    parallelFor (orders, o.threads, [&] (int m) { solveCoupledBlock (ctx, m, coupled[static_cast<std::size_t> (m)]); });
    const auto t2 = std::chrono::steady_clock::now();
    parallelFor (orders, o.threads, [&] (int m) { solveHighModes (ctx, m, high[static_cast<std::size_t> (m)]); });
    const auto t3 = std::chrono::steady_clock::now();

    auto model = std::make_unique<ModalModel>();
    for (int m = 0; m < orders; ++m)
    {
        auto& block = coupled[static_cast<std::size_t> (m)];
        if (! block.empty())
            ++model->coupledBlocks;
        const int base = static_cast<int> (model->modes.size());
        for (auto& mode : block)
        {
            for (auto& link : mode.damping)
                link.first += base;
            model->dampingLinks += static_cast<int> (mode.damping.size());
            model->modes.push_back (std::move (mode));
        }
    }
    model->dampingLinks /= 2;
    for (auto& set : high)
        for (auto& mode : set)
        {
            model->modes.push_back (std::move (mode));
            ++model->highModes;
        }
    for (int m = 0; m < orders && o.interior; ++m)
        model->airModes += static_cast<int> (airModes (m, ctx.cross[static_cast<std::size_t> (m)], ctx.air, ctx.omegaAir).size());
    parallelFor (static_cast<int> (model->modes.size()), o.threads, [&] (int i)
    {
        computeMicTransfers (ctx, mics, model->modes[static_cast<std::size_t> (i)]);
        computeBodyTransfers (ctx, model->modes[static_cast<std::size_t> (i)]);
    });
    const auto t4 = std::chrono::steady_clock::now();

    for (const auto& set : ctx.batter)
        model->batterCount += static_cast<int> (set.size());
    for (const auto& set : ctx.rear)
        model->rearCount += static_cast<int> (set.size());
    if (ctx.body)
    {
        model->bodyPanels = ctx.body->panels();
        model->bodyGridPoints = static_cast<int> (ctx.body->grid().size());
        model->bodyMaxHz = ctx.body->maxOmega() / (2.0 * pi);
        model->bodySeconds = ctx.body->seconds();
        model->bodyAsymmetry = ctx.body->maxAsymmetry();
    }
    model->radius = drum.radius;
    model->batterMaterial = ctx.batterMaterial;
    model->rearMaterial = ctx.rearMaterial;
    model->batter = std::move (ctx.batter);
    model->rear = std::move (ctx.rear);
    model->batterGradient = std::move (ctx.batterGradient);
    model->rearGradient = std::move (ctx.rearGradient);
    const auto secs = [] (auto x, auto y) { return std::chrono::duration<double> (y - x).count(); };
    model->seconds = secs (t0, t4);
    std::ostringstream log;
    log.precision (3);
    log << std::fixed << "precompute " << secs (t0, t4) << " s (heads/catalogue/body exterior " << secs (t0, t1)
        << ", coupled blocks " << secs (t1, t2) << ", lone batter modes " << secs (t2, t3)
        << ", microphone transfers " << secs (t3, t4) << ") on " << o.threads << " threads";
    model->log = log.str();
    return model;
}

BlockSummary solveBlockForTest (const DrumSpec& drum, const BuildOptions& o, int m)
{
    BlockSummary summary;
    DrumContext ctx;
    std::string error;
    if (! prepareContext (drum, o, ctx, error))
        return summary;
    std::vector<CoupledMode> modes;
    solveCoupledBlock (ctx, m, modes);
    std::sort (modes.begin(), modes.end(), [] (const auto& x, const auto& y) { return x.omega < y.omega; });
    for (const auto& mode : modes)
    {
        summary.omega.push_back (mode.omega);
        summary.alpha.push_back (mode.alpha);
        summary.shareBatter.push_back (mode.shareBatter);
        summary.shareRear.push_back (mode.shareRear);
        summary.shareAir.push_back (mode.shareAir);
    }
    return summary;
}

} // namespace taikor::candidates::modal
