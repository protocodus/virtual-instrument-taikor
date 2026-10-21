// Building the finite-difference drum: grids, coupling, footprint, threads.

#include "FdtdEngine.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <thread>

namespace taikor::candidates::fdtd
{
namespace
{
inline std::size_t mod3Of (long long n) noexcept { return static_cast<std::size_t> (((n % 3) + 3) % 3); }
} // namespace

template <typename Real>
FdtdEngine<Real>::FdtdEngine (std::unique_ptr<NonlinearPotential> batterPotential, FdtdOptions opts)
    : options (opts)
{
    heads[0].potential = std::move (batterPotential);
    heads[0].struck = true;
}

template <typename Real>
FdtdEngine<Real>::~FdtdEngine() = default;

template <typename Real>
std::string FdtdEngine<Real>::description() const
{
    const std::string batter = heads[0].potential ? heads[0].potential->name() : std::string ("none");
    // A batter whose stretching the energy does not bound (von Karman) gets
    // its grid from sizing runs of the stroke (FdtdEngineSizing.cpp); say so.
    // description() is read before prepare(), so this asks the potential.
    const auto* probe = dynamic_cast<const StretchingProbe*> (heads[0].potential.get());
    const std::string sizing = probe != nullptr && ! probe->boundedByEnergy()
                                   ? "; the batter's grid is sized at the stroke for the peak stretching tension "
                                     "measured by runs of that stroke (not by the energy bound the Berger letter's "
                                     "grid uses)"
                                   : "";
    const std::string air = options.wallLoss && ! options.lossless
                                ? "a sealed rigid barrel whose interior side wall carries Kirchhoff's viscothermal "
                                  "boundary layer (G. Kirchhoff 1868; A. D. Pierce, Acoustics, sec. 10.5) as a passive "
                                  "sqrt(j w) boundary of diffusive first-order branches, energy-consistent with the "
                                  "leapfrog air (S. Bilbao, B. Hamilton, J. Botts & L. Savioja, IEEE/ACM TASLP 24(1) "
                                  "2016); the air is otherwise lossless (no bulk atmospheric absorption: it needs a "
                                  "humidity the drum spec does not carry)"
                                : "lossless air in a sealed rigid barrel";
    return "3-D finite-difference drum (S. Bilbao, JASA 131(1) 2012; S. Bilbao & C. Webb, JAES 61(10) 2013; "
           "A. Torin, PhD thesis, Edinburgh 2016): two clamped stiff lossy hide heads on staircase grids in "
           "leapfrog 3-D air (barrel interior and exterior, rigid staircase barrel, Grote-Sim PML), "
           "energy-consistent head-air coupling, batter nonlinearity " + batter
           + " and rear Berger via explicit SAV (Bilbao, Ducceschi & Zama, JCP 472 2023), energy-conserving "
             "Hertz bachi collision (Bilbao, Torin & Chatziioannou, Acta Acust. 101 2015), pressure capsules; "
             "output limited to the air grid's valid band (flat to ~4.9 kHz, nothing above its ~8.9 kHz on-axis "
             "cut-off); " + air + sizing + ".";
}

template <typename Real>
void FdtdEngine<Real>::fail (const std::string& what, const std::string& error)
{
    if (failed.empty())
        failed = what + ": " + (error.empty() ? std::string ("unknown error") : error);
}

// The largest tension the SAV step can apply to head `index` with headEnergy
// joules in the drum (FdtdDrumConfig.h bergerTensionBound). It sets both
// heads' spacing whether or not their nonlinearity is on, so a run with the
// nonlinearity off differs from the letter in that alone.
//
// A batter potential without such a bound (von Karman, StretchingProbe.h)
// has its grid built for the design tension its sizing runs measured
// (FdtdEngineSizing.cpp) instead.
template <typename Real>
double FdtdEngine<Real>::tensionBound (int index) const
{
    if (index == 0 && batterDesign > 0.0)
        return batterDesign;
    return bergerTensionBound (index == 0 ? drumSpec.batter : drumSpec.rear, drumSpec.radius, headEnergy);
}

// batterOnly: rebuild the batter alone (its sizing), keeping the rear head.
template <typename Real>
void FdtdEngine<Real>::buildHeads (const DrumSpec& drum, std::string& error, bool batterOnly)
{
    const HeadSpec* specs[2] = { &drum.batter, &drum.rear };
    const double known[2] = { batterKnownInclusion > 0.0 ? batterKnownInclusion : options.batterInclusion,
                              options.rearInclusion };
    const double forced[2] = { options.batterSpacing, options.rearSpacing };
    double spacing[2];
    for (int index = 0; index < 2; ++index)
        spacing[index] = forced[index] > 0.0
                             ? forced[index]
                             : std::max (headSpacing (*specs[index], k, tensionBound (index)), options.minimumHeadSpacing);
    // The staircase correction of each head (HeadCalibration.h), both heads at once.
    std::thread rear ([&]
    {
        if (batterOnly)
            return;
        auto& head = heads[1];
        if (known[1] > 0.0)
            head.calibration.inclusionRadius = known[1];
        else
            head.calibration = calibrateInclusionRadius (*specs[1], drum.radius, k, spacing[1]);
    });
    {
        auto& head = heads[0];
        if (known[0] > 0.0)
            head.calibration.inclusionRadius = known[0];
        else
            head.calibration = calibrateInclusionRadius (*specs[0], drum.radius, k, spacing[0]);
    }
    rear.join();
    for (int index = 0; index < (batterOnly ? 1 : 2); ++index)
    {
        auto& head = heads[index];
        const HeadSpec& spec = *specs[index];
        head.grid.build (spec, drum.radius, k, spacing[index], head.calibration.inclusionRadius, options.lossless);
        const auto n = static_cast<std::size_t> (head.grid.size());
        for (auto& level : head.w)
            level.assign (n, 0.0);
        for (auto& level : head.lap)
            level.assign (n, 0.0);
        head.rExpl.assign (n, 0.0);
        head.gradient.assign (n, 0.0);
        head.g.assign (n, 0.0);
        head.footprint.assign (n, 0.0);
        head.xc.assign (n, 0.0);
        head.nonlinear = index == 0 ? options.withBatterNonlinearity : options.bergerOnRear;
        if (index == 1)
            head.potential = makeBergerPotential();
        if (head.nonlinear && head.potential
            && ! head.potential->prepare (head.grid.grid(), spec, drum.radius, error))
            return;
        if (! head.potential)
            head.nonlinear = false;
        head.probe = head.nonlinear ? dynamic_cast<const StretchingProbe*> (head.potential.get()) : nullptr;
        head.monitorSav = head.probe != nullptr && ! head.probe->boundedByEnergy();
        head.consistency = {};
        head.stretchTension = head.vPrevious = 0.0;
        // A batter potential a team of threads can share (StretchingProbe.h):
        // set its team after every prepare(), so the sizing runs, the self-check
        // twin and the take all evaluate with the same team (same rounding).
        head.team = index == 0 && head.nonlinear ? dynamic_cast<TeamPotential*> (head.potential.get()) : nullptr;
        head.teamSize = head.team != nullptr ? head.team->setTeam (batterTeam) : 1;
        head.savTension = head.peakSavTension = head.peakSavTensionTime = head.windowSavTension = 0.0;
        head.alpha = head.grid.alpha();
        const double h = head.grid.spacing();
        head.beta = k * k / (4.0 * spec.arealDensity * h * h);
        head.gamma = k * k / (spec.arealDensity * h * h);
        head.psi[0] = head.psi[1] = std::sqrt (numerics.savEpsilon);
        // Row chunks of the explicit part: P1 work items, fixed at prepare().
        head.chunks = std::max (2, threadCount / 2);
    }
}

template <typename Real>
void FdtdEngine<Real>::buildAir (const DrumSpec& drum)
{
    AirPlaneLayout& layout = airLayout;
    std::vector<std::array<int, 2>>& columns = airColumns;
    layout = {};
    columns.clear();
    if (options.withAir)
    {
        AirGeometry geometry;
        geometry.radius = drum.radius;
        geometry.depth = drum.depth;
        geometry.wall = drum.shellWallThickness;
        geometry.delta = drum.depth / numerics.airCellsPerDepth;
        geometry.margin = options.exteriorAir ? numerics.airMargin : 0.0;
        geometry.marginAbove = options.exteriorAir ? std::max (0.0, numerics.airMarginAbove - numerics.airMargin) : 0.0;
        geometry.pmlCells = options.rigidBox || ! options.exteriorAir ? 0 : numerics.pmlCells;
        geometry.pmlReflection = numerics.pmlReflection;
        domain.build (geometry, drum.soundSpeed, drum.airDensity, k);
        columns = domain.headColumns (drum.radius);
        domain.closeHeadFaces (domain.planeAbove (0.0), columns, true);
        domain.closeHeadFaces (domain.planeAbove (-drum.depth), columns, true);
        domain.finalise();
        domain.listInterior (drum.radius, drum.depth);
        layout = { domain.x0, domain.y0, domain.delta, domain.nx, domain.ny };
    }
    // The barrel wall's boundary layer (orchestrator decision: every engine
    // applies Kirchhoff's side-wall loss; WallLoss.h), not in a lossless drum.
    wall = WallLoss<Real>();
    if (options.withAir && options.wallLoss && ! options.lossless)
    {
        wallFit = fitSquareRoot (k);
        wall.build (domain, drum.radius, drum.depth, wallFit);
        wall.setTracking (options.trackWallLoss);
    }
    wallActive = wall.active();
}

// Each head's air planes and its coupling operators (Coupling.h); rebuilt
// with the heads when a stroke needs wider head grids.
template <typename Real>
void FdtdEngine<Real>::buildHeadCoupling (const DrumSpec& drum)
{
    const AirPlaneLayout& layout = airLayout;
    const std::vector<std::array<int, 2>>& columns = airColumns;
    const double zFace[2] = { 0.0, -drum.depth };
    for (int index = 0; index < 2; ++index)
    {
        auto& head = heads[index];
        const HeadSpec& spec = index == 0 ? drum.batter : drum.rear;
        const int planes = static_cast<int> (columns.size());
        head.cellAbove.assign (static_cast<std::size_t> (planes), 0);
        head.cellBelow.assign (static_cast<std::size_t> (planes), 0);
        head.hasAbove = head.hasBelow = false;
        if (options.withAir)
        {
            const int kz = domain.planeAbove (zFace[index]);
            for (int c = 0; c < planes; ++c)
            {
                head.cellAbove[static_cast<std::size_t> (c)] = domain.index (columns[static_cast<std::size_t> (c)][0], columns[static_cast<std::size_t> (c)][1], kz);
                head.cellBelow[static_cast<std::size_t> (c)] = domain.index (columns[static_cast<std::size_t> (c)][0], columns[static_cast<std::size_t> (c)][1], kz - 1);
            }
            head.hasAbove = planes > 0 && domain.type (head.cellAbove[0]) == CellType::Slice;
            head.hasBelow = planes > 0 && domain.type (head.cellBelow[0]) == CellType::Slice;
            // Each air side contributes rho c^2 k^2 / (4 sigma Delta) (Coupling.h).
            const int sides = (head.hasAbove ? 1 : 0) + (head.hasBelow ? 1 : 0);
            head.mu = sides * drum.airDensity * drum.soundSpeed * drum.soundSpeed * k * k
                      / (4.0 * spec.arealDensity * domain.delta);
            head.loadScale = k * drum.airDensity / (2.0 * spec.arealDensity);
        }
        head.nu = head.mu / head.alpha;
        head.coupling.build (head.grid, layout, columns, head.nu, numerics.neumannOrder);
        for (auto* v : { &head.jwPrev, &head.y, &head.gy, &head.br, &head.bg, &head.zr, &head.zg, &head.zc, &head.jc })
            v->assign (static_cast<std::size_t> (planes), 0.0);
        // The wall's ring of plane cells: the implicit thermal term scales the
        // head's flux into each by 1 / (1 + beta), and the head solve sees the
        // mean over its air sides (FdtdEngine.h HeadState::sliceScale).
        head.wallAbove.clear();
        head.wallBelow.clear();
        head.factorAbove.clear();
        head.factorBelow.clear();
        head.sliceScale.clear();
        if (wallActive)
        {
            head.wallAbove.assign (static_cast<std::size_t> (planes), -1);
            head.wallBelow.assign (static_cast<std::size_t> (planes), -1);
            head.factorAbove.assign (static_cast<std::size_t> (planes), 1.0);
            head.factorBelow.assign (static_cast<std::size_t> (planes), 1.0);
            head.sliceScale.assign (static_cast<std::size_t> (planes), 1.0);
            const int sides = (head.hasAbove ? 1 : 0) + (head.hasBelow ? 1 : 0);
            for (int c = 0; c < planes; ++c)
            {
                const auto u = static_cast<std::size_t> (c);
                if (head.hasAbove && (head.wallAbove[u] = wall.sliceEntry (head.cellAbove[u])) >= 0)
                    head.factorAbove[u] = wall.sliceFactor (head.wallAbove[u]);
                if (head.hasBelow && (head.wallBelow[u] = wall.sliceEntry (head.cellBelow[u])) >= 0)
                    head.factorBelow[u] = wall.sliceFactor (head.wallBelow[u]);
                head.sliceScale[u] = ((head.hasAbove ? head.factorAbove[u] : 0.0) + (head.hasBelow ? head.factorBelow[u] : 0.0))
                                     / std::max (1, sides);
            }
            head.coupling.setSliceScale (head.sliceScale, numerics.neumannOrder);
        }
    }
}

// The footprint of the stick tip (a uniform disk, shared decision 3) as
// coverage weights of the nodes' cells, summing to one; then the fixed
// solves A^-1 c, J c and (I + nu G)^-1 J c the contact needs every step.
template <typename Real>
void FdtdEngine<Real>::precompute()
{
    auto& head = heads[0];
    const auto& g = head.grid.grid();
    const double h = head.grid.spacing(), b = strokeSpec.footprintRadius;
    // Numerical: 16 x 16 sub-samples per cell (coverage to ~1/256 of a cell).
    constexpr int sub = 16;
    double total = 0.0;
    int rowMin = g.ny, rowMax = -1;
    std::fill (head.footprint.begin(), head.footprint.end(), 0.0);
    for (int j = 0; j < g.ny; ++j)
        for (int i = head.grid.rowBegin (j); i < head.grid.rowEnd (j); ++i)
        {
            int inside = 0;
            for (int sy = 0; sy < sub; ++sy)
                for (int sx = 0; sx < sub; ++sx)
                {
                    const double x = g.x (i) + h * ((sx + 0.5) / sub - 0.5) - strokeSpec.strikeX;
                    const double y = g.y (j) + h * ((sy + 0.5) / sub - 0.5) - strokeSpec.strikeY;
                    inside += x * x + y * y < b * b ? 1 : 0;
                }
            if (inside > 0)
            {
                const double weight = static_cast<double> (inside) / (sub * sub);
                head.footprint[static_cast<std::size_t> (g.index (i, j))] = weight;
                total += weight;
                rowMin = std::min (rowMin, j);
                rowMax = std::max (rowMax, j);
            }
        }
    for (auto& v : head.footprint)
        v /= total;
    head.footprintCells = total;
    head.footprintBegin = rowMin;
    head.footprintEnd = rowMax + 1;
    const int planes = head.coupling.planeSize();
    head.coupling.spread.multiply (head.footprint.data(), head.jc.data(), 0, planes);
    head.coupling.inverse.multiply (head.jc.data(), head.zc.data(), 0, planes);
    // xc = (c - nu I S zc) / alpha (S = I without the wall loss).
    std::vector<double> izc (head.footprint.size(), 0.0);
    if (head.sliceScale.empty())
        head.coupling.interpolate.multiply (head.zc.data(), izc.data(), 0, g.size());
    else
    {
        std::vector<double> szc (head.zc.size());
        for (std::size_t c = 0; c < szc.size(); ++c)
            szc[c] = head.sliceScale[c] * head.zc[c];
        head.coupling.interpolate.multiply (szc.data(), izc.data(), 0, g.size());
    }
    for (std::size_t l = 0; l < izc.size(); ++l)
        head.xc[l] = head.grid.isInterior (static_cast<int> (l))
                         ? (head.footprint[l] - head.nu * izc[l]) / head.alpha : 0.0;
    head.cXc = head.grid.dot (head.footprint.data(), head.xc.data());
}

template <typename Real>
bool FdtdEngine<Real>::prepare (const DrumSpec& drum, const MicSpec& mics, double sampleRate, std::string& error)
{
    const auto t0 = std::chrono::steady_clock::now();
    numerics = FdtdNumerics::fromEnvironment();
    drumSpec = drum;
    micSpec = mics;
    outputRate = sampleRate;
    k = 1.0 / (sampleRate * numerics.oversample);
    headEnergy = options.strokeEnergy;
    batterDesign = options.batterDesignTension;
    batterKnownInclusion = 0.0;
    threadCount = options.threads > 0 ? options.threads : candidateThreadCount();
    // Every member of the batter's team is a P1 item held by its own thread
    // until the team finishes, so the team can never exceed the pool.
    batterTeam = std::clamp (options.potentialTeam > 0 ? options.potentialTeam : automaticTeamSize (threadCount), 1,
                             threadCount);
    buildHeads (drum, error);
    if (! error.empty())
        return false;
    buildAir (drum);
    buildHeadCoupling (drum);
    const int threads = threadCount;
    pool = std::make_unique<ThreadPool> (threads);
    barrier = std::make_unique<Barrier> (threads);
    partials.assign (static_cast<std::size_t> (threads), {});
    wallLost.assign (static_cast<std::size_t> (threads), 0.0);
    profiling = std::getenv ("TAIKOR_FDTD_PROFILE") != nullptr;
    profile.assign (static_cast<std::size_t> (threads), {});
    for (auto& head : heads)
        head.scratch.assign (static_cast<std::size_t> (threads), std::vector<double> (head.w[0].size(), 0.0));
    // P1 work items {kind, a, b}: potentials first (the longest; a team's
    // members {0, head, member} consecutively, so distinct threads take them
    // and start together), then the air planes, the heads' explicit chunks
    // and the air.
    items.clear();
    for (int index = 0; index < 2; ++index)
        if (heads[index].nonlinear)
            for (int member = 0; member < heads[index].teamSize; ++member)
                items.push_back ({ 0, index, member });
    if (options.withAir)
        for (int index = 0; index < 2; ++index)
            items.push_back ({ 1, index, 0 });
    for (int index = 0; index < 2; ++index)
        for (int chunk = 0; chunk < heads[index].chunks; ++chunk)
            items.push_back ({ 2, index, chunk });
    if (options.withAir)
    {
        const int planeChunks = 3 * threads, pmlChunks = 2 * threads;
        for (int c = 0; c < planeChunks; ++c)
        {
            int a, b;
            splitRange (domain.nz, planeChunks, c, a, b);
            items.push_back ({ 3, a, b });
        }
        for (int c = 0; c < threads; ++c)
        {
            int a, b;
            splitRange (domain.irregularCount(), threads, c, a, b);
            items.push_back ({ 4, a, b });
        }
        if (domain.hasPml())
        {
            const auto bounds = domain.pmlPlaneSplit (pmlChunks);
            for (int c = 0; c < pmlChunks; ++c)
                if (bounds[static_cast<std::size_t> (c + 1)] > bounds[static_cast<std::size_t> (c)])
                    items.push_back ({ 5, bounds[static_cast<std::size_t> (c)], bounds[static_cast<std::size_t> (c + 1)] });
            pmlSplit = domain.pmlPlaneSplit (threads);
        }
    }
    itemCount = static_cast<int> (items.size());
    // The output band ends at the air grid's on-axis cut-off (Decimator.h).
    const double stopHz = options.withAir ? std::asin (domain.courant()) / (3.14159265358979323846 * k)
                                          : 0.45 * sampleRate;
    for (auto& d : decimator)
        d.design (numerics.oversample, 1.0 / k, std::min (stopHz, 0.45 * sampleRate));
    hertz = HertzContact();
    step = 0;
    struckYet = false;
    prepared = true;
    prepareSeconds = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
    return true;
}

template <typename Real>
void FdtdEngine<Real>::strike (const StrokeSpec& stroke)
{
    strokeSpec = stroke;
    stickMass = stroke.tipMass;
    hertz.stiffness = stroke.hertzStiffness;
    hertz.exponent = stroke.hertzExponent;
    // The head grids must admit the SAV tension this stroke's energy allows
    // (FdtdDrumConfig.h headSpacing). The drum holds the stroke's kinetic
    // energy plus each head's psi0^2 / 2 = eps / 2; eps (1e-18 J) is left out
    // so the regulariser does not move the grids (test 3 halves it), which the
    // spacing's 1e-6 widening covers many times over. A drum still at rest has
    // its heads rebuilt; a moving one cannot be, and says so in diagnostics().
    const double energy = strokeEnergy (stroke);
    if (energy > headEnergy)
    {
        if (step == 0 && ! struckYet)
        {
            const auto t0 = std::chrono::steady_clock::now();
            headEnergy = energy;
            std::string error;
            buildHeads (drumSpec, error);
            if (! error.empty())
                fail ("head rebuild for the stroke's energy", error);
            buildHeadCoupling (drumSpec);
            for (auto& head : heads)
                head.scratch.assign (static_cast<std::size_t> (pool->size()), std::vector<double> (head.w[0].size(), 0.0));
            rebuildSeconds = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
        }
        else
            strokeBeyondBound = true;
    }
    // A batter whose stretching the energy does not bound gets its grid from
    // sizing runs of this stroke (FdtdEngineSizing.cpp), on a drum at rest.
    if (step == 0 && ! struckYet && needsBatterSizing() && failed.empty())
        sizeBatter (stroke);
    for (auto& head : heads)
    {
        head.savTension = head.peakSavTension = head.peakSavTensionTime = head.windowSavTension = 0.0;
        head.consistency = {};
    }
    precompute();
    // On the head surface, moving in at the impact speed (shared decision 3).
    const auto& head = heads[0];
    const double surface = head.grid.dot (head.footprint.data(), head.w[step % 3].data());
    stickY[step % 3] = surface;
    stickY[(step + 2) % 3] = surface + stroke.impactSpeed * k;
    stats = ContactStatistics();
    struckYet = true;
    outputPreroll = true;
    if (options.selfChecks && checkReport.empty() && failed.empty())
        runSelfChecks (drumSpec, micSpec, outputRate);
}

template <typename Real>
void FdtdEngine<Real>::setHeadAtRest (int index, const std::vector<double>& shape)
{
    auto& head = heads[index];
    for (std::size_t l = 0; l < shape.size(); ++l)
        head.w[mod3Of (step)][l] = head.w[mod3Of (step - 1)][l] = head.grid.isInterior (static_cast<int> (l)) ? shape[l] : 0.0;
    head.grid.laplacian (head.w[mod3Of (step - 1)].data(), head.lap[static_cast<std::size_t> ((step + 1) & 1)].data(), 0, head.grid.rows());
}

template class FdtdEngine<float>;
template class FdtdEngine<double>;

} // namespace taikor::candidates::fdtd
