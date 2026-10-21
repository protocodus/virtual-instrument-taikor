// Foppl-von Karman stretching potential with an immovable edge; the model,
// discretisation and references are in VonKarmanPotential.h.

#include "VonKarmanPotential.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>

namespace taikor::candidates
{

bool VonKarmanPotential::prepare (const MembraneGrid& grid, const HeadSpec& head,
                                  double headRadius, std::string& error)
{
    if (grid.nx != grid.ny || grid.nx < 8 || ! (grid.h > 0.0)
        || static_cast<int> (grid.interior.size()) != grid.size())
    {
        error = "von Karman potential: grid must be a square MembraneGrid";
        return false;
    }
    const double E = head.youngsModulus;
    const double thickness = head.thickness;
    const double nu = head.poissonRatio;
    if (! (E > 0.0) || ! (thickness > 0.0) || ! (nu > -1.0 && nu < 0.5) || ! (headRadius > 0.0))
    {
        error = "von Karman potential: head needs E > 0, h > 0, -1 < nu < 1/2 and a > 0";
        return false;
    }
    // headRadius is not used: the rim is the grid's staircase (interior
    // flags), the same boundary the linear head and the Berger letter see.
    grid_ = grid;
    tension_ = head.tension > 0.0 ? head.tension : 0.0;
    n_ = grid.nx;
    N_ = n_ * n_;
    h_ = grid.h;
    c_.c11 = E * thickness / (1.0 - nu * nu);
    c_.c12 = nu * c_.c11;
    c_.c66 = 0.5 * (1.0 - nu) * c_.c11;
    vk::InPlaneSolver::Settings settings = solver_.settings();
    if (! solver_.build (n_, grid.interior, c_, settings, error))
        return false;
    pad_ = solver_.fine().pad;
    const auto nodeSize = static_cast<std::size_t> (N_ + 2 * pad_);
    for (auto* v : { &w_, &ax_, &by_, &sx_, &sy_, &qxx_, &qyy_, &qxy_, &nxx_, &nyy_, &nxy_, &fx_, &fy_ })
        v->assign (nodeSize, 0.0);
    const auto levelSize = static_cast<std::size_t> (solver_.fine().size());
    rhs_.assign (levelSize, 0.0);
    u0_.assign (levelSize, 0.0);
    u1_.assign (levelSize, 0.0);
    history_ = 0;
    stats_ = {};
    return true;
}

void VonKarmanPotential::slopes (double* wp, double* ax, double* by, double* sx, double* sy) const
{
    slopesRows (wp, ax, by, 0, N_);
    cornerSlopesRows (ax, by, sx, sy, 0, N_);
}

void VonKarmanPotential::slopesRows (const double* wp, double* ax, double* by, int k0, int k1) const
{
    const int n = n_;
    const double invH = 1.0 / h_;
    for (int k = k0; k < k1; ++k)
    {
        ax[k] = (wp[k + 1] - wp[k]) * invH;
        by[k] = (wp[k + n] - wp[k]) * invH;
    }
}

void VonKarmanPotential::cornerSlopesRows (const double* ax, const double* by, double* sx, double* sy, int k0, int k1) const
{
    const int n = n_;
    for (int k = k0; k < k1; ++k)
    {
        sx[k] = 0.5 * (ax[k] + ax[k + n]);
        sy[k] = 0.5 * (by[k] + by[k + 1]);
    }
}

// Stages of computeStrainSources over the node range [k0, k1) (whole rows);
// a team (Team.h) meets between them, since each reads its neighbours' rows
// of the stage before.
void VonKarmanPotential::sourcesMask (const double* w, int k0, int k1)
{
    const std::uint8_t* interior = grid_.interior.data();
    double* wp = w_.data() + pad_;
    for (int k = k0; k < k1; ++k)
        wp[k] = interior[k] != 0 ? w[k] : 0.0;
}

double VonKarmanPotential::sourcesStrain (int k0, int k1)
{
    const int n = n_;
    double* ax = ax_.data() + pad_;
    double* by = by_.data() + pad_;
    double* sx = sx_.data() + pad_;
    double* sy = sy_.data() + pad_;
    double* qxx = qxx_.data() + pad_;
    double* qyy = qyy_.data() + pad_;
    double* qxy = qxy_.data() + pad_;
    cornerSlopesRows (ax, by, sx, sy, k0, k1);
    for (int k = k0; k < k1; ++k)
    {
        qxx[k] = 0.25 * (ax[k] * ax[k] + ax[k - 1] * ax[k - 1]);
        qyy[k] = 0.25 * (by[k] * by[k] + by[k - n] * by[k - n]);
        qxy[k] = sx[k] * sy[k];
    }
    // The in-plane stresses of the stretching strain with U = 0.
    double* sxx = nxx_.data() + pad_;
    double* syy = nyy_.data() + pad_;
    double* sxy = nxy_.data() + pad_;
    const double c11 = c_.c11, c12 = c_.c12, c66 = c_.c66;
    double peak = 0.0;
    for (int k = k0; k < k1; ++k)
    {
        sxx[k] = c11 * qxx[k] + c12 * qyy[k];
        syy[k] = c12 * qxx[k] + c11 * qyy[k];
        sxy[k] = c66 * qxy[k];
        peak = std::max (peak, sxx[k] + syy[k]);
    }
    return peak;
}

void VonKarmanPotential::sourcesRhs (int k0, int k1)
{
    // Right-hand side of K U = -h D^T C q (the in-plane forces of the
    // stretching strain with U = 0), on the free dofs only.
    const int n = n_;
    const double* sxx = nxx_.data() + pad_;
    const double* syy = nyy_.data() + pad_;
    const double* sxy = nxy_.data() + pad_;
    const vk::InPlaneLevel& L = solver_.fine();
    double* bu = rhs_.data() + L.base (0);
    double* bv = rhs_.data() + L.base (1);
    const std::uint8_t* fu = L.freeU.data();
    const std::uint8_t* fv = L.freeV.data();
    const double mh = -h_;
    for (int k = k0; k < k1; ++k)
    {
        const double gu = sxx[k] - sxx[k + 1] + sxy[k - n] - sxy[k];
        const double gv = syy[k] - syy[k + n] + sxy[k - 1] - sxy[k];
        bu[k] = fu[k] != 0 ? mh * gu : 0.0;
        bv[k] = fv[k] != 0 ? mh * gv : 0.0;
    }
}

void VonKarmanPotential::computeStrainSources (const double* w)
{
    const int N = N_;
    sourcesMask (w, 0, N);
    slopesRows (w_.data() + pad_, ax_.data() + pad_, by_.data() + pad_, 0, N);
    sourceStress_ = 0.5 * sourcesStrain (0, N);
    sourcesRhs (0, N);
}

double VonKarmanPotential::energyRows (const double* inPlane, int k0, int k1)
{
    const int n = n_;
    const double invH = 1.0 / h_;
    const vk::InPlaneLevel& L = solver_.fine();
    const double* u = inPlane + L.base (0);
    const double* v = inPlane + L.base (1);
    const double* qxx = qxx_.data() + pad_;
    const double* qyy = qyy_.data() + pad_;
    const double* qxy = qxy_.data() + pad_;
    double* sxx = nxx_.data() + pad_;
    double* syy = nyy_.data() + pad_;
    double* sxy = nxy_.data() + pad_;
    const double c11 = c_.c11, c12 = c_.c12, c66 = c_.c66;
    double sum = 0.0;
    for (int k = k0; k < k1; ++k)
    {
        const double exx = (u[k] - u[k - 1]) * invH + qxx[k];
        const double eyy = (v[k] - v[k - n]) * invH + qyy[k];
        const double gxy = (u[k + n] - u[k] + v[k + 1] - v[k]) * invH + qxy[k];
        sxx[k] = c11 * exx + c12 * eyy;
        syy[k] = c12 * exx + c11 * eyy;
        sxy[k] = c66 * gxy;
        sum += sxx[k] * exx + syy[k] * eyy + sxy[k] * gxy;
    }
    return sum;
}

double VonKarmanPotential::energyAndGradient (const double* inPlane, double* gradient)
{
    const double h2 = h_ * h_;
    const double energy = 0.5 * h2 * energyRows (inPlane, 0, N_);
    if (gradient == nullptr)
        return energy;

    stressDivergence (ax_.data() + pad_, by_.data() + pad_, sx_.data() + pad_, sy_.data() + pad_, gradient);
    return energy;
}

void VonKarmanPotential::stressDivergence (const double* ax, const double* by, const double* sx,
                                           const double* sy, double* gradient)
{
    divergenceFlux (ax, by, sx, sy, 0, N_);
    divergenceGradient (gradient, 0, N_);
}

void VonKarmanPotential::divergenceFlux (const double* ax, const double* by, const double* sx, const double* sy,
                                         int k0, int k1)
{
    // dV/d(slope) on every edge, h^2 N : dq/d(slope), then dV/dw by the
    // transpose of the forward differences (divergenceGradient).
    const int n = n_;
    const double* sxx = nxx_.data() + pad_;
    const double* syy = nyy_.data() + pad_;
    const double* sxy = nxy_.data() + pad_;
    double* fx = fx_.data() + pad_;
    double* fy = fy_.data() + pad_;
    const double half = 0.5 * h_ * h_;
    for (int k = k0; k < k1; ++k)
    {
        fx[k] = half * (ax[k] * (sxx[k] + sxx[k + 1]) + sxy[k] * sy[k] + sxy[k - n] * sy[k - n]);
        fy[k] = half * (by[k] * (syy[k] + syy[k + n]) + sxy[k] * sx[k] + sxy[k - 1] * sx[k - 1]);
    }
}

void VonKarmanPotential::divergenceGradient (double* gradient, int k0, int k1)
{
    const int n = n_;
    const double invH = 1.0 / h_;
    const double* fx = fx_.data() + pad_;
    const double* fy = fy_.data() + pad_;
    const std::uint8_t* interior = grid_.interior.data();
    for (int k = k0; k < k1; ++k)
    {
        const double g = (fx[k - 1] - fx[k] + fy[k - n] - fy[k]) * invH;
        gradient[k] = interior[k] != 0 ? g : 0.0;
    }
}

void VonKarmanPotential::applyStressOperator (const double* d, double* out)
{
    const auto nodeSize = w_.size();
    for (auto* v : { &dw_, &dax_, &dby_, &dsx_, &dsy_ })
        if (v->size() != nodeSize)
            v->assign (nodeSize, 0.0);
    double* dp = dw_.data() + pad_;
    const std::uint8_t* interior = grid_.interior.data();
    for (int k = 0; k < N_; ++k)
        dp[k] = interior[k] != 0 ? d[k] : 0.0;
    slopes (dp, dax_.data() + pad_, dby_.data() + pad_, dsx_.data() + pad_, dsy_.data() + pad_);
    stressDivergence (dax_.data() + pad_, dby_.data() + pad_, dsx_.data() + pad_, dsy_.data() + pad_, out);
}

double VonKarmanPotential::peakStretchingTension() const
{
    if (peakTension_ >= 0.0)
        return peakTension_;
    // Over the whole square (the stresses at the clamped nodes next to the
    // rim count too: the edges joining them to interior nodes carry slopes).
    const int n = n_, N = N_;
    const double* sxx = nxx_.data() + pad_;
    const double* syy = nyy_.data() + pad_;
    const double* sxy = nxy_.data() + pad_;
    double peak = 0.0;
    for (int k = n; k < N - n - 1; ++k)
    {
        const double xEdge = 0.5 * (sxx[k] + sxx[k + 1]) + 0.5 * (std::abs (sxy[k]) + std::abs (sxy[k - n]));
        const double yEdge = 0.5 * (syy[k] + syy[k + n]) + 0.5 * (std::abs (sxy[k]) + std::abs (sxy[k - 1]));
        peak = std::max (peak, std::max (xEdge, yEdge));
    }
    peakTension_ = peak;
    return peak;
}

double VonKarmanPotential::evaluate (const double* w, double* gradient)
{
    peakTension_ = -1.0;
    computeStrainSources (w);
    // Initial guess for the in-plane solve, written into u1_ (the older
    // history slot), then swapped to the front.
    const std::size_t size = u0_.size();
    const int order = std::min (warmOrder_, history_);
    if (order <= 0)
        std::fill (u1_.begin(), u1_.end(), 0.0);
    else if (order == 1)
        std::copy (u0_.begin(), u0_.end(), u1_.begin());
    else
        for (std::size_t k = 0; k < size; ++k)
            u1_[k] = 2.0 * u0_[k] - u1_[k];
    double residual = 0.0;
    double tolerance = solver_.settings().relativeTolerance;
    if (tensionRelative_ && sourceStress_ > 0.0)
        tolerance *= 1.0 + tension_ / sourceStress_;
    const int iterations = solver_.solve (rhs_.data(), u1_.data(), residual, tolerance);
    stats_.toleranceSum += tolerance;
    std::swap (u0_, u1_);
    history_ = std::min (history_ + 1, 2);

    ++stats_.calls;
    stats_.iterations += iterations;
    stats_.maxIterations = std::max (stats_.maxIterations, iterations);
    stats_.lastIterations = iterations;
    stats_.lastResidual = residual;
    stats_.maxResidual = std::max (stats_.maxResidual, residual);
    stats_.initialResidualSum += solver_.lastInitialResidual();
    stats_.initialResidualMax = std::max (stats_.initialResidualMax, solver_.lastInitialResidual());
    return energyAndGradient (u0_.data(), gradient);
}

int VonKarmanPotential::setTeam (int size)
{
    size = std::clamp (size, 1, vk::maxTeam);
    teamBarrier_ = size > 1 ? std::make_unique<fdtd::Barrier> (size) : nullptr;
    teamSize_ = size;
    for (int m = 0; m < vk::maxTeam; ++m)
    {
        auto& team = teams_[static_cast<std::size_t> (m)];
        team = {};
        team.member = m;
        team.size = size;
        team.barrier = teamBarrier_.get();
        team.shared = &teamShared_;
    }
    return size;
}

double VonKarmanPotential::evaluateMember (const double* w, double* gradient, int member)
{
    if (teamSize_ == 1)
        return evaluate (w, gradient);
    return evaluateTeam (w, gradient, teams_[static_cast<std::size_t> (member)]);
}

double VonKarmanPotential::evaluateTeam (const double* w, double* gradient, vk::Team& team)
{
    // The serial evaluate() stage by stage, each member on its rows; the team
    // meets where a stage reads rows another member wrote.
    int j0, j1;
    team.rows (n_, j0, j1);
    const int k0 = j0 * n_, k1 = j1 * n_;
    if (team.member == 0)
        peakTension_ = -1.0;
    sourcesMask (w, k0, k1);
    team.sync();
    slopesRows (w_.data() + pad_, ax_.data() + pad_, by_.data() + pad_, k0, k1);
    team.sync();
    const double peak = team.max (sourcesStrain (k0, k1)); // also the barrier before sourcesRhs
    sourcesRhs (k0, k1);
    // The warm start into u1_ (the older history slot): the members' rows of
    // both components (the pads stay zero).
    const vk::InPlaneLevel& L = solver_.fine();
    const int order = std::min (warmOrder_, history_);
    for (int c = 0; c < 2; ++c)
    {
        const long b0 = L.base (c) + k0, b1 = L.base (c) + k1;
        double* __restrict older = u1_.data();
        const double* __restrict newer = u0_.data();
        for (long k = b0; k < b1; ++k)
            older[k] = order <= 0 ? 0.0 : order == 1 ? newer[k] : 2.0 * newer[k] - older[k];
    }
    double tolerance = solver_.settings().relativeTolerance;
    if (tensionRelative_ && peak > 0.0)
        tolerance *= 1.0 + tension_ / (0.5 * peak);
    team.sync(); // rhs and the initial guess complete
    double residual = 0.0;
    double* solution = u1_.data();
    const int iterations = solver_.solve (rhs_.data(), solution, residual, tolerance, team);
    // The energy and its gradient at the solution (complete: solve() ends
    // with a barrier after its last update).
    const double energy = 0.5 * h_ * h_ * team.sum (energyRows (solution, k0, k1)); // barrier: N complete
    divergenceFlux (ax_.data() + pad_, by_.data() + pad_, sx_.data() + pad_, sy_.data() + pad_, k0, k1);
    team.sync();
    divergenceGradient (gradient, k0, k1);
    team.sync(); // every member done with u0_ / u1_ and the gradient complete
    if (team.member == 0)
    {
        sourceStress_ = 0.5 * peak;
        std::swap (u0_, u1_);
        history_ = std::min (history_ + 1, 2);
        stats_.toleranceSum += tolerance;
        ++stats_.calls;
        stats_.iterations += iterations;
        stats_.maxIterations = std::max (stats_.maxIterations, iterations);
        stats_.lastIterations = iterations;
        stats_.lastResidual = residual;
        stats_.maxResidual = std::max (stats_.maxResidual, residual);
        stats_.initialResidualSum += solver_.lastInitialResidual();
        stats_.initialResidualMax = std::max (stats_.initialResidualMax, solver_.lastInitialResidual());
    }
    return energy;
}

double VonKarmanPotential::energyAt (const double* w, const double* inPlane)
{
    peakTension_ = -1.0;
    computeStrainSources (w);
    return energyAndGradient (inPlane, nullptr);
}

std::string VonKarmanPotential::diagnostics() const
{
    long free = 0;
    const auto& L = solver_.fine();
    for (int k = 0; k < L.N; ++k)
        free += (L.freeU[static_cast<std::size_t> (k)] != 0) + (L.freeV[static_cast<std::size_t> (k)] != 0);
    const double calls = stats_.calls > 0 ? static_cast<double> (stats_.calls) : 1.0;
    char text[512];
    std::snprintf (text, sizeof (text),
                   "von Karman in-plane solve: %ld calls, %ld free in-plane dofs on %d multigrid levels, "
                   "PCG iterations mean %.2f (max %d), relative residual max %.1e (tolerance %.0e), warm-start residual mean %.2e",
                   stats_.calls, free, solver_.levelCount(), static_cast<double> (stats_.iterations) / calls,
                   stats_.maxIterations, stats_.maxResidual, solver_.settings().relativeTolerance,
                   stats_.initialResidualSum / calls);
    return text;
}

std::unique_ptr<NonlinearPotential> makeVonKarmanPotential()
{
    return std::make_unique<VonKarmanPotential>();
}

} // namespace taikor::candidates
