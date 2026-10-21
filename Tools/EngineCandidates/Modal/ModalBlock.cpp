// The coupled block of one angular order, and the lone high batter modes.
//
// SELF-CONSISTENT EXTERIOR LOADING. The exterior added mass M_ext(w) depends
// on frequency, so each coupled mode k must satisfy
//   [K - w_k^2 (M_0 + M_ext(w_k))] x_k = 0
// (the nonlinear eigenproblem Christian et al. 1984 and Tiwari & Gupta 2017
// solve by iterating each mode). Sorting the eigenvalues of the pencil
// (K, M_0 + M_ext(w)) makes each lambda_k(w) a continuous function of w
// (Courant-Fischer), so the k-th coupled mode is the fixed point
// w^2 = lambda_k(w). The pencil is solved on a log-spaced frequency grid
// (gridPerOctave points per octave); each fixed point is found by linear
// interpolation of sqrt(lambda_k(w)) - w between the two grid points that
// bracket it, the eigenvector is the Rayleigh-Ritz vector of the pencil with
// M_0 + M_ext interpolated to the fixed point, projected on the neighbouring
// eigenvectors of both bracketing grid points (ritzRefine), and the frequency
// is its Rayleigh quotient.
// The grid density is a numerical parameter (tests/modal build_grid_convergence).
//
// DAMPING. The damping operator C (mass-normalised x_k, decay-rate matrix
// A_kl = x_k^T C x_l / 2) collects
//  - hide and rim loss exactly as the shared finite-difference form
//    -2 sigma sigma0 w_t + 2 sigma sigma1 lap w_t, i.e. C_h = 2 sigma0 M_hide +
//    2 sigma sigma1 Gamma (Gamma = int grad phi_i . grad phi_j, ModalMath.h);
//  - radiation: the exterior resistance matrix R_ext(w_k) of the drum body
//    (both heads, ModalBodyExterior.h), per head baffled above its table;
//  - barrel air: the Kirchhoff viscothermal boundary-layer loss on the rigid
//    side wall (the shared decision; ModalCavity.h wallLoss), of the coupled
//    mode's whole wall field (its air coordinates and, for m = 0, the uniform
//    compression), at the coupled mode's own frequency.
// Each mode's own decay is the first-order projection alpha_k = A_kk. The
// couplings A_kl between modes closer than dampingBandwidth (alpha_k + alpha_l)
// are kept as well (see the end of solveCoupledBlock): a hide-damped head mode
// and a lightly damped air mode that the conservative problem mixes 50/50 do
// not hybridise once their decay rates differ by more than their coupling, and
// only the dense damping of their cluster reproduces that (tests/modal
// damping_cluster_pair against the exact damped two-degree-of-freedom solution).

#include "ModalEigen.h"
#include "ModalInternal.h"

#include <algorithm>
#include <cmath>
#include <functional>

namespace taikor::candidates::modal
{
namespace
{
double quadratic (const Matrix& m, const std::vector<double>& x)
{
    double sum = 0.0;
    for (int i = 0; i < m.n; ++i)
    {
        double row = 0.0;
        for (int j = 0; j < m.n; ++j)
            row += m (i, j) * x[static_cast<std::size_t> (j)];
        sum += x[static_cast<std::size_t> (i)] * row;
    }
    return sum;
}

// Rayleigh-Ritz refinement of a fixed point's eigenvector: the pencil
// (K, M_star) at the interpolated mass is projected on the eigenvectors k-1,
// k, k+1 of the two bracketing grid points (they span the vector's first-order
// variation across the grid step), and the Ritz vector that overlaps most with
// the nearer grid vector replaces it. The error in the eigenvector, and hence
// in the damping projection, becomes second order in the grid step
// (tests/modal build_grid_convergence).
void ritzRefine (const Matrix& k, const Matrix& mass, const Matrix& left, const Matrix& right, int index,
                 std::vector<double>& x)
{
    const int n = k.n;
    std::vector<std::vector<double>> basis;
    const auto mdot = [&] (const std::vector<double>& p, const std::vector<double>& q)
    {
        double sum = 0.0;
        for (int i = 0; i < n; ++i)
        {
            double row = 0.0;
            for (int j = 0; j < n; ++j)
                row += mass (i, j) * q[static_cast<std::size_t> (j)];
            sum += p[static_cast<std::size_t> (i)] * row;
        }
        return sum;
    };
    for (const Matrix* src : { &left, &right })
        for (int c = std::max (0, index - 1); c <= std::min (n - 1, index + 1); ++c)
        {
            std::vector<double> v (static_cast<std::size_t> (n));
            for (int i = 0; i < n; ++i)
                v[static_cast<std::size_t> (i)] = (*src) (i, c);
            const double original = std::sqrt (mdot (v, v));
            for (int pass = 0; pass < 2; ++pass)
                for (const auto& u : basis)
                {
                    const double proj = mdot (u, v);
                    for (int i = 0; i < n; ++i)
                        v[static_cast<std::size_t> (i)] -= proj * u[static_cast<std::size_t> (i)];
                }
            const double norm = std::sqrt (std::max (mdot (v, v), 0.0));
            if (norm <= 1.0e-8 * original)
                continue;
            for (auto& value : v)
                value /= norm;
            basis.push_back (std::move (v));
        }
    const int p = static_cast<int> (basis.size());
    if (p < 2)
        return;
    Matrix small (p);
    for (int i = 0; i < p; ++i)
        for (int j = i; j < p; ++j)
        {
            double sum = 0.0;
            for (int r = 0; r < n; ++r)
            {
                double row = 0.0;
                for (int s = 0; s < n; ++s)
                    row += k (r, s) * basis[static_cast<std::size_t> (j)][static_cast<std::size_t> (s)];
                sum += basis[static_cast<std::size_t> (i)][static_cast<std::size_t> (r)] * row;
            }
            small (i, j) = small (j, i) = sum;
        }
    std::vector<double> theta;
    Matrix z;
    if (! symmetricEigen (small, theta, z))
        return;
    std::vector<double> overlap (static_cast<std::size_t> (p));
    for (int i = 0; i < p; ++i)
        overlap[static_cast<std::size_t> (i)] = mdot (basis[static_cast<std::size_t> (i)], x);
    int best = 0;
    double bestValue = -1.0;
    for (int r = 0; r < p; ++r)
    {
        double dot = 0.0;
        for (int i = 0; i < p; ++i)
            dot += z (i, r) * overlap[static_cast<std::size_t> (i)];
        if (std::abs (dot) > bestValue)
        {
            bestValue = std::abs (dot);
            best = r;
        }
    }
    std::fill (x.begin(), x.end(), 0.0);
    for (int i = 0; i < p; ++i)
        for (int r = 0; r < n; ++r)
            x[static_cast<std::size_t> (r)] += z (i, best) * basis[static_cast<std::size_t> (i)][static_cast<std::size_t> (r)];
}
} // namespace

void solveCoupledBlock (const DrumContext& ctx, int m, std::vector<CoupledMode>& out)
{
    const BuildOptions& o = ctx.options;
    const auto um = static_cast<std::size_t> (m);
    std::vector<const HeadMode*> bModes, rModes;
    if (um < ctx.batter.size())
        for (const auto& mode : ctx.batter[um])
            if (mode.omega < ctx.omegaCoupled)
                bModes.push_back (&mode);
    if (um < ctx.rear.size())
        for (const auto& mode : ctx.rear[um])
            if (mode.omega < ctx.omegaCoupled)
                rModes.push_back (&mode);
    std::vector<AirMode> air;
    const std::vector<CrossMode>* cross = um < ctx.cross.size() ? &ctx.cross[um] : nullptr;
    if (o.interior && cross != nullptr)
        air = airModes (m, *cross, ctx.air, ctx.omegaAir);
    const int nb = static_cast<int> (bModes.size()), nr = static_cast<int> (rModes.size());
    const int nh = nb + nr, na = static_cast<int> (air.size()), n = nh + na;
    if (nh == 0)
        return; // air modes with no head to couple to are never excited

    const double sb = ctx.batterLoss.arealDensity, sr = ctx.rearLoss.arealDensity;
    const auto hideMass = [&] (int i) { return i < nb ? sb * bModes[static_cast<std::size_t> (i)]->normSq
                                                      : sr * rModes[static_cast<std::size_t> (i - nb)]->normSq; };
    Matrix k (n), m0 (n);
    for (int i = 0; i < nh; ++i)
    {
        const double w = i < nb ? bModes[static_cast<std::size_t> (i)]->omega : rModes[static_cast<std::size_t> (i - nb)]->omega;
        k (i, i) = w * w * hideMass (i);
        m0 (i, i) = hideMass (i);
    }
    // m = 0: the volume each head mode displaces outward (Lam_0 W(0)), for the
    // barrel's uniform compression pressure in the wall loss (ModalCavity.h).
    std::vector<double> volume;
    if (o.interior && cross != nullptr)
    {
        const double a = ctx.air.radius;
        std::vector<std::vector<double>> bh (static_cast<std::size_t> (nb)), rh (static_cast<std::size_t> (nr));
        for (int i = 0; i < nb; ++i)
            for (const auto& q : *cross)
                bh[static_cast<std::size_t> (i)].push_back (hankelTransform (*bModes[static_cast<std::size_t> (i)], a, q.gamma, q.jmEdge, 0.0));
        for (int j = 0; j < nr; ++j)
            for (const auto& q : *cross)
                rh[static_cast<std::size_t> (j)].push_back (hankelTransform (*rModes[static_cast<std::size_t> (j)], a, q.gamma, q.jmEdge, 0.0));
        if (m == 0 && ! cross->empty() && (*cross)[0].gamma == 0.0)
        {
            for (int i = 0; i < nb; ++i)
                volume.push_back (angularNorm (0) * bh[static_cast<std::size_t> (i)][0]);
            for (int j = 0; j < nr; ++j)
                volume.push_back (angularNorm (0) * rh[static_cast<std::size_t> (j)][0]);
        }
        const InteriorBlock ib = interiorBlock (m, *cross, air, ctx.air, bh, rh);
        for (int i = 0; i < nh; ++i)
            for (int j = 0; j < nh; ++j)
            {
                k (i, j) += ib.staticStiffness[static_cast<std::size_t> (i * nh + j)];
                m0 (i, j) += ib.residualMass[static_cast<std::size_t> (i * nh + j)];
            }
        for (int s = 0; s < na; ++s)
        {
            for (int i = 0; i < nh; ++i)
            {
                const double c = ib.coupling[static_cast<std::size_t> (i) * static_cast<std::size_t> (na) + static_cast<std::size_t> (s)];
                k (i, nh + s) = c;
                k (nh + s, i) = c;
            }
            k (nh + s, nh + s) = air[static_cast<std::size_t> (s)].omega * air[static_cast<std::size_t> (s)].omega;
            m0 (nh + s, nh + s) = 1.0;
        }
    }

    // Exterior added mass at w: the per-head baffled King form (ModalExterior.h)
    // for every pair, with the sub-block of the head modes below the body table's
    // top replaced, at its grid points, by the drum body's boundary-element
    // solution (both heads in one field, batter-rear terms included;
    // ModalBodyExterior.h). Table rows: the first nbc batter then the first nrc
    // rear modes of this block (both lists ascend in frequency).
    const BodyExterior* body = ctx.body.get();
    const int nbc = body != nullptr ? body->batterColumns (m) : 0;
    const int nrc = body != nullptr ? body->rearColumns (m) : 0;
    const auto blockIndex = [&] (int t) { return t < nbc ? t : nb + (t - nbc); };
    const auto overlay = [&] (const std::vector<double>& sub, Matrix& target)
    {
        const int nc = nbc + nrc;
        for (int i = 0; i < nc; ++i)
            for (int j = 0; j < nc; ++j)
                target (blockIndex (i), blockIndex (j)) = sub[static_cast<std::size_t> (i * nc + j)];
    };
    const auto exteriorAt = [&] (int g, double omega, Matrix& ext)
    {
        ext = Matrix (nh);
        if (! o.exterior)
            return;
        std::vector<double> mb, mr;
        const double kw = omega / ctx.air.soundSpeed;
        ctx.exterior->addedMass (bModes, kw, ctx.air.density, mb);
        ctx.exterior->addedMass (rModes, kw, ctx.air.density, mr);
        for (int i = 0; i < nb; ++i)
            for (int j = 0; j < nb; ++j)
                ext (i, j) = mb[static_cast<std::size_t> (i * nb + j)];
        for (int i = 0; i < nr; ++i)
            for (int j = 0; j < nr; ++j)
                ext (nb + i, nb + j) = mr[static_cast<std::size_t> (i * nr + j)];
        if (body != nullptr && nbc + nrc > 0)
        {
            const auto& z = body->impedance (m, g);
            if (z.size() == static_cast<std::size_t> ((nbc + nrc) * (nbc + nrc)))
            {
                std::vector<double> sub (z.size());
                for (std::size_t e = 0; e < z.size(); ++e)
                    sub[e] = z[e].imag() / omega;
                overlay (sub, ext);
            }
        }
    };
    const auto massAt = [&] (const Matrix& ext)
    {
        Matrix mm = m0;
        for (int i = 0; i < nh; ++i)
            for (int j = 0; j < nh; ++j)
                mm (i, j) += ext (i, j);
        return mm;
    };

    // A fixed point, normalised against the pencil's mass at that frequency.
    struct Found
    {
        double omega = 0.0, shareBatter = 0.0, shareRear = 0.0, shareAir = 0.0;
        std::vector<double> x;
    };
    std::vector<Found> found (static_cast<std::size_t> (n));
    std::vector<char> done (static_cast<std::size_t> (n), 0);
    std::vector<double> mx (static_cast<std::size_t> (n));
    const auto finish = [&] (Found& f, const Matrix& mass)
    {
        for (int i = 0; i < n; ++i)
        {
            double row = 0.0;
            for (int j = 0; j < n; ++j)
                row += mass (i, j) * f.x[static_cast<std::size_t> (j)];
            mx[static_cast<std::size_t> (i)] = row;
        }
        double norm = 0.0;
        for (int i = 0; i < n; ++i)
            norm += f.x[static_cast<std::size_t> (i)] * mx[static_cast<std::size_t> (i)];
        if (! (norm > 0.0))
            return false;
        const double scale = 1.0 / std::sqrt (norm);
        for (int i = 0; i < n; ++i)
        {
            f.x[static_cast<std::size_t> (i)] *= scale;
            const double e = f.x[static_cast<std::size_t> (i)] * mx[static_cast<std::size_t> (i)] * scale;
            (i < nb ? f.shareBatter : i < nh ? f.shareRear : f.shareAir) += e;
        }
        f.omega = std::sqrt (std::max (quadratic (k, f.x), 0.0));
        return true;
    };

    if (! o.exterior)
    {
        std::vector<double> values;
        Matrix vectors;
        generalizedEigen (k, m0, values, vectors);
        for (int c = 0; c < n; ++c)
        {
            auto& f = found[static_cast<std::size_t> (c)];
            f.x.resize (static_cast<std::size_t> (n));
            for (int i = 0; i < n; ++i)
                f.x[static_cast<std::size_t> (i)] = vectors (i, c);
            done[static_cast<std::size_t> (c)] = finish (f, m0) ? 1 : 0;
        }
    }
    else
    {
        double wLo = 1e300, wHi = 0.0;
        for (const auto* mode : bModes) { wLo = std::min (wLo, mode->omega); wHi = std::max (wHi, mode->omega); }
        for (const auto* mode : rModes) { wLo = std::min (wLo, mode->omega); wHi = std::max (wHi, mode->omega); }
        for (const auto& am : air) { wLo = std::min (wLo, am.omega); wHi = std::max (wHi, am.omega); }
        // Air loading lowers a head mode by at most ~35 % here ((0,1): 80 -> 55 Hz),
        // so half the lowest uncoupled frequency brackets every fixed point from below.
        wLo *= 0.5;
        wHi *= 1.3;
        // The global grid w_g = 2 pi 2^(g/G) (shared with the body table) over [wLo, wHi].
        const int gLo = static_cast<int> (std::floor (o.gridPerOctave * std::log2 (wLo / (2.0 * pi))));
        const int gHi = static_cast<int> (std::ceil (o.gridPerOctave * std::log2 (wHi / (2.0 * pi))));
        const int points = std::max (2, gHi - gLo + 1);
        std::vector<double> prevValues, values;
        Matrix prevVectors, vectors, prevExt, ext;
        double prevOmega = 0.0;
        for (int g = 0; g < points; ++g)
        {
            const double omega = bodyGridOmega (gLo + g, o.gridPerOctave);
            exteriorAt (gLo + g, omega, ext);
            if (! generalizedEigen (k, massAt (ext), values, vectors))
                continue;
            for (int c = 0; c < n; ++c)
            {
                if (done[static_cast<std::size_t> (c)])
                    continue;
                const double h = std::sqrt (std::max (values[static_cast<std::size_t> (c)], 0.0)) - omega;
                const bool last = g == points - 1;
                if (h > 0.0 && ! last)
                    continue;
                auto& f = found[static_cast<std::size_t> (c)];
                f.x.resize (static_cast<std::size_t> (n));
                if (g == 0 || (h > 0.0 && last))
                {
                    // Below the first or above the last grid point: that point's pencil.
                    for (int i = 0; i < n; ++i)
                        f.x[static_cast<std::size_t> (i)] = vectors (i, c);
                    done[static_cast<std::size_t> (c)] = finish (f, massAt (ext)) ? 1 : 2;
                }
                else
                {
                    const double hp = std::sqrt (std::max (prevValues[static_cast<std::size_t> (c)], 0.0)) - prevOmega;
                    const double star = prevOmega + hp * (omega - prevOmega) / (hp - h);
                    const double t = std::log (star / prevOmega) / std::log (omega / prevOmega);
                    const Matrix& src = t < 0.5 ? prevVectors : vectors;
                    for (int i = 0; i < n; ++i)
                        f.x[static_cast<std::size_t> (i)] = src (i, c);
                    Matrix mix (nh);
                    for (std::size_t e = 0; e < mix.a.size(); ++e)
                        mix.a[e] = (1.0 - t) * prevExt.a[e] + t * ext.a[e];
                    const Matrix mass = massAt (mix);
                    ritzRefine (k, mass, prevVectors, vectors, c, f.x);
                    done[static_cast<std::size_t> (c)] = finish (f, mass) ? 1 : 2;
                }
            }
            prevValues = values;
            prevVectors = vectors;
            prevExt = ext;
            prevOmega = omega;
        }
    }

    const auto& bForm = ctx.batterGradient[um];
    const GradientForm* rForm = um < ctx.rearGradient.size() ? &ctx.rearGradient[um] : nullptr;
    const std::size_t first = out.size();
    std::vector<std::vector<double>> blockX, blockCx, blockHalf;
    for (int c = 0; c < n; ++c)
    {
        if (done[static_cast<std::size_t> (c)] != 1)
            continue;
        const auto& f = found[static_cast<std::size_t> (c)];
        CoupledMode mode;
        mode.m = m;
        mode.omega = f.omega;
        mode.shareBatter = f.shareBatter;
        mode.shareRear = f.shareRear;
        mode.shareAir = f.shareAir;
        for (int i = 0; i < nb; ++i)
        {
            mode.batterIndex.push_back (i);
            mode.batterCoef.push_back (f.x[static_cast<std::size_t> (i)]);
        }
        for (int j = 0; j < nr; ++j)
        {
            mode.rearIndex.push_back (j);
            mode.rearCoef.push_back (f.x[static_cast<std::size_t> (nb + j)]);
        }
        // The damping operator applied to the mode, C x (C = 2 sigma0 M_hide +
        // 2 sigma sigma1 Gamma + R_ext(w_k) + the side-wall loss), so that
        // alpha_k = x^T C x / 2 and the couplings x_l^T C x_k / 2 follow; cx holds
        // the hide part, radiation and wall loss enter as Gram half-vectors.
        std::vector<double> cx (static_cast<std::size_t> (n), 0.0);
        const auto applyGradient = [&] (const GradientForm& form, int offset, int count, double scale)
        {
            for (int i = 0; i < count; ++i)
                cx[static_cast<std::size_t> (offset + i)] += scale * form.diagonal[static_cast<std::size_t> (i)] * f.x[static_cast<std::size_t> (offset + i)];
            for (const auto& l : form.lowRank)
            {
                double dot = 0.0;
                for (int i = 0; i < count; ++i)
                    dot += l[static_cast<std::size_t> (i)] * f.x[static_cast<std::size_t> (offset + i)];
                for (int i = 0; i < count; ++i)
                    cx[static_cast<std::size_t> (offset + i)] -= scale * l[static_cast<std::size_t> (i)] * dot;
            }
        };
        for (int i = 0; i < nh; ++i)
            cx[static_cast<std::size_t> (i)] += 2.0 * (i < nb ? ctx.batterLoss.sigma0 : ctx.rearLoss.sigma0) * hideMass (i) * f.x[static_cast<std::size_t> (i)];
        applyGradient (bForm, 0, nb, 2.0 * sb * ctx.batterLoss.sigma1);
        if (nr > 0 && rForm != nullptr)
            applyGradient (*rForm, nb, nr, 2.0 * sr * ctx.rearLoss.sigma1);
        double hide = 0.0;
        for (int i = 0; i < nh; ++i)
            hide += 0.5 * f.x[static_cast<std::size_t> (i)] * cx[static_cast<std::size_t> (i)];
        mode.alphaHide = hide;
        // Radiation resistance at the mode's own frequency.
        if (o.exterior)
        {
            const double kw = mode.omega / ctx.air.soundSpeed;
            std::vector<double> rb, rr, rBody;
            ctx.exterior->resistance (bModes, kw, ctx.air.density, ctx.air.soundSpeed, rb);
            ctx.exterior->resistance (rModes, kw, ctx.air.density, ctx.air.soundSpeed, rr);
            std::vector<double> sub;
            if (body != nullptr && nbc + nrc > 0)
                body->resistance (m, mode.omega, sub);
            if (sub.size() == static_cast<std::size_t> ((nbc + nrc) * (nbc + nrc)))
            {
                // The full matrix: King per head, the body's sub-block overlaid.
                Matrix full (nh);
                for (int i = 0; i < nb; ++i)
                    for (int j = 0; j < nb; ++j)
                        full (i, j) = rb[static_cast<std::size_t> (i * nb + j)];
                for (int i = 0; i < nr; ++i)
                    for (int j = 0; j < nr; ++j)
                        full (nb + i, nb + j) = rr[static_cast<std::size_t> (i * nr + j)];
                overlay (sub, full);
                rBody = full.a;
            }
            // R(w_k)^{1/2} x per head: the radiation part of the couplings is the Gram
            // matrix of these vectors, positive semi-definite by construction although
            // each mode's resistance is taken at its own frequency (ModalBlock.cpp end).
            std::vector<double> half;
            double r = 0.0;
            const auto rootTimes = [&] (const std::vector<double>& res, int count, int offset)
            {
                Matrix rm (count);
                rm.a = res;
                std::vector<double> values;
                Matrix vectors;
                symmetricEigen (rm, values, vectors);
                std::vector<double> proj (static_cast<std::size_t> (count), 0.0);
                for (int e = 0; e < count; ++e)
                {
                    double dot = 0.0;
                    for (int i = 0; i < count; ++i)
                        dot += vectors (i, e) * f.x[static_cast<std::size_t> (offset + i)];
                    proj[static_cast<std::size_t> (e)] = std::sqrt (std::max (values[static_cast<std::size_t> (e)], 0.0)) * dot;
                }
                for (int i = 0; i < count; ++i)
                {
                    double v = 0.0;
                    for (int e = 0; e < count; ++e)
                        v += vectors (i, e) * proj[static_cast<std::size_t> (e)];
                    half.push_back (v);
                    r += v * v;
                }
            };
            if (! rBody.empty())
                rootTimes (rBody, nh, 0);   // both heads in one field below the table's top
            else
            {
                rootTimes (rb, nb, 0);
                if (nr > 0)
                    rootTimes (rr, nr, nb);
            }
            mode.alphaRadiation = 0.5 * r;
            blockHalf.push_back (std::move (half));
        }
        else
            blockHalf.emplace_back (static_cast<std::size_t> (nh), 0.0);
        // Kirchhoff side-wall loss of the mode's whole wall field (its air
        // coordinates and, for m = 0, the barrel's uniform compression) at the
        // mode's own frequency (ModalCavity.h wallLoss), in the same Gram form.
        if (cross != nullptr && (na > 0 || ! volume.empty()))
        {
            double compression = 0.0;
            if (! volume.empty())
            {
                double displaced = 0.0;
                for (int i = 0; i < nh; ++i)
                    displaced += volume[static_cast<std::size_t> (i)] * f.x[static_cast<std::size_t> (i)];
                const double rhoC2 = ctx.air.density * ctx.air.soundSpeed * ctx.air.soundSpeed;
                const double barrel = pi * ctx.air.radius * ctx.air.radius * ctx.air.length;
                compression = -std::sqrt (rhoC2) * displaced / barrel;   // P_0 / sqrt(rho c^2)
            }
            mode.alphaWall = wallLoss (m, *cross, air, ctx.air, f.x.data() + nh, compression, mode.omega, &blockHalf.back());
        }
        blockX.push_back (f.x);
        blockCx.push_back (std::move (cx));
        mode.alpha = mode.alphaHide + mode.alphaRadiation + mode.alphaWall;
        out.push_back (std::move (mode));
    }
    // Off-diagonal decay-rate couplings A_kl = x_l^T C x_k / 2 between close
    // modes (radiation and wall loss, both taken at each mode's own frequency:
    // (C(w_k)^{1/2} x_k) . (C(w_l)^{1/2} x_l) / 2, so the matrix is a Gram matrix
    // with each diagonal the mode's own alpha). First-order projection keeps only A_kk. In the
    // quadratic eigenproblem (lambda^2 + lambda Ct + W^2) v = 0 (Ct = 2A) a
    // coupling to a mode Delta w away shifts lambda_k by
    //   sum_l j w_k Ct_kl^2 / (2 (w_l^2 - w_k^2 + j w_k Ct_ll)):
    // a frequency shift of order Ct_kl^2 / (4 Delta w) and a decay change of
    // order Ct_kl^2 Ct_ll / (8 Delta w^2), which with |Ct_kl| <= 2 sqrt(A_kk A_ll)
    // is below 1 / (4 B^2) of the decay once Delta w > B (A_kk + A_ll). Pairs
    // closer than dampingBandwidth B are kept (when |A_kl| exceeds dampingFloor
    // sqrt(A_kk A_ll)) and stepped as dense damping in their cluster
    // (ModalRuntime.h): there the damping can decouple modes that the
    // conservative problem mixes. Close pairs also see one consistent C.
    const std::size_t count = out.size() - first;
    const auto coupling = [&] (std::size_t k, std::size_t l)
    {
        // The hide operator is frequency independent (symmetric); radiation and
        // wall loss enter through the Gram of C(w)^{1/2} x.
        double s1 = 0.0;
        for (int i = 0; i < n; ++i)
            s1 += blockX[l][static_cast<std::size_t> (i)] * blockCx[k][static_cast<std::size_t> (i)];
        double s2 = 0.0;
        for (std::size_t i = 0; i < blockHalf[k].size() && i < blockHalf[l].size(); ++i)
            s2 += blockHalf[k][i] * blockHalf[l][i];
        return 0.5 * (s1 + s2);
    };
    // Clusters: connected components of the close, non-negligible pairs. Every
    // pair inside a cluster then keeps its coupling, so each cluster's decay
    // matrix is a principal submatrix of the block's Gram matrix: positive
    // semi-definite, and so is any sub-cluster left after residue pruning.
    std::vector<std::size_t> parent (count);
    for (std::size_t k = 0; k < count; ++k)
        parent[k] = k;
    const std::function<std::size_t (std::size_t)> find = [&] (std::size_t x) { return parent[x] == x ? x : parent[x] = find (parent[x]); };
    for (std::size_t k = 0; k < count; ++k)
        for (std::size_t l = k + 1; l < count; ++l)
        {
            const CoupledMode& mk = out[first + k];
            const CoupledMode& ml = out[first + l];
            if (std::abs (mk.omega - ml.omega) > o.dampingBandwidth * (mk.alpha + ml.alpha))
                continue;
            if (std::abs (coupling (k, l)) >= o.dampingFloor * std::sqrt (mk.alpha * ml.alpha))
                parent[find (k)] = find (l);
        }
    std::vector<std::vector<std::size_t>> groups (count);
    for (std::size_t k = 0; k < count; ++k)
        groups[find (k)].push_back (k);
    for (const auto& members : groups)
        for (std::size_t i = 0; i < members.size(); ++i)
            for (std::size_t j = i + 1; j < members.size(); ++j)
            {
                const double a = coupling (members[i], members[j]);
                out[first + members[i]].damping.emplace_back (static_cast<int> (members[j]), a);
                out[first + members[j]].damping.emplace_back (static_cast<int> (members[i]), a);
            }
}

void solveHighModes (const DrumContext& ctx, int m, std::vector<CoupledMode>& out)
{
    const auto um = static_cast<std::size_t> (m);
    if (um >= ctx.batter.size())
        return;
    const auto& list = ctx.batter[um];
    const auto& form = ctx.batterGradient[um];
    const double a = ctx.air.radius, rho = ctx.air.density, c = ctx.air.soundSpeed;
    for (std::size_t i = 0; i < list.size(); ++i)
    {
        const HeadMode& hm = list[i];
        if (hm.omega < ctx.omegaCoupled)
            continue;
        const double hide = ctx.batterLoss.arealDensity * hm.normSq;
        double omega = hm.omega, extra = 0.0;
        for (int iteration = 0; iteration < 8; ++iteration)
        {
            const double kw = omega / c;
            extra = 0.0;
            if (ctx.options.exterior)
                extra += ctx.exterior->selfAddedMass (hm, kw, rho);
            if (ctx.options.interior && um < ctx.cross.size())
                extra += evanescentInteriorMass (hm, a, kw, ctx.air.length, rho, ctx.cross[um]);
            const double next = hm.omega * std::sqrt (hide / (hide + extra));
            const bool converged = std::abs (next - omega) < 1e-12 * omega;
            omega = next;
            if (converged)
                break;
        }
        const double total = hide + extra;
        CoupledMode mode;
        mode.m = m;
        mode.high = true;
        mode.omega = omega;
        mode.batterIndex.push_back (static_cast<int> (i));
        mode.batterCoef.push_back (1.0 / std::sqrt (total));
        mode.shareBatter = 1.0;
        double gamma = form.diagonal[i];
        for (const auto& l : form.lowRank)
            gamma -= l[i] * l[i];
        mode.alphaHide = (ctx.batterLoss.sigma0 * hide + ctx.batterLoss.arealDensity * ctx.batterLoss.sigma1 * gamma) / total;
        if (ctx.options.exterior)
            mode.alphaRadiation = 0.5 * ctx.exterior->selfResistance (hm, omega / c, rho, c) / total;
        mode.alpha = mode.alphaHide + mode.alphaRadiation;
        out.push_back (std::move (mode));
    }
}

} // namespace taikor::candidates::modal
