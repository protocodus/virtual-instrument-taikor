// Microphone transfer of every coupled mode, evaluated at its complex pole.
//
// BATTER. The batter radiates into the half-space z > 0 as a baffled plane
// (Rayleigh's first integral): p(x, t) = rho int w_tt(y, t - R/c) / (2 pi R) dS.
// A free-decaying mode w = Re{A e^{s t}} phi(y), s = -alpha + j omega, gives
// exactly p = Re{rho s^2 G(s) A e^{s t}} with
//   G(s) = int phi(y) e^{-s R/c} / (2 pi R) dS,
// the Laplace residue of the pole (the e^{alpha R/c} factor is the earlier,
// larger amplitude the listener hears). The integral keeps the near field
// (the evanescent content of high orders decays as e^{-b1 z} at the capsule).
// For phi(r) cos(m (theta - theta_s)) and a capsule at (rho_m, theta_m, z):
//   G = cos(m (theta_m - theta_s)) int_0^a phi(r) r A_m(r) dr,
//   A_m(r) = 1/(2 pi) int_0^{2 pi} cos(m p) e^{-s R(p)/c} / R(p) dp,
//   R(p) = sqrt(r^2 + rho_m^2 - 2 r rho_m cos p + z^2).
// Quadrature (numerical, tests/modal mics_convergence): Gauss-Legendre in r
// with 30 + ceil((b1 + |s|/c) a) nodes, and the periodic trapezoid in p with
// 64 + m + ceil(|s|/c r rho_m / z) points (spectrally accurate for a smooth
// periodic integrand), using its symmetry in p.
//
// REAR HEAD. Only its net volume displacement V reaches the batter-side
// capsules in this model: a monopole p = rho s^2 V e^{-s R'/c} / (4 pi R') at
// the rear head's centre, R' the straight distance to the capsule. This is the
// leading term in ka of the field of a compact source of net volume velocity
// on a rigid compact body: the body alters it only through its scattered
// wave, of relative order ka and above (A. D. Pierce, Acoustics, 3rd ed.,
// Springer/ASA 2019, sec. 4.3 (monopole) and sec. 9.1 (scattering by small
// objects)). Above ka ~ 1 the barrel's shadowing and the rear head's higher
// multipoles are not represented (a stated approximation; only m = 0 rear
// components have a volume displacement at all).

#include "ModalInternal.h"

#include <cmath>

namespace taikor::candidates::modal
{

void rayleighNodes (int m, double radius, const std::function<double (double)>& profile, double maxWavenumber,
                    double kAbs, double rhoMic, double zMic, int refine,
                    const std::function<void (double weight, double distance)>& visit)
{
    const double a = radius;
    const int nr = refine * (30 + static_cast<int> (std::ceil ((maxWavenumber + kAbs) * a)));
    std::vector<double> r, w;
    special::gaussLegendre (nr, r, w);
    const double zEff = std::max (std::abs (zMic), 1e-3);
    for (int n = 0; n < nr; ++n)
    {
        const double rr = 0.5 * a * (r[static_cast<std::size_t> (n)] + 1.0);
        const double wr = 0.5 * a * w[static_cast<std::size_t> (n)];
        const int half = refine * (32 + (m + static_cast<int> (std::ceil (kAbs * rr * rhoMic / zEff))) / 2 + 1);
        const int count = 2 * half;
        // A_m(r) = (1/2 pi) int cos(m p) e^{-sR/c}/R dp by the periodic trapezoid,
        // using the integrand's symmetry in p.
        const double radial = wr * profile (rr) * rr / static_cast<double> (count);
        for (int j = 0; j <= half; ++j)
        {
            const double phi = 2.0 * pi * j / count;
            const double dist = std::sqrt (rr * rr + rhoMic * rhoMic - 2.0 * rr * rhoMic * std::cos (phi) + zMic * zMic);
            const double multiplicity = (j == 0 || j == half) ? 1.0 : 2.0;
            visit (radial * multiplicity * std::cos (m * phi) / dist, dist);
        }
    }
}

std::complex<double> rayleighIntegral (int m, double radius, const std::function<double (double)>& profile,
                                       double maxWavenumber, std::complex<double> s, double soundSpeed,
                                       double rhoMic, double zMic, int refine)
{
    std::complex<double> g = 0.0;
    rayleighNodes (m, radius, profile, maxWavenumber, std::abs (s) / soundSpeed, rhoMic, zMic, refine,
                   [&] (double weight, double dist) { g += weight * std::exp (-s * (dist / soundSpeed)); });
    return g;
}

double rearDistance (double length, const std::array<double, 3>& capsule)
{
    return std::sqrt (capsule[0] * capsule[0] + capsule[1] * capsule[1] + (capsule[2] + length) * (capsule[2] + length));
}

std::function<double (double)> batterProfile (const CoupledMode& mode, const std::vector<std::vector<HeadMode>>& batter,
                                              double radius, double& maxWavenumber)
{
    const auto& list = batter[static_cast<std::size_t> (mode.m)];
    maxWavenumber = 0.0;
    for (int idx : mode.batterIndex)
        maxWavenumber = std::max (maxWavenumber, list[static_cast<std::size_t> (idx)].b1);
    return [&mode, &list, radius] (double r)
    {
        double u = 0.0;
        for (std::size_t i = 0; i < mode.batterIndex.size(); ++i)
            u += mode.batterCoef[i] * evaluateShape (list[static_cast<std::size_t> (mode.batterIndex[i])], radius, r).value;
        return u;
    };
}

double rearVolume (const CoupledMode& mode, const std::vector<std::vector<HeadMode>>& rear)
{
    double volume = 0.0;
    if (mode.m == 0 && ! rear.empty())
        for (std::size_t j = 0; j < mode.rearIndex.size(); ++j)
            volume += mode.rearCoef[j] * rear[0][static_cast<std::size_t> (mode.rearIndex[j])].volume;
    return volume;
}

void computeMicTransfers (const DrumContext& ctx, const MicSpec& mics, CoupledMode& mode)
{
    const double a = ctx.air.radius, rho = ctx.air.density, c = ctx.air.soundSpeed, len = ctx.air.length;
    const std::complex<double> s (-mode.alpha, mode.omega);
    double b1Max = 0.0;
    const auto profile = batterProfile (mode, ctx.batter, a, b1Max);
    const double volume = rearVolume (mode, ctx.rear);
    for (int mic = 0; mic < 2; ++mic)
    {
        const auto& p = mics.capsule[static_cast<std::size_t> (mic)];
        const double rhoM = std::hypot (p[0], p[1]);
        const std::complex<double> g = rayleighIntegral (mode.m, a, profile, b1Max, s, c, rhoM, p[2], 1);
        mode.micBatter[static_cast<std::size_t> (mic)] = rho * s * s * g;
        if (volume != 0.0)
        {
            const double path = rearDistance (len, p);
            mode.micRear[static_cast<std::size_t> (mic)] = rho * s * s * volume * std::exp (-s * (path / c)) / (4.0 * pi * path);
        }
    }
}

// BODY. For the modes of the body table's band the boundary-element solution of
// the whole drum (ModalBodyExterior.h) gives the capsule pressure per unit
// normal velocity of each head column; the table keeps it minus the transfer
// above (baffled Rayleigh batter, compact-monopole rear) at the same frequency,
// so that summing the mode's columns gives dF(w), the body's addition, and
// micBody = s dF(w_k) (velocity s A). Evaluated at the real w_k (not at the
// pole: a relative error of order alpha tau, below 3 % for these modes).
// Head components above the table (none carry much of a mode in this band)
// keep the plain transfer.
//
// It is a diagnostic and is not added to the capsule signal: the output's
// exactness while the stick drives the head (causality, and a pressure whose
// time integral over the stroke vanishes; ModalFeedthrough.cpp) needs the
// correction's impulse response, of which the table knows only the frequency
// response below its top. Added as a pole-residue term alone (heard through a
// delay line of its first arrival at the rim corner, exact in free decay) it
// left 18 % of the uncorrected pressure integral, dominated by the lowest
// modes. Its size is in the diagnostics: the baffled batter plus compact rear
// monopole overstate the capsule level of the (0,2) and (1,2) modes by about
// 3-4 dB and the (0,1) pair by 1-2 dB.
void computeBodyTransfers (const DrumContext& ctx, CoupledMode& mode)
{
    const BodyExterior* body = ctx.body.get();
    if (body == nullptr || mode.high || mode.omega > body->maxOmega())
        return;
    std::vector<std::complex<double>> correction;
    body->fieldCorrection (mode.m, mode.omega, correction);
    const int nbc = body->batterColumns (mode.m), nrc = body->rearColumns (mode.m);
    const int nh = nbc + nrc;
    if (nh == 0 || correction.size() != 2 * static_cast<std::size_t> (nh))
        return;
    const std::complex<double> s (-mode.alpha, mode.omega);
    for (std::size_t mic = 0; mic < 2; ++mic)
    {
        std::complex<double> sum = 0.0;
        for (std::size_t i = 0; i < mode.batterIndex.size(); ++i)
            if (mode.batterIndex[i] < nbc)
                sum += mode.batterCoef[i] * correction[mic * static_cast<std::size_t> (nh) + static_cast<std::size_t> (mode.batterIndex[i])];
        for (std::size_t j = 0; j < mode.rearIndex.size(); ++j)
            if (mode.rearIndex[j] < nrc)
                sum += mode.rearCoef[j] * correction[mic * static_cast<std::size_t> (nh) + static_cast<std::size_t> (nbc + mode.rearIndex[j])];
        mode.micBody[mic] = s * sum;
    }
}

} // namespace taikor::candidates::modal
