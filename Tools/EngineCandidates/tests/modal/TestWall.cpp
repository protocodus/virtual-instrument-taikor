// The Kirchhoff side-wall loss of the shared decision (ModalCavity.h): the
// plane-wave law with the decision's constants, the per-mode closed form that
// follows from the local law, the whole-wall-field form against a direct
// quadrature of the local law, and the block's use of it with rigid heads.

#include "ModalBuild.h"
#include "ModalCavity.h"
#include "ModalFixtures.h"
#include "ModalTestUtil.h"

#include <algorithm>
#include <random>

using namespace taikor::candidates;
using namespace taikor::candidates::modal;

namespace
{
// The decision's air constants, restated so the engine is checked against them:
// dynamic viscosity 1.81e-5 Pa s (20 degC), Prandtl number 0.71, gamma 1.4.
constexpr double viscosity = 1.81e-5, prandtl = 0.71, gammaAir = 1.4;
constexpr double dbPerNeper = 8.685889638065035;   // 20 log10(e)

AirSpec nagadoAir()
{
    const DrumSpec d = modaltest::nagadoDrum();
    return { d.radius, d.depth, d.airDensity, d.soundSpeed };
}

// The decision's law: a plane wave in a duct of radius a decays at
// (1/a) sqrt(w nu/2) (1 + (gamma - 1)/sqrt Pr) per unit time.
double planeWaveLaw (const AirSpec& s, double w)
{
    const double nu = viscosity / s.density;
    return std::sqrt (w * nu / 2.0) * (1.0 + (gammaAir - 1.0) / std::sqrt (prandtl)) / s.radius;
}

// The local law integrated over the side wall for one rigid-cylinder mode
// (m, j'_mq, p) at its own frequency, by hand: with gamma = j'/a and k = p pi/L,
//   alpha = sqrt(w nu/2) [ (m^2/a^2 + k^2)/(gamma^2 + k^2) + (gamma_air - 1)/sqrt Pr ] / (a (1 - m^2/j'^2)).
// (1 - m^2/j'^2 is 1 for m = 0.) Only the velocity tangential to the wall enters
// the viscous part, so the radial modes (m = 0, p = 0) keep only the thermal one.
double modeLaw (const AirSpec& s, int m, double jPrime, int p, double w)
{
    const double a = s.radius, k = p * pi / s.length, g = jPrime / a;
    const double nu = viscosity / s.density;
    const double shape = m == 0 ? 1.0 : 1.0 - double (m * m) / (jPrime * jPrime);
    const double viscous = (double (m * m) / (a * a) + k * k) / (g * g + k * k);
    return std::sqrt (w * nu / 2.0) * (viscous + (gammaAir - 1.0) / std::sqrt (prandtl)) / (a * shape);
}
} // namespace

MODAL_TEST (wall_loss_kirchhoff)
{
    const AirSpec spec = nagadoAir();
    // Plane waves (m = 0, q = 0) to 8 kHz: the engine's rate is the decision's law.
    {
        const auto cross = crossModes (0, spec.radius, 1.0);
        const auto air = airModes (0, cross, spec, 2.0 * pi * 8000.0);
        double worst = 0.0;
        int count = 0;
        for (const auto& s : air)
            if (s.q == 0)
            {
                const double rate = wallLossRate (0, cross[0], s, spec);
                worst = std::max (worst, std::abs (rate / planeWaveLaw (spec, s.omega) - 1.0));
                ++count;
            }
        std::printf ("  %d plane-wave modes to 8 kHz\n", count);
        modaltest::below ("plane-wave modes: max relative deviation from the decision's law", worst, 1e-12);
        modaltest::nearRel ("the engine's plane-wave law at 3 kHz (its constants)", kirchhoffPlaneWaveRate (spec, 2.0 * pi * 3000.0),
                            planeWaveLaw (spec, 2.0 * pi * 3000.0), 1e-14);
        // The decision quotes 12.4 / 14.3 / 16.0 / 17.5 dB/s at 3 / 4 / 5 / 6 kHz for
        // this barrel; the plane-wave mode nearest each, moved there by sqrt(w).
        const double target[4][2] = { { 3000.0, 12.4 }, { 4000.0, 14.3 }, { 5000.0, 16.0 }, { 6000.0, 17.5 } };
        for (const auto& t : target)
        {
            const AirMode* best = nullptr;
            for (const auto& s : air)
                if (s.q == 0 && (best == nullptr || std::abs (s.omega - 2.0 * pi * t[0]) < std::abs (best->omega - 2.0 * pi * t[0])))
                    best = &s;
            const double db = dbPerNeper * wallLossRate (0, cross[0], *best, spec) * std::sqrt (2.0 * pi * t[0] / best->omega);
            char what[96];
            std::snprintf (what, sizeof (what), "plane-wave decay at %.0f Hz (dB/s)", t[0]);
            modaltest::near (what, db, t[1], 0.05);
        }
    }
    // Every rigid-cylinder mode to 8 kHz for several orders: the side-wall
    // integrals of the engine against the closed form above.
    double worst = 0.0, lowest = 1e9, highest = 0.0;
    int count = 0;
    for (int m : { 0, 1, 2, 5, 12, 30 })
    {
        const auto cross = crossModes (m, spec.radius, 2.0 * pi * 8000.0 / spec.soundSpeed + 1.0);
        for (const auto& s : airModes (m, cross, spec, 2.0 * pi * 8000.0))
        {
            const double jPrime = cross[static_cast<std::size_t> (s.q)].gamma * spec.radius;
            const double law = modeLaw (spec, m, jPrime, s.p, s.omega);
            worst = std::max (worst, std::abs (wallLossRate (m, cross[static_cast<std::size_t> (s.q)], s, spec) / law - 1.0));
            lowest = std::min (lowest, law / planeWaveLaw (spec, s.omega));
            highest = std::max (highest, law / planeWaveLaw (spec, s.omega));
            ++count;
        }
    }
    std::printf ("  %d modes of orders 0, 1, 2, 5, 12, 30 to 8 kHz: rate / plane-wave law at the same frequency %.3f .. %.3f\n",
                 count, lowest, highest);
    modaltest::below ("all modes: max relative deviation from the closed form", worst, 1e-12);
    // The radial modes' thermal-only share: (gamma - 1)/sqrt Pr / (1 + (gamma - 1)/sqrt Pr).
    modaltest::nearRel ("smallest ratio = the radial modes' thermal share", lowest,
                        (gammaAir - 1.0) / std::sqrt (prandtl) / (1.0 + (gammaAir - 1.0) / std::sqrt (prandtl)), 1e-12);
}

MODAL_TEST (wall_loss_field_quadrature)
{
    // A coupled mode's wall field from arbitrary air coordinates (and, for m = 0,
    // a uniform compression) at a frequency w off every air mode: the engine's
    // wallLoss against the local law integrated directly over the side wall,
    //   alpha = (1/w^2) oint 1/2 sqrt(w nu/2) [ |grad_t p|^2/(w^2 rho) + (gamma-1) |p|^2/(rho c^2 sqrt Pr) ] dS,
    // p = sqrt(rho c^2) sum_s J_m(gamma_q a) (w^2/Omega_s) x_s/sqrt(Lam_s) cos(m theta) cos(p_s pi z/L) + P_0,
    // by trapezoid in theta (exact for the trigonometric integrand) and Gauss-Legendre in z.
    const AirSpec spec = nagadoAir();
    const double rhoC2 = spec.density * spec.soundSpeed * spec.soundSpeed;
    const double nu = viscosity / spec.density;
    const double w = 2.0 * pi * 1234.5;
    std::mt19937 rng (7);
    std::normal_distribution<double> normal;
    for (int m : { 0, 1, 3 })
    {
        const auto cross = crossModes (m, spec.radius, 2.0 * pi * 3000.0 / spec.soundSpeed + 1.0);
        const auto air = airModes (m, cross, spec, 2.0 * pi * 3000.0);
        std::vector<double> x (air.size());
        for (auto& v : x)
            v = normal (rng) * 1e-3;
        const double compression = m == 0 ? normal (rng) * 1e-3 : 0.0;
        std::vector<double> half;
        const double alpha = wallLoss (m, cross, air, spec, x.data(), compression, w, &half);
        int maxP = 0;
        for (const auto& s : air)
            maxP = std::max (maxP, s.p);
        std::vector<double> zn, zw;
        special::gaussLegendre (2 * maxP + 40, zn, zw);
        const int nTheta = 4 * m + 16;
        double power = 0.0;
        for (std::size_t iz = 0; iz < zn.size(); ++iz)
        {
            const double z = -0.5 * spec.length * (1.0 + zn[iz]);
            const double dz = 0.5 * spec.length * zw[iz];
            double f = m == 0 ? compression : 0.0, fz = 0.0;   // radial profile at the wall: value and d/dz
            for (std::size_t s = 0; s < air.size(); ++s)
            {
                const double k = air[s].p * pi / spec.length;
                const double u = cross[static_cast<std::size_t> (air[s].q)].jmEdge * (w * w / air[s].omega) * x[s]
                               / std::sqrt (air[s].lambda);
                f += u * std::cos (k * z);
                fz -= u * k * std::sin (k * z);
            }
            for (int it = 0; it < nTheta; ++it)
            {
                const double theta = 2.0 * pi * it / nTheta;
                const double p = std::sqrt (rhoC2) * f * std::cos (m * theta);
                const double pTheta = -std::sqrt (rhoC2) * f * m * std::sin (m * theta) / spec.radius;
                const double pZ = std::sqrt (rhoC2) * fz * std::cos (m * theta);
                const double density = 0.5 * std::sqrt (w * nu / 2.0)
                                      * ((pTheta * pTheta + pZ * pZ) / (w * w * spec.density)
                                         + (gammaAir - 1.0) * p * p / (rhoC2 * std::sqrt (prandtl)));
                power += density * spec.radius * (2.0 * pi / nTheta) * dz;
            }
        }
        const double direct = power / (w * w);
        double diagonal = 0.0, halfNorm = 0.0;
        for (std::size_t s = 0; s < air.size(); ++s)
            diagonal += wallLossRate (m, cross[static_cast<std::size_t> (air[s].q)], air[s], spec, w) * x[s] * x[s];
        for (double h : half)
            halfNorm += h * h;
        std::printf ("  m %d: %zu air modes (p to %d); field %.6g 1/s, direct %.6g, sum of per-mode rates %.6g\n", m,
                     air.size(), maxP, alpha, direct, diagonal);
        modaltest::nearRel ("whole-field wall loss vs direct quadrature of the local law", alpha, direct, 1e-10);
        modaltest::nearRel ("|half|^2 / 2 = alpha", 0.5 * halfNorm, alpha, 1e-13);
        // One coordinate alone: the per-mode rate at w.
        for (std::size_t s : { std::size_t (0), air.size() / 2, air.size() - 1 })
        {
            std::vector<double> e (air.size(), 0.0);
            e[s] = 1.0;
            const double one = wallLoss (m, cross, air, spec, e.data(), 0.0, w);
            modaltest::nearRel ("single coordinate = wallLossRate(w)", one,
                                wallLossRate (m, cross[static_cast<std::size_t> (air[s].q)], air[s], spec, w), 1e-13);
        }
    }
}

MODAL_TEST (wall_loss_block_rigid_heads)
{
    // Heads a million times heavier (cavity_rigid_heads): each air-dominated
    // coupled mode is one rigid-cylinder mode at its own frequency, and its decay
    // is that mode's Kirchhoff rate (closed form above; the heads' hide loss is
    // negligible at their tiny share and the exterior is off).
    DrumSpec d = modaltest::nagadoDrum();
    d.batter.arealDensity *= 1e6;
    d.rear.arealDensity *= 1e6;
    BuildOptions o;
    o.exterior = false;
    o.headMaxHz = o.coupledMaxHz = 3.0;
    o.airMaxHz = 1500.0;
    const AirSpec spec = nagadoAir();
    for (int m : { 0, 1, 3 })
    {
        const auto block = solveBlockForTest (d, o, m);
        const auto cross = crossModes (m, d.radius, 2.0 * pi * o.airMaxHz / d.soundSpeed + 40.0 / d.radius);
        auto air = airModes (m, cross, spec, 2.0 * pi * o.airMaxHz);
        std::sort (air.begin(), air.end(), [] (const auto& x, const auto& y) { return x.omega < y.omega; });
        std::vector<std::size_t> coupled;
        for (std::size_t k = 0; k < block.omega.size(); ++k)
            if (block.shareAir[k] > 0.5)
                coupled.push_back (k);
        double worst = 0.0;
        int compared = 0;
        for (std::size_t i = 0; i < std::min (coupled.size(), air.size()); ++i)
        {
            // Skip near-degenerate pairs, which the heads' tiny coupling may mix.
            const bool isolated = (i == 0 || air[i].omega - air[i - 1].omega > 1e-3 * air[i].omega)
                               && (i + 1 == air.size() || air[i + 1].omega - air[i].omega > 1e-3 * air[i].omega);
            if (! isolated)
                continue;
            const double jPrime = cross[static_cast<std::size_t> (air[i].q)].gamma * d.radius;
            const double law = modeLaw (spec, m, jPrime, air[i].p, air[i].omega);
            worst = std::max (worst, std::abs (block.alpha[coupled[i]] / law - 1.0));
            ++compared;
        }
        std::printf ("  m %d: %d isolated air modes compared\n", m, compared);
        modaltest::below ("air-dominated decay vs the mode's Kirchhoff rate (relative)", worst, 1e-4);
    }
}
