// Test 7: two identical heads coupled through the barrel's air, the
// Suzuki & Hwang fixture (H. Suzuki & Y.-F. Hwang, "Coupling between two
// membranes of a Japanese drum", Acoust. Sci. Tech. 29(3), 2008, 215-223,
// sec. 4: R 0.20 m, L 0.50 m, 4 kg/m^2, T 14 kN/m, air 1.21 kg/m^3 at
// 343 m/s, no exterior air, rigid body, ideal membranes). Their bulged barrel
// gives the push-pull (in-phase) and push-push (out-of-phase) pair at
// 110.5 / 120.0 Hz (their FEM check, with a 1 MPa skin: 110.6 / 119.5 Hz);
// the shipping engine's CoupledCavity solver gives 110.1 / 122.6 Hz for the
// straight cylinder modelled here. The paper neglects bending stiffness, so
// the heads use its FEM check's 1 MPa skin (D ~ 7e-4 N m), on a 4 mm grid
// (above the stability bound of ~1 mm; the (0,1) then has beta h ~ 0.05).
//
// The gate is an independent semi-analytic solution of the same straight
// drum (ductPairReference below), not the shipping solver: the heads'
// axisymmetric clamped stiff modes coupled through the exact rigid-duct field
// of each radial order.

#include "TestSupport.h"

#include "Analytic.h"
#include "FdtdEngine.h"
#include "HeadCalibration.h"
#include "SpecialFunctions.h"

#include <cmath>

namespace
{
using namespace taikor::candidates;
using namespace taikor::candidates::fdtd;
using namespace taikor::candidates::fdtd::tests;

// Two clamped stiff heads (batter at z = 0, rear at z = -L, w along +z) on a
// rigid closed cylinder, axisymmetric motion. Heads: w = sum_n a_n Phi_n(r)
// with the in-vacuo modes Phi_n (Analytic.h), modal mass M_n = sigma int
// Phi_n^2 dA. Air: p = sum_q Psi_q(r) P_q(z), Psi_q = J0(k_q r), k_q =
// j'_0q / a (rigid wall), N_q = int Psi_q^2 dA, and P_q'' + kappa_q^2 P_q = 0
// with kappa_q^2 = (w/c)^2 - k_q^2 and the heads' normal acceleration as the
// end conditions dp/dz = rho w^2 w. Solving each duct order exactly gives the
// pressures on the two heads and, with C_nq = int Phi_n Psi_q dA, the
// symmetric system (e.g. P. M. Morse & K. U. Ingard, Theoretical Acoustics,
// 1968, sec. 9.2 for the duct field)
//   M_n (w_n^2 - w^2) a_n + sum_q rho w^2 / N_q C_nq (cot_q B_q - csc_q R_q) = 0
//   M_n (w_n^2 - w^2) b_n + sum_q rho w^2 / N_q C_nq (cot_q R_q - csc_q B_q) = 0
// B_q = sum_m C_mq a_m, R_q = sum_m C_mq b_m, cot_q = cos(kappa L) /
// (kappa sin(kappa L)), csc_q = 1 / (kappa sin(kappa L)) (their hyperbolic
// forms below cut-off). Its determinant vanishes at the coupled frequencies.
// Numerical truncation: head modes below 4 kHz, duct orders up to k_q a < 120,
// radial quadrature of 10-point Gauss panels graded into the rim's bending
// boundary layer (~1 / gamma = 0.2 mm); halving each changes the pair by
// < 1e-4 Hz.
std::vector<double> ductPairReference (const HeadSpec& head, double a, double L, double rho, double c,
                                       double lowHz, double highHz)
{
    std::vector<double> r, weight, gx, gw;
    special::gaussLegendre (10, gx, gw);
    std::vector<double> edges;
    for (int i = 0; i <= 300; ++i)
        edges.push_back (0.95 * a * i / 300.0);
    for (int i = 1; i <= 200; ++i)
        edges.push_back (a - 0.05 * a * std::pow (1.0e-6, i / 200.0));
    edges.push_back (a);
    for (std::size_t e = 0; e + 1 < edges.size(); ++e)
        for (std::size_t g = 0; g < gx.size(); ++g)
        {
            const double lo = edges[e], hi = edges[e + 1];
            r.push_back (lo + 0.5 * (hi - lo) * (gx[g] + 1.0));
            weight.push_back (0.5 * (hi - lo) * gw[g] * 2.0 * special::pi * r.back());
        }
    const auto modes = clampedStiffModes (head, a, 0, 4000.0);
    const auto zeros = special::besselJPrimeZeros (0, 120.0)[0];
    const std::size_t nm = modes.size(), nq = zeros.size();
    std::vector<double> omega (nm), mass (nm, 0.0), norm (nq, 0.0), C (nm * nq, 0.0);
    std::vector<double> phi (r.size()), psi (r.size());
    for (std::size_t q = 0; q < nq; ++q)
        for (std::size_t i = 0; i < r.size(); ++i)
            norm[q] += std::pow (special::besselJ (0, zeros[q] / a * r[i]), 2) * weight[i];
    for (std::size_t n = 0; n < nm; ++n)
    {
        omega[n] = 2.0 * special::pi * modes[n].frequency;
        for (std::size_t i = 0; i < r.size(); ++i)
        {
            phi[i] = clampedStiffShape (modes[n], r[i], a);
            mass[n] += head.arealDensity * phi[i] * phi[i] * weight[i];
        }
        for (std::size_t q = 0; q < nq; ++q)
            for (std::size_t i = 0; i < r.size(); ++i)
                C[n * nq + q] += phi[i] * special::besselJ (0, zeros[q] / a * r[i]) * weight[i];
    }
    // Sign of det(system) at w, by LU with partial pivoting, rows and columns
    // scaled by 1 / sqrt(M_n w_n^2) so the pivots stay O(1).
    const auto sign = [&] (double hz)
    {
        const double w = 2.0 * special::pi * hz, kk = w / c;
        const std::size_t size = 2 * nm;
        std::vector<double> A (size * size, 0.0);
        for (std::size_t q = 0; q < nq; ++q)
        {
            const double kz2 = kk * kk - std::pow (zeros[q] / a, 2);
            const double kz = std::sqrt (std::abs (kz2));
            const double cot = kz2 > 0.0 ? std::cos (kz * L) / (kz * std::sin (kz * L)) : -1.0 / (kz * std::tanh (kz * L));
            const double csc = kz2 > 0.0 ? 1.0 / (kz * std::sin (kz * L)) : -1.0 / (kz * std::sinh (kz * L));
            const double f = rho * w * w / norm[q];
            for (std::size_t n = 0; n < nm; ++n)
                for (std::size_t m = 0; m < nm; ++m)
                {
                    const double cc = f * C[n * nq + q] * C[m * nq + q];
                    A[n * size + m] += cc * cot;
                    A[(nm + n) * size + nm + m] += cc * cot;
                    A[n * size + nm + m] -= cc * csc;
                    A[(nm + n) * size + m] -= cc * csc;
                }
        }
        for (std::size_t n = 0; n < nm; ++n)
            for (const std::size_t d : { n, nm + n })
                A[d * size + d] += mass[n] * (omega[n] * omega[n] - w * w);
        for (std::size_t i = 0; i < size; ++i)
            for (std::size_t j = 0; j < size; ++j)
                A[i * size + j] /= std::sqrt (mass[i % nm] * omega[i % nm] * omega[i % nm]
                                              * mass[j % nm] * omega[j % nm] * omega[j % nm]);
        double s = 1.0;
        for (std::size_t col = 0; col < size; ++col)
        {
            std::size_t pivot = col;
            for (std::size_t i = col + 1; i < size; ++i)
                if (std::abs (A[i * size + col]) > std::abs (A[pivot * size + col]))
                    pivot = i;
            if (pivot != col)
            {
                for (std::size_t j = 0; j < size; ++j)
                    std::swap (A[col * size + j], A[pivot * size + j]);
                s = -s;
            }
            const double p = A[col * size + col];
            s *= p < 0.0 ? -1.0 : 1.0;
            for (std::size_t i = col + 1; i < size; ++i)
            {
                const double factor = A[i * size + col] / p;
                for (std::size_t j = col; j < size; ++j)
                    A[i * size + j] -= factor * A[col * size + j];
            }
        }
        return s;
    };
    std::vector<double> roots;
    const double step = 0.05;
    double previous = sign (lowHz);
    for (double f = lowHz + step; f <= highHz; f += step)
    {
        const double current = sign (f);
        if (current != previous)
        {
            double lo = f - step, hi = f;
            for (int it = 0; it < 60; ++it)
            {
                const double mid = 0.5 * (lo + hi);
                (sign (mid) == previous ? lo : hi) = mid;
            }
            roots.push_back (0.5 * (lo + hi));
        }
        previous = current;
    }
    return roots;
}
} // namespace

FDTD_TEST (suzukiHwangPair)
{
    DrumSpec drum;
    drum.radius = 0.20;
    drum.depth = 0.50;
    drum.shellWallThickness = 0.02; // rigid; its thickness does not enter without exterior air
    drum.airDensity = 1.21;
    drum.soundSpeed = 343.0;
    HeadSpec head;
    head.tension = 14000.0;
    head.arealDensity = 4.0;
    head.thickness = 0.002;
    head.youngsModulus = 1.0e6;
    head.poissonRatio = 0.3;
    head.bendingStiffness = head.youngsModulus * std::pow (head.thickness, 3) / (12.0 * (1.0 - 0.09));
    drum.batter = drum.rear = head;
    MicSpec mics;
    FdtdOptions options;
    options.lossless = true;
    options.exteriorAir = false;
    options.rigidBox = true;
    options.withBatterNonlinearity = false;
    options.bergerOnRear = false;
    options.selfChecks = false;
    options.threads = 2;
    options.minimumHeadSpacing = 0.004;
    FdtdEngine<double> engine (makeBergerPotential(), options);
    std::string error;
    if (! engine.prepare (drum, mics, 48000.0, error))
    {
        report.check ("prepare", false, error);
        return;
    }
    const auto fundamental = clampedStiffModes (head, drum.radius, 0, 300.0).front();
    const auto shape = modeShapeOnGrid (engine.batterState().grid, fundamental, false);
    engine.setHeadAtRest (0, shape);
    const int steps = static_cast<int> (0.45 * engine.internalRate());
    std::vector<double> sum, difference;
    const auto& grid = engine.batterState().grid;
    for (int n = 0; n < steps; ++n)
    {
        engine.advance (1);
        const long long s = engine.stepsTaken();
        const double b = grid.dot (engine.batterState().w[s % 3].data(), shape.data());
        const double r = grid.dot (engine.rearState().w[s % 3].data(), shape.data());
        sum.push_back (b + r);
        difference.push_back (b - r);
    }
    const double k = 1.0 / engine.internalRate();
    const double inPhase = spectralPeakFrequency (sum, k, 85.0, 118.0);
    const double outOfPhase = spectralPeakFrequency (difference, k, 112.0, 150.0);
    const auto reference = ductPairReference (head, drum.radius, drum.depth, drum.airDensity, drum.soundSpeed, 80.0, 150.0);
    if (reference.size() != 2)
    {
        report.check ("duct-mode reference", false, format ("expected 2 coupled frequencies in 80-150 Hz, found %zu", reference.size()));
        return;
    }
    report.note (format ("in vacuo (0,1) %.2f Hz (analytic clamped stiff); FDTD pair %.2f / %.2f Hz; independent "
                         "duct-mode solution %.3f / %.3f Hz; shipping CoupledCavity (cylinder) 110.1 / 122.6 Hz; paper "
                         "(bulged barrel) 110.5 / 120.0 Hz",
                         fundamental.frequency, inPhase, outOfPhase, reference[0], reference[1]));
    // Limit 1 %: the staircase heads (after their (0,1) correction) and the
    // staircase barrel (cavity modes within 1 %, test 2).
    report.check ("push-pull (in-phase) pair member", std::abs (inPhase / reference[0] - 1.0) < 0.01,
                  format ("%.2f Hz vs %.3f Hz (%+.2f %%), limit 1 %%", inPhase, reference[0], 100.0 * (inPhase / reference[0] - 1.0)));
    report.check ("push-push (out-of-phase) pair member", std::abs (outOfPhase / reference[1] - 1.0) < 0.01,
                  format ("%.2f Hz vs %.3f Hz (%+.2f %%), limit 1 %%", outOfPhase, reference[1],
                          100.0 * (outOfPhase / reference[1] - 1.0)));
}
