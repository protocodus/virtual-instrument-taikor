#pragma once

// In-vacuo modes of a clamped stiff membrane (a hide under tension with
// bending stiffness), and the closed-form integrals of their shapes that the
// coupled-mode drum needs.
//
// The hide obeys sigma w_tt = T lap w - D lap^2 w (Suzuki, Yamaguchi & Izumi,
// "Theoretical and experimental studies on the resonance frequencies of a
// stretched circular plate", Acoust. Sci. & Tech. 30(5):348, 2009, Eq. 1;
// Leissa, "Vibration of Plates", NASA SP-160, 1969, sec. 2.1 with in-plane
// tension). Separating w = phi(r) cos(m (theta - theta_s)) e^{j w t} gives
// (lap + b1^2)(lap - b2^2) phi = 0 with
//   b1^2 = (sqrt(T^2 + 4 D sigma w^2) - T) / (2 D),  b2^2 = b1^2 + T / D,
// so phi = J_m(b1 r) - c I_m(b2 r), and the clamped edge (phi = phi' = 0 at
// r = a) fixes c = J_m(b1 a) / I_m(b2 a) and the frequency equation
//   b1 J_m'(b1 a) I_m(b2 a) - b2 I_m'(b2 a) J_m(b1 a) = 0
// (Suzuki 2009 Eq. 12). As D -> 0, b2 -> infinity and the roots b1 a tend to
// the zeros of J_m (the ideal clamped membrane).
//
// Only the cos(m (theta - theta_s)) partner of each degenerate pair is kept:
// a uniform head struck at azimuth theta_s never excites the sin partner.

#include "SpecialFunctions.h"

#include <vector>

namespace taikor::candidates::modal
{

inline constexpr double pi = special::pi;

// Lambda_m = int_0^{2 pi} cos^2(m theta) d theta.
[[nodiscard]] inline double angularNorm (int m) noexcept { return m == 0 ? 2.0 * pi : pi; }

struct HeadMaterial
{
    double tension = 0.0;          // T, N/m
    double arealDensity = 0.0;     // sigma, kg/m^2
    double bending = 0.0;          // D, N m
    double radius = 0.0;           // a, m
};

struct HeadMode
{
    int m = 0, n = 0;              // angular order, radial index (1-based)
    double b1 = 0.0, b2 = 0.0;     // rad/m
    double omega = 0.0;            // in-vacuo angular frequency, rad/s
    double jEdge = 0.0;            // J_m(b1 a)
    double iRatio = 0.0;           // I_m'(b2 a) / I_m(b2 a)
    double iEdgeScaled = 0.0;      // e^{-b2 a} I_m(b2 a)
    double q = 0.0;                // b1 J_m'(b1 a) (= c b2 I_m'(b2 a) by the edge condition)
    double normSq = 0.0;           // int phi^2 dS over the head (angular factor included), m^2
    double gradDiag = 0.0;         // b1^2 normSq: the diagonal part of int grad phi_i . grad phi_j dS
    double volume = 0.0;           // int phi dS (m = 0 only), m^2
};

// phi(r) and phi'(r). The I_m part is evaluated scaled, so b2 a ~ 200 cannot
// overflow: c I_m(b2 r) = J_m(b1 a) [e^{-b2 r} I_m(b2 r) / e^{-b2 a} I_m(b2 a)] e^{b2 (r - a)}.
struct ShapeValue
{
    double value = 0.0, slope = 0.0;
};
[[nodiscard]] ShapeValue evaluateShape (const HeadMode& mode, double radius, double r);

// Solves every cos-partner mode of one head with in-vacuo frequency below
// omegaMax. Roots are bracketed by a scan in b1 a at step 0.05 (consecutive
// roots are about pi apart, never closer than ~2) and refined to 1e-14.
// Mode norms and volumes are in closed form (Lommel integrals; ModalMath.cpp),
// checked against quadrature in tests/modal heads_closed_forms.
[[nodiscard]] std::vector<std::vector<HeadMode>> solveHeadModes (const HeadMaterial&, double omegaMax);

// W(kappa) = int_0^a phi(r) J_m(kappa r) r dr, the order-m Hankel transform of
// the radial shape, in closed form by the Lommel integrals (J with J, and J
// with I through J_m(i x) = i^m I_m(x)):
//   W = a (b1^2 + b2^2) [kappa J_m'(kappa a) J_m(b1 a) - Q J_m(kappa a)]
//       / ((b1^2 - kappa^2)(b2^2 + kappa^2)),   Q = b1 J_m'(b1 a).
// jm = J_m(kappa a), jmPrime = J_m'(kappa a). The singularity at kappa = b1 is
// removable: its limit is a^2/2 [J_m'(b1 a)^2 + (1 - m^2/(b1 a)^2) J_m(b1 a)^2],
// and within 1e-6 b1 of it the transform is the quadratic through that limit
// and the closed form 2e-6 b1 either side.
[[nodiscard]] double hankelTransform (const HeadMode&, double radius, double kappa, double jm, double jmPrime);

// Exact mean of phi cos(m(theta - theta_s)) over a uniform disk of radius b
// centred at (rs, theta_s): the mean-value theorem for (lap + b1^2) u = 0 and
// (lap - b2^2) u = 0 gives
//   J_m(b1 rs) 2 J_1(b1 b)/(b1 b) - c I_m(b2 rs) 2 I_1(b2 b)/(b2 b).
[[nodiscard]] double diskAverage (const HeadMode&, double radius, double rs, double b);

// The Berger gradient form of one head at one angular order,
//   Gamma_ij = int grad(phi_i cos) . grad(phi_j cos) dS = d_i delta_ij - sum_r l_ri l_rj,
// exactly Gamma = diag(b1_i^2 normSq_i) - C K C with C = diag(c_i (b1_i^2 + b2_i^2)),
// c_i = J_m(b1_i a)/I_m(b2_i a), K_ij = P_ij / (b1_i^2 + b1_j^2 + T/D) and
// P_ij = int I_m(b2_i r) I_m(b2_j r) dS. (Derivation: lap phi = -b1^2 phi - c A I_m(b2 r),
// and Green's identity with the clamped edge gives int phi_i I_m(b2_j r) dS in terms
// of P.) K is positive semi-definite (Schur product of a Gram matrix and a Cauchy
// matrix), so a pivoted Cholesky K ~ L L^T with PSD residual keeps
// D - C L L^T C >= Gamma >= 0. Its eigenvalues fall about tenfold per rank.
struct GradientForm
{
    std::vector<double> diagonal;              // d_i
    std::vector<std::vector<double>> lowRank;  // l_r (each of length N)
    std::vector<double> pivots;                // relative pivot of each l_r; any prefix leaves a PSD residual
    double truncation = 0.0;                   // largest discarded pivot, relative to max d_i
};
[[nodiscard]] GradientForm gradientForm (const std::vector<HeadMode>& modes, const HeadMaterial&,
                                         double relativeTolerance);

} // namespace taikor::candidates::modal
