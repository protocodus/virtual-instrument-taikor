#pragma once

// The air inside the barrel: a rigid closed cylinder of radius a and length L
// (heads at z = 0 and z = -L), coupled to both heads by the acoustoelastic
// modal equations of Dowell, Gorman & Smith, "Acoustoelasticity: general
// theory, acoustic natural modes and forced response to sinusoidal
// excitation", J. Sound Vib. 52(4):519-542 (1977):
//   P_s'' + Omega_s^2 P_s = -(rho c^2 / Lam_s) oint psi_s w_n'' dS,
//   M_i q_i'' + K_i q_i = sum_s C_is P_s + ...,  C_is = int_head phi_i psi_s dS,
// with rigid-wall modes psi_s = J_m(gamma_q r) cos(m theta) cos(p pi z / L),
// gamma_q = j'_mq / a (zeros of J_m'), Omega_s = c sqrt(gamma_q^2 + (p pi/L)^2),
// Lam_s = int psi_s^2 dV = N_q L eps_p (eps_0 = 1, eps_p = 1/2),
// N_q = Lam_m a^2/2 (1 - m^2/j'^2) J_m(j')^2 (Lam_0 = 2 pi, else pi), and
// w_n the OUTWARD normal displacement of each head (batter +z, rear -z).
//
// Symmetric form used here. Eliminating P_s in the frequency domain gives
// P_s = (rho c^2/Lam_s) lambda/(Omega_s^2 - lambda) (C^T q)_s with lambda = w^2.
// Writing lambda/(Omega^2 - lambda) = Omega^2/(Omega^2 - lambda) - 1 and
// introducing a unit-mass air coordinate z_s for the first part gives the
// real symmetric pencil
//   K = [ K_h + sum_s (rho c^2/Lam_s) c_s c_s^T ,  -c_s Omega_s sqrt(rho c^2/Lam_s) ]
//       [ symmetric                            ,   diag(Omega_s^2)                   ]
//   M = [ M_h + M_res , 0 ; 0 , I ],
// exact for every retained mode. The uniform (m = 0, q = 0, p = 0) mode has
// Omega = 0: it is the static compression of the closed barrel and enters as
// the stiffness (rho c^2 / V) c_0 c_0^T alone. The truncated modes enter
// through their quasi-static limit lambda/Omega^2 as the residual mass
//   M_res = sum_{s truncated} rho c^2 c_s c_s^T / (Lam_s Omega_s^2),
// summed over p in closed form per cross-mode q:
//   sum_p 1/(eps_p (gamma^2 + (p pi/L)^2)) = L coth(gamma L)/gamma  (same head),
//   sum_p (-1)^p/(eps_p (gamma^2 + (p pi/L)^2)) = L csch(gamma L)/gamma  (across),
// minus the retained terms; for q = 0 the zeta sums pi^2/6 and -pi^2/12.

#include "ModalMath.h"

#include <vector>

namespace taikor::candidates::modal
{

struct CrossMode
{
    double gamma = 0.0;            // rad/m
    double normQ = 0.0;            // N_q, m^2
    double jmEdge = 0.0;           // J_m(gamma a)
};

// Zeros of J_m' below xMax (for m = 0 the plane-wave zero at 0 first), by a
// scan at step 0.25 (consecutive zeros are never closer than ~2.4) and
// safeguarded Newton with J_m'' from Bessel's equation.
[[nodiscard]] std::vector<double> besselPrimeZeros (int m, double xMax);

[[nodiscard]] std::vector<CrossMode> crossModes (int m, double radius, double gammaMax);

struct AirMode
{
    int q = 0, p = 0;              // cross-mode index into the catalogue, axial index
    double omega = 0.0;            // rad/s
    double lambda = 0.0;           // Lam_s = N_q L eps_p
};

struct AirSpec
{
    double radius = 0.0, length = 0.0, density = 0.0, soundSpeed = 0.0;
};

// Rigid-cylinder modes of order m with 0 < Omega <= omegaMax (the m = 0, q = 0,
// p = 0 compression mode is excluded; see above).
[[nodiscard]] std::vector<AirMode> airModes (int m, const std::vector<CrossMode>&, const AirSpec&, double omegaMax);

// Head-side interior terms of one angular block. Head unknowns are ordered
// batter modes then rear modes; batterHankel[i][q] = W_i(gamma_q) (ModalMath.h)
// for every cross-mode q of the catalogue (for m = 0, q = 0 this is
// int phi r dr, so Lam_0 W = the mode's volume).
struct InteriorBlock
{
    int heads = 0, air = 0;
    std::vector<double> staticStiffness;   // heads x heads, row-major
    std::vector<double> residualMass;      // heads x heads
    std::vector<double> coupling;          // heads x air: -c_s Omega_s sqrt(rho c^2/Lam_s)
};
[[nodiscard]] InteriorBlock interiorBlock (int m, const std::vector<CrossMode>& cross,
                                           const std::vector<AirMode>& air, const AirSpec& spec,
                                           const std::vector<std::vector<double>>& batterHankel,
                                           const std::vector<std::vector<double>>& rearHankel);

// KIRCHHOFF SIDE-WALL LOSS (the shared decision for every candidate: the
// viscothermal boundary layer on the barrel's rigid interior side wall; the end
// caps are the heads, whose own boundary layers are not included; no bulk
// atmospheric absorption). The time-averaged power lost per unit area of a
// rigid wall under an acoustic field of pressure p and tangential particle
// velocity v_t just outside the boundary layer is
//   1/2 sqrt(w nu/2) [rho |v_t|^2 + (gamma - 1) |p|^2 / (rho c^2 sqrt(Pr))]
// (G. Kirchhoff, "Ueber den Einfluss der Waermeleitung in einem Gase auf die
// Schallbewegung", Ann. Phys. Chem. 134:177-193, 1868; A. D. Pierce, Acoustics,
// 3rd ed., Springer/ASA 2019, ch. 10, sec. 10.5; Morse & Ingard, Theoretical
// Acoustics, 1968, sec. 6.4), nu = mu/rho. For a plane wave in a duct of radius
// a (|v_t| = |p|/(rho c)) this is the decision's law, an amplitude decay per
// unit time of
//   alpha(w) = (1/a) sqrt(w nu/2) (1 + (gamma - 1)/sqrt(Pr)),
// 12.4 / 14.3 / 16.0 / 17.5 dB/s at 3 / 4 / 5 / 6 kHz for the nagado barrel
// (tests/modal wall_loss_kirchhoff). The viscous term follows the TANGENTIAL
// velocity, so a barrel mode whose velocity at the side wall is normal (the
// radial modes m = 0, p = 0) loses only the thermal part there.
//
// One air mode s on its own, with v_t = grad_t p / (w rho):
//   alpha = (c^2/2) sqrt(w nu/2) oint_side (|grad_t psi|^2/w^2 + (gamma-1) psi^2/(c^2 sqrt Pr)) dS / Lam_s.
// A coupled mode carries the air at its own frequency w, not at Omega_s:
// with omega > 0 the result is half the damping coefficient on the unit-mass
// air coordinate at w (the boundary layer's sqrt(w), its tangential velocity
// grad p/(w rho), and the air mode's pressure amplitude (w^2/Omega) sqrt(rho c^2/Lam) z
// in the symmetric pencil above); at w = Omega_s it is the mode's own decay rate.
[[nodiscard]] double wallLossRate (int m, const CrossMode&, const AirMode&, const AirSpec&, double omega = 0.0);

// The decision's plane-wave law alone: (1/a) sqrt(w nu/2) (1 + (gamma - 1)/sqrt Pr).
[[nodiscard]] double kirchhoffPlaneWaveRate (const AirSpec&, double omega);

// The same law applied to the WHOLE wall field of a coupled mode of order m at
// its frequency omega (what the engine uses). The loss is quadratic in the
// total field at the wall, and the side-wall traces of different cross-modes of
// the same axial index p are not orthogonal (J_m(gamma_q a) J_m(gamma_q' a) != 0),
// so the per-mode rates above do not add. With the field at the wall written as
//   p(a, theta, z) = sqrt(rho c^2) cos(m theta) sum_p S_p cos(p pi z / L),
//   S_p = sum_{q} J_m(gamma_q a) (w^2/Omega_s) x_s / sqrt(Lam_s)   (air coordinates x_s, s = (q, p))
//         + [m = 0, p = 0] P_0 / sqrt(rho c^2),
// where P_0 = -(rho c^2 / V) (volume displaced outward by both heads) is the
// barrel's uniform compression (the Omega = 0 mode, which the pencil carries as
// stiffness), the axial functions are orthogonal on the wall, and
//   alpha = (1/w^2) sum_p (c^2/2) sqrt(w nu/2) (tau_p/w^2 + (gamma-1) pi_p/(c^2 sqrt Pr)) S_p^2,
//   pi_p  = oint cos^2(m theta) cos^2(p pi z/L) dS = a Lam_m L eps_p,
//   tau_p = a [ (m/a)^2 Lam'_m L eps_p + (p pi/L)^2 Lam_m L (1 - eps_p) ]   (Lam'_m = pi, 0 for m = 0).
// With one coordinate this is wallLossRate(omega) times x_s^2 (tests/modal
// wall_loss_field_quadrature checks the whole form against a direct quadrature
// of the local law over the wall). The truncated air modes above the explicit
// band would add their quasi-static wall pressure; they are not carried (a
// 6 kHz band instead of 4.5 kHz changes the 30 lowest modes' wall loss by at
// most 3.6 %, less than 1e-3 of their decay: tests/modal build_coupled_band).
// Returns alpha; with half given, appends one entry
// sqrt(2 w_p) S_p per axial index present (alpha = |half|^2 / 2), whose dot
// products between modes are the wall part of the damping couplings.
[[nodiscard]] double wallLoss (int m, const std::vector<CrossMode>& cross, const std::vector<AirMode>& air,
                               const AirSpec& spec, const double* airCoordinates, double compression, double omega,
                               std::vector<double>* half = nullptr);

// Interior loading of one head mode treated on its own above the explicit
// air band: the evanescent cross-modes (gamma_q > k) of the closed duct,
//   m_int = rho sum_q coth(kappa_q L) W(gamma_q)^2 Lam_m^2 / (kappa_q N_q),  kappa_q = sqrt(gamma_q^2 - k^2),
// the same-face term of the exact duct impedance (j w rho/k_z) cot(k_z L) C^2/N_q
// at imaginary k_z. Propagating cross-modes (gamma_q < k) are left out: their
// reactance changes sign across every barrel resonance and has zero mean, and a
// resolved treatment would need the air modes above the explicit band, which
// this engine does not carry (see ModalBuild.cpp).
[[nodiscard]] double evanescentInteriorMass (const HeadMode&, double radius, double k, double length,
                                             double density, const std::vector<CrossMode>&);

} // namespace taikor::candidates::modal
