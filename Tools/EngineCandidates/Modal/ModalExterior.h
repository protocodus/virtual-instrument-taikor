#pragma once

// Exterior air loading of a head that radiates into its own baffled
// half-space: the King (Hankel-transform) form of the baffled-plane radiation
// impedance between two modes of the same angular order m,
//   Z_ij(k) = rho c k Lam_m int_0^inf W_i(kappa) W_j(kappa) kappa dkappa / sqrt(k^2 - kappa^2),
// with W the order-m Hankel transform of the radial shape (ModalMath.h) and,
// for kappa > k, sqrt(k^2 - kappa^2) = -j sqrt(kappa^2 - k^2) (e^{jwt}).
// (L. V. King, "On the acoustic radiation field of the piezo-electric
// oscillator...", Can. J. Res. 11:135, 1934; E. G. Williams, Fourier
// Acoustics, Academic Press 1999, ch. 2 and 4.) Its real part is the radiation
// resistance R_ij (power to the far field); its imaginary part over omega is
// the added mass
//   m_ij = rho Lam_m int_k^inf W_i W_j kappa dkappa / sqrt(kappa^2 - k^2),
//   R_ij = rho c k Lam_m int_0^{pi/2} W_i W_j (k sin t) k sin t dt   (kappa = k sin t).
// Modes of different m, and the two heads (separate half-spaces), do not couple.
// That neglect is only acceptable where a mode's near field stays on its own
// head: this form serves the lone batter modes above the coupled band and,
// inside the coupled blocks, every pair that involves a head mode above the
// body table; below it the drum body's exterior replaces it (ModalBodyExterior.h:
// one field for both heads, 2.8 s against 0.87 s T60 for the push-pull (0,1)).
//
// QUADRATURE (numerical parameters, checked in tests/modal/TestExterior.cpp
// against the exact baffled piston and by doubling each):
//  - [0, k]: Gauss-Legendre in t with 24 + ceil(1.5 k a) nodes;
//  - [k, kappa_p]: kappa = k cosh u, Gauss-Legendre in u with 32 nodes, where
//    kappa_p is the first panel boundary at or beyond k + 2 pi/a;
//  - [kappa_p, kappaMax]: fixed panels of width pi/a (half a period of
//    J_m(kappa a)) with 8 Gauss nodes each, shared by every frequency so the
//    Bessel table is built once;
//  - kappaMax = 12 max(b1) + 60/a: beyond it the transforms fall as kappa^-2 or
//    faster and the tail is below 1e-8 of the self terms (TestExterior).

#include "ModalMath.h"

#include <functional>
#include <vector>

namespace taikor::candidates::modal
{

class ExteriorTable
{
public:
    ExteriorTable (double radius, double kappaMax, int maxOrder);

    // Added-mass and resistance matrices (row-major N x N) of modes that all
    // share one angular order, at wavenumber k.
    void addedMass (const std::vector<const HeadMode*>& modes, double k, double density, std::vector<double>& out) const;
    void resistance (const std::vector<const HeadMode*>& modes, double k, double density, double soundSpeed,
                     std::vector<double>& out) const;
    // The self terms of one mode (the diagonal entries of the above).
    [[nodiscard]] double selfAddedMass (const HeadMode&, double k, double density) const;
    [[nodiscard]] double selfResistance (const HeadMode&, double k, double density, double soundSpeed) const;

    [[nodiscard]] double radius() const noexcept { return a; }
    [[nodiscard]] double maxKappa() const noexcept { return kappaMax; }
    [[nodiscard]] int maxOrder() const noexcept { return orders; }

private:
    double a, kappaMax;
    int orders;
    std::vector<double> panelKappa, panelWeight;   // fixed panel nodes
    std::vector<double> panelJ, panelJp;           // [node * (orders + 1) + m]

    struct Node
    {
        double kappa, weight, j, jp;
    };
    // Nodes and weights of the added-mass integral int f kappa dkappa / sqrt(kappa^2 - k^2)
    // for order m (weights include the kernel), and of the resistance integral.
    void reactanceNodes (int m, double k, std::vector<Node>& nodes) const;
    void resistanceNodes (int m, double k, std::vector<Node>& nodes) const;
};

// Generic King integrals for a transform W(kappa) given as a function (used by
// the tests with the exact piston transform a J_1(kappa a)/kappa). Returns the
// added mass (kg) and resistance (kg/s) for order m.
void kingIntegrals (int m, double radius, double k, double kappaMax, double density, double soundSpeed,
                    const std::function<double (double)>& transform, double& addedMassOut, double& resistanceOut);

} // namespace taikor::candidates::modal
