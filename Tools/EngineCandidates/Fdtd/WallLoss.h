#pragma once

// The viscothermal boundary layer on the barrel's interior side wall (rigid
// wood), as a passive frequency-dependent boundary of the leapfrog air.
//
// Physics. At a rigid wall the no-slip and isothermal conditions confine the
// air's viscous shear and heat conduction to layers of thickness
// sqrt(2 nu / w) and sqrt(2 chi / w) (nu = mu / rho, chi = nu / Pr). Seen from
// the acoustic field outside them they act as two boundary conditions
// (G. Kirchhoff, "Ueber den Einfluss der Waermeleitung in einem Gase auf die
// Schallbewegung", Ann. Phys. 210(6), 1868, 177-193; A. D. Pierce, Acoustics,
// 3rd ed., Springer/ASA 2019, sec. 10.5; M. Berggren, A. Bernland & D.
// Noreland, "Acoustic boundary layers as boundary conditions", J. Comput.
// Phys. 371, 2018, 633-650):
//   thermal: a volume flux into the wall per unit area u = C_th sqrt(s) p,
//            C_th = (gamma - 1) sqrt(chi) / (rho c^2);
//   viscous: a shear stress on the tangential velocity tau = rho sqrt(nu) sqrt(s) v_t
//            (Stokes' second problem),
// whose time-averaged power per unit area is
//   1/2 sqrt(w nu / 2) [rho |v_t|^2 + (gamma - 1) |p|^2 / (rho c^2 sqrt(Pr))],
// the loss the modal engine integrates over each barrel mode
// (Modal/ModalCavity.h wallLossRate); for a plane wave in a duct of radius a
// it is Kirchhoff's amplitude decay (1/a) sqrt(w nu / 2)(1 + (gamma-1)/sqrt(Pr)),
// 12.4 / 14.3 / 16.0 / 17.5 dB/s at 3 / 4 / 5 / 6 kHz for this barrel. Bulk
// (classical + molecular-relaxation) absorption is not included: it needs a
// humidity the drum spec does not carry.
//
// sqrt(s) is realised by its diffusive representation, a positive sum of
// first-order sections sqrt(s) ~ sum_i mu_i s / (s + xi_i) (G. Montseny,
// "Diffusive representation of pseudo-differential time-operators", ESAIM
// Proc. 5, 1998, 159-175; T. Helie & D. Matignon, "Diffusive representations
// for the analysis and simulation of flared acoustic pipes with visco-thermal
// losses", Math. Models Methods Appl. Sci. 16(4), 2006, 503-536): each section
// is a series resistance-compliance branch, positive real for mu_i >= 0, so the
// boundary is passive for any weights the fit returns (WallLossFit.cpp).
//
// Discretisation (WallLoss.cpp), energy-consistent with the leapfrog air in
// the manner of S. Bilbao, B. Hamilton, J. Botts & L. Savioja, "Finite volume
// time domain room acoustics simulation under general impedance boundary
// conditions", IEEE/ACM TASLP 24(1), 2016, 161-173 (every branch stepped by
// d_t. and mu_t., its storage and dissipation exact):
//  - thermal, per wall cell of true wall area A: the cell equation gains
//    -(c^2 k^2 / Delta^3) q^n, q^n = A C_th sum_i mu_i d_t.theta_i^n, with
//    d_t.theta_i + xi_i mu_t.theta_i = p^n = rho d_t.Psi^n. That is implicit
//    in the cell's own Psi^{n+1} only (a scalar per cell);
//  - viscous, per friction element e of area A_e: a linear combination of
//    face fluxes V_e = sum_f w_ef G_f estimates the tangential flow at the
//    wall, and each face carries G_f = (Psi_nb - Psi) - sum_e w_ef D_e,
//    D_e = sum_i a_i phi_i, a_i = A_e sqrt(nu) mu_i / Delta^3, with
//    d_t.phi_i + xi_i mu_t.phi_i = V_e^n. The cells use G^n, known at step n,
//    so this part is explicit, and it is passive for any weights w_ef (the
//    force goes back along the same combination that senses the flow).
// The scheme's energy is the air's (with G on those faces) plus the branch
// storage (A C_th mu_i xi_i / 4)(theta^{n+1}^2 + theta^n^2) and
// (rho Delta a_i xi_i / 4)(phi^{n+1}^2 + phi^n^2); it falls by
// k A C_th mu_i (d_t.theta_i)^2 and k rho Delta a_i (d_t.phi_i)^2 per step.
//
// Staircase. The barrel is a staircase of x and y faces. A wall face whose
// centre sees the true wall's normal n carries the true area Delta^2 /
// (|n_x| + |n_y|) (the x faces of a monotone staircase project onto the
// true wall with |n_x|, the y faces with |n_y|); a uniform factor then makes
// them add up to the barrel's 2 pi a L exactly (WallLoss.cpp). Axial shear:
// one element per z face beside the wall, of half the wall area of the two
// cells it joins, sensing that face's flux (exact for axial flow: the z faces
// are tangential to every staircase face). Azimuthal shear: one element per
// wall cell, of the cell's wall area, sensing the tangential component
// -sin(theta) <G_x> + cos(theta) <G_y> of the cell centre's flux (the mean of
// its two x and of its two y faces, a wall face counting as zero; theta the
// centre's azimuth): on a staircase a face's own flux mixes the normal flow
// with the tangential one, the projection separates them (tests/fdtd
// wallDuctModes measures the resulting decay per duct-mode family).

#include "AirDomain.h"

#include <complex>
#include <utility>
#include <vector>

namespace taikor::candidates::fdtd
{

// Air at 20 degC (the temperature of the spec's rho = 1.2041 kg/m^3): dynamic
// viscosity 1.81e-5 Pa s, Prandtl number 0.71, ratio of specific heats 1.40
// (Kinsler, Frey, Coppens & Sanders, Fundamentals of Acoustics, 4th ed.,
// 2000, Appendix A10; the modal engine uses the same values).
struct Viscothermal
{
    double viscosity = 1.81e-5;
    double prandtl = 0.71;
    double gamma = 1.40;
};

// Weights of sqrt(s) ~ sum_i mu_i s / (s + xi_i) as the discrete branches
// realise it at step k (WallLossFit.cpp), one set per boundary kind.
struct SquareRootFit
{
    std::vector<double> xi;       // poles, 1/s
    std::vector<double> thermal;  // mu_i matching sqrt(j w)
    std::vector<double> viscous;  // mu_i matching cos^2(w k / 2) sqrt(j w)
    double bandLow = 0.0, bandHigh = 0.0;           // fit band, Hz
    // Worst |relative error| of the real part (the loss) from bandLow up to
    // errorHigh, the air grid's valid band (on-axis cut-off 8.87 kHz).
    double errorHigh = 0.0, thermalError = 0.0, viscousError = 0.0;
};

[[nodiscard]] SquareRootFit fitSquareRoot (double k);

// The discrete branch sum at angular frequency w (the admittance the scheme
// realises, relative to C_th or rho sqrt(nu)): sum_i mu_i B_i(w) with
// B_i = j sin(wk)/k / (j sin(wk)/k + xi_i cos(wk)).
[[nodiscard]] std::complex<double> realisedSquareRoot (const SquareRootFit& fit, const std::vector<double>& weights,
                                                       double k, double w);

template <typename Real>
class WallLoss
{
public:
    // Finds the wall cells (interior air cells with a face on the barrel's
    // solid wall) and the friction elements in `air` (after closeHeadFaces and
    // finalise). Slice cells (the air planes the head coupling steps) are
    // listed apart: their terms run inside the coupling's phases.
    void build (const AirDomain<Real>& air, double radius, double depth, const SquareRootFit& fit,
                const Viscothermal& gas = {});
    [[nodiscard]] bool active() const noexcept { return ! cells.empty(); }

    // Geometry, for diagnostics and tests.
    int wallCellCount() const noexcept { return wallCells; }
    int frictionFaceCount() const noexcept { return static_cast<int> (faces.size()); }
    int frictionCount() const noexcept { return static_cast<int> (elements.size()); } // friction elements
    double wallArea() const noexcept { return thermalArea; }   // sum of the wall faces' areas, m^2 (= 2 pi a L)
    double rawWallArea() const noexcept { return staircaseArea; } // the same before the normalisation
    double axialArea() const noexcept { return axialFriction; } // friction area on z faces, m^2
    double azimuthalArea() const noexcept { return azimuthalFriction; }
    int poles() const noexcept { return static_cast<int> (xi.size()); }

    // One step n -> n+1 (after the stencil kernels have written level n+1):
    // the friction elements [begin, end)'s branches to n+1 (reads level n only),
    double stepFriction (long long n, int begin, int end) noexcept;
    // and the ordinary (non-slice) cells [begin, end): viscous corrections,
    // the thermal implicit update and the thermal branches.
    // With tracking on (tests), both return the energy their range
    // dissipated in the step, J (0 otherwise).
    double finishCells (long long n, int begin, int end) noexcept;
    int ordinaryCellCount() const noexcept { return ordinaryCells; }

    // Slice cells: entry of an air cell (-1 when it carries no wall term),
    // then per step the part before the head solve (viscous corrections and
    // the thermal partial update, after AirDomain::stepMasked), the factor
    // 1 / (1 + beta) the head's own flux into the cell is scaled by, and the
    // thermal branches once the cell's level n+1 is final.
    [[nodiscard]] int sliceEntry (int cell) const;
    void slicePartial (long long n, int entry) noexcept;
    [[nodiscard]] double sliceFactor (int entry) const noexcept;
    double sliceFinish (long long n, int entry) noexcept;

    // Energy the wall adds to AirDomain::energy (n) at n + 1/2 (after step n):
    // the sheared faces' flux correction and the branches' storage, J.
    [[nodiscard]] double energy (long long n) const;
    void setTracking (bool on) noexcept { tracking = on; }

    // Tests: the strength of both parts (1 = Kirchhoff); 0 leaves the wall
    // rigid. It changes sliceFactor(), so an engine would have to rebuild its
    // coupling; only the standalone barrel of tests/fdtd calls it.
    void scale (double thermalScale, double viscousScale);

private:
    struct Cell
    {
        int cell = 0;            // air cell index
        int faceBegin = 0, faceEnd = 0; // its sheared faces in cellFaces
        int thermal = -1;        // index into the thermal arrays, -1: no wall face
    };
    struct Face
    {
        int low = 0, high = 0;   // cells on the - and + sides
    };
    struct Element
    {
        int begin = 0, end = 0;  // its faces in elementFaces / elementWeights
        double area = 0.0;       // friction area, m^2
        bool axial = false;
    };
    const AirDomain<Real>* domain = nullptr;
    std::vector<double> xi, muThermal, muViscous;
    std::vector<Cell> cells;     // ordinary cells first, then slice cells
    int ordinaryCells = 0, wallCells = 0;
    std::vector<int> cellFaces;  // 2 * face + (1 if the cell is the face's + side)
    std::vector<Face> faces;
    std::vector<Element> elements;
    std::vector<int> elementFaces;       // CSR of each element's faces
    std::vector<double> elementWeights;
    std::vector<int> faceStart, faceElements; // CSR of each face's elements (transpose)
    std::vector<double> faceWeights;
    std::vector<double> faceGain; // sqrt(nu) A_e / Delta^3 per element (a_i = faceGain mu_i)
    std::vector<double> cellArea; // true wall area per thermal cell, m^2
    std::vector<double> beta;     // per thermal cell
    std::vector<double> theta[2]; // [level & 1][thermal * poles + i]
    std::vector<double> phi[2];   // [level & 1][element * poles + i]
    std::vector<double> drag[2];  // [level & 1][element]: D_e = sum_i a_i phi_i
    std::vector<std::pair<int, int>> sliceLookup; // (cell, entry), sorted
    // Per pole: (1 - k xi) / (1 + k xi), 2k / (1 + k xi), 1 / (1 + k xi).
    std::vector<double> decay, gain, relief;
    double thermalArea = 0.0, staircaseArea = 0.0, axialFriction = 0.0, azimuthalFriction = 0.0;
    double lambda2 = 0.0, k = 0.0, rho = 0.0, c = 0.0, delta = 0.0, cTh = 0.0, rootNu = 0.0;
    double thermalScale = 1.0, viscousScale = 1.0;
    bool tracking = false;
    void setConstants();
    // sum_e w_ef D_e^n on face f, and the viscous correction of cell entry e at step n.
    double faceDrag (int f, long long n) const noexcept;
    double correction (const Cell& e, long long n) const noexcept;
    // Thermal: Psi^{n+1} from the explicit value, and the branches to n+1.
    double thermalUpdate (int t, long long n, double explicitValue, double previous) const noexcept;
    double thermalBranches (int t, long long n, double next, double previous) noexcept;
};

} // namespace taikor::candidates::fdtd
