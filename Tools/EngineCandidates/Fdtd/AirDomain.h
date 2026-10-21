#pragma once

// The air inside and around the drum: a velocity potential Psi on a
// cell-centred cubic grid, stepped by the standard leapfrog 7-point scheme
//   Psi^{n+1} = 2 Psi^n - Psi^{n-1} + lambda^2 sum_{open faces} (Psi_nb - Psi),
// lambda = c k / Delta <= 1/sqrt(3), p = rho d_t Psi, v = -grad Psi
// (Bilbao, "Time domain simulation and sound synthesis for the snare drum",
// JASA 131(1), 2012, eqs. 17-20; Bilbao & Webb, JAES 61(10), 2013, eq. 20).
// A face is open when the cells on both sides are air and it is not a head.
// Dropping a face is a rigid (Neumann) wall that keeps the scheme's energy
//   H = rho Delta^3 / (2 c^2) sum_cells (d_t- Psi)^2
//       + rho Delta / 2 sum_{open faces} (dPsi^{n+1})(dPsi^n)
// exactly conserved, so the barrel is a staircase of dropped faces.
//
// Geometry (head-centre frame): batter at z = 0, rear at z = -L, barrel wall
// a <= r < a + wall for -L < z < 0, all rigid. Both heads lie on cell faces
// (L / Delta integer) with an air plane half a cell on each side. Faces sit
// at x = 0 and y = 0 so the grid has the head grid's square symmetry.
//
// Outer boundary: either a rigid box (energy tests) or a perfectly matched
// layer for the second-order wave equation (M. J. Grote & I. Sim, "Efficient
// PML for the wave equation", arXiv:1001.0319, 2010):
//   u_tt + (z1+z2+z3) u_t + (z1z2+z2z3+z3z1) u + z1z2z3 psi = c^2 Lap u + div phi
//   psi_t = u
//   phi_j,t = -z_j phi_j + c^2 (z_i + z_k - z_j) d_j u + c^2 z_i z_k d_j psi
// (derived here from the complex stretch s_j = 1 + z_j / (i w): multiply
// s^2 u = c^2 sum_j s_j^-1 d_j (s_j^-1 d_j u) by s1 s2 s3). u and psi live at
// cell centres, phi_j on j-faces at half steps; the profile is
// z(d) = zmax (d / thickness)^3.

#include <array>
#include <cstdint>
#include <vector>

namespace taikor::candidates::fdtd
{

struct AirGeometry
{
    double radius = 0.0, depth = 0.0, wall = 0.0; // drum: a, L, barrel wall
    double delta = 0.0;                           // cell size (depth / integer)
    double margin = 0.0;                          // air beyond the drum, m
    int pmlCells = 0;                             // 0: rigid box
    double pmlReflection = 1.0e-4;
    // Extra air above the batter (the microphone side), m, added to margin.
    double marginAbove = 0.0;
};

enum class CellType : std::uint8_t { Solid, Fixed, Regular, Irregular, Slice, Pml };

// Face bits of a cell: -x, +x, -y, +y, -z, +z.
enum FaceBit : std::uint8_t { FaceXm = 1, FaceXp = 2, FaceYm = 4, FaceYp = 8, FaceZm = 16, FaceZp = 32 };

template <typename Real>
class AirDomain
{
public:
    void build (const AirGeometry& geometry, double soundSpeed, double density, double k);

    int nx = 0, ny = 0, nz = 0;
    long long cells() const noexcept { return static_cast<long long> (nx) * ny * nz; }
    int index (int i, int j, int kz) const noexcept { return (kz * ny + j) * nx + i; }
    double x (int i) const noexcept { return x0 + (i + 0.5) * delta; }
    double y (int j) const noexcept { return y0 + (j + 0.5) * delta; }
    double z (int kz) const noexcept { return z0 + (kz + 0.5) * delta; }
    double spacing() const noexcept { return delta; }
    double courant() const noexcept { return lambda; }
    CellType type (int cell) const noexcept { return types[static_cast<std::size_t> (cell)]; }
    std::uint8_t faces (int cell) const noexcept { return masks[static_cast<std::size_t> (cell)]; }

    // Psi at time level n (levels rotate through three buffers).
    Real* level (long long n) noexcept { return psiBase[static_cast<std::size_t> (((n % 3) + 3) % 3)]; }
    const Real* level (long long n) const noexcept { return psiBase[static_cast<std::size_t> (((n % 3) + 3) % 3)]; }

    // One step n -> n+1, in pieces that may run concurrently (they write
    // disjoint cells of level n+1 and read levels n, n-1 only):
    // the regular cells of z-planes [k0, k1),
    void stepRegular (long long n, int k0, int k1) noexcept;
    // irregular cells [begin, end) of irregularCells(),
    void stepIrregular (long long n, int begin, int end) noexcept;
    // one cell with its own face mask (the head slices use this),
    void stepMasked (long long n, int cell) noexcept;
    // and the PML cells of z-planes [k0, k1): u to n+1, psi to n+1 and
    // their own faces' phi to n+1/2, in one pass (see AirDomainPml.cpp).
    void stepPml (long long n, int k0, int k1) noexcept;

    int irregularCount() const noexcept { return static_cast<int> (irregular.size()); }
    // Plane boundaries splitting the PML cells into `parts` ranges of nearly
    // equal cell count (the layer is full planes at the top and bottom and
    // thin strips in between, so equal plane counts would not balance).
    std::vector<int> pmlPlaneSplit (int parts) const;
    long long pmlCellCount() const;
    bool hasPml() const noexcept { return pml > 0; }

    // Acoustic pressure rho (Psi^{n+1} - Psi^{n-1}) / 2k at a point, trilinear.
    // Reads levels n+1 and n-1: call it before step n+1 starts.
    double pressure (long long n, const std::array<double, 3>& point) const noexcept;
    // Psi at level n at a point, trilinear.
    double potential (long long n, const std::array<double, 3>& point) const noexcept;

    // Discrete energy at n+1/2 (rigid box only; with a PML it is the energy
    // of the cells outside the layer), J.
    double energy (long long n) const;

    // Declares a head face: the z-face between planes kz - 1 and kz, for the
    // columns listed (i, j). With markSlices both adjacent cells become Slice
    // cells (stepped by the head coupling); otherwise the face is simply rigid.
    void closeHeadFaces (int kzAbove, const std::vector<std::array<int, 2>>& columns, bool markSlices = true);
    // The columns (i, j) whose centres lie inside the head, r < radius.
    std::vector<std::array<int, 2>> headColumns (double radius) const;
    // Rebuilds the regular runs and the irregular list after closeHeadFaces.
    void finalise();

    // The sealed barrel's air: every air cell with r_c < radius between the
    // head planes (listed once, after closeHeadFaces).
    void listInterior (double radius, double depth);
    int interiorCount() const noexcept { return static_cast<int> (interiorCells.size()); }
    // Adds `shift` to Psi at levels n and n-1 of interior cells [begin, end).
    // The interior touches the outside only through the heads, whose coupling
    // reads time differences of Psi, and its walls are rigid or carry the
    // boundary layer (WallLoss.h), whose thermal part reads time differences
    // and whose shear reads differences between interior cells, so a constant
    // shift of both levels changes nothing in exact arithmetic.
    void shiftInterior (long long n, double shift, int begin, int end) noexcept;
    double interiorSample (long long n) const noexcept;

    // z-plane index whose lower face is at height z (z must be a face plane).
    int planeAbove (double zFace) const noexcept;

    double density = 1.2041, c = 343.0, k = 0.0, delta = 0.0, lambda = 0.0;
    double x0 = 0.0, y0 = 0.0, z0 = 0.0;

private:
    struct Run
    {
        int offset, length;
    };

    int pml = 0;
    // Every field is one allocation with a distinct leading pad, so the many
    // same-index streams of a layer cell do not fall into the same cache
    // sets (the arrays are equal in size; unpadded they alias). Numerical
    // performance only.
    std::vector<Real> storage;
    Real* psiBase[3] = {};
    Real* auxBase[2] = {};
    Real* phiBase[2][3] = {};
    std::vector<CellType> types;
    std::vector<std::uint8_t> masks;
    std::vector<Run> runs;
    std::vector<int> runsOfPlane;   // runs of plane kz are [runsOfPlane[kz], runsOfPlane[kz + 1])
    std::vector<int> irregular;
    std::vector<int> interiorCells;
    // PML state (full-size arrays; only the layer is touched).
    // auxBase[n & 1]: psi = int u dt at level n; phiBase[n & 1][d]: phi^{n-1/2}
    // on the +d face of each cell.
    std::vector<double> zetaCell[3], zetaFace[3];
    // Face update coefficients (1 - k z/2) / (1 + k z/2) and k c^2 / (Delta (1 + k z/2)).
    std::vector<Real> faceDecay[3], faceGain[3];
    std::vector<Real> zetaCellR[3], zetaFaceR[3]; // the profiles in the air's precision
};

} // namespace taikor::candidates::fdtd
