#pragma once

// The air outside the drum: an axisymmetric (body-of-revolution) boundary
// element solution of the exterior Helmholtz problem of the whole rigid body
// (both heads, the rim annuli and the barrel), so that both heads load, damp
// and couple to each other through one exterior field.
//
// FORMULATION (e^{j w t}, G = e^{-j k R} / (4 pi R), n the body's outward normal,
// dp/dn = -j w rho v_n). For x in the air, the Kirchhoff-Helmholtz
// representation p(x) = int_S [p dG/dn_y - G dp/dn_y] dS_y; on a smooth surface
// point the left side is p/2; inside the body it is 0. A field of angular
// order m, p = p_m(r, z) cos(m theta), reduces it to the generatrix with the
// ring kernels
//   g_m(x, y) = (1/(4 pi)) int_0^{2 pi} cos(m phi) e^{-j k R} / R dphi,  h_m = d g_m / dn_y,
//   R^2 = r_x^2 + r_y^2 - 2 r_x r_y cos phi + (z_x - z_y)^2
// (A. F. Seybert, B. Soenarko, F. J. Rizzo & D. J. Shippy, "A special integral
// equation formulation for acoustic radiation and scattering for axisymmetric
// bodies and boundary conditions", JASA 80(4):1241, 1986; B. Soenarko, "A
// boundary element formulation for radiation of acoustic waves from
// axisymmetric bodies with arbitrary boundary conditions", JASA 93(2):631,
// 1993). The generatrix is cut into straight panels carrying a constant
// pressure, collocated at their midpoints:
//   p_i / 2 - sum_j A_ij p_j = j w rho sum S_i v,   A_ij = int_panel j h_m r dt,
// with the normal velocity v represented by its values at each panel's four
// Gauss nodes (a cubic per panel).
//
// IRREGULAR FREQUENCIES. The exterior equation is singular at the interior
// Dirichlet eigenfrequencies of the body. CHIEF rows (the representation at
// interior points, which must vanish) are appended and the system is solved by
// least squares (H. A. Schenck, "Improved integral formulation for acoustic
// radiation problems", JASA 44(1):41, 1968; Seybert et al. 1986).
//
// QUADRATURE (numerical; tests/modal bem_*): the phi integral by Gauss-Legendre
// (M + k sqrt(r_x r_y) + 16 nodes for orders up to M), with the sinh
// transformation phi = eps sinh(u) on [0, phi1] when the 1/R peak (width eps =
// d / sqrt(r_x r_y), d the points' distance in the (r, z) plane) is narrower
// than phi1 = 2 / (1 + M + k sqrt(r_x r_y)) (D. J. Johnston & D. Elliott, "A sinh
// transformation for evaluating nearly singular boundary element integrals",
// IJNME 62:564, 2005), all orders at once by the Chebyshev recurrence of
// cos(m phi), orders past the kernel's e^{-m eta} decay left at zero; the panel integral by
// four Gauss nodes when the point is far, by sinh-graded Gauss towards the
// nearest point when it is near, and by a cubic grading towards the midpoint
// on the panel itself (log singularity).

#include <array>
#include <complex>
#include <vector>

namespace taikor::candidates::modal
{

struct BemPanel
{
    double r0 = 0.0, z0 = 0.0, r1 = 0.0, z1 = 0.0;
    double length = 0.0, tr = 0.0, tz = 0.0, nr = 0.0, nz = 0.0; // unit tangent and outward normal
    int part = 0;                  // which surface it belongs to (drum: 0 batter, 1 rear, 2 rigid)
    [[nodiscard]] double r (double t) const noexcept { return r0 + t * tr; }
    [[nodiscard]] double z (double t) const noexcept { return z0 + t * tz; }
};

struct BemGeometry
{
    std::vector<BemPanel> panels;
    std::vector<std::array<double, 2>> chief; // interior points (r, z)
};

// A generatrix polyline from nodes (r, z), ordered so that the body lies to its
// left (n = (-t_z, t_r) points into the air); segment s belongs to parts[s] and
// is cut into panels no longer than steps[s].
[[nodiscard]] BemGeometry polylineGeometry (const std::vector<std::array<double, 2>>& nodes, const std::vector<int>& parts,
                                            const std::vector<double>& steps, const std::vector<std::array<double, 2>>& chief);

// The drum body: batter head (part 0) at z = 0 for r < a, rim annulus a..b,
// barrel r = b down to z = -L, rim, rear head (part 1) at z = -L; b = a + wall.
// Panels on the heads no longer than headStep, elsewhere bodyStep; CHIEF points
// spread through the interior (a fixed low-discrepancy set).
[[nodiscard]] BemGeometry drumGeometry (double a, double wall, double length, double headStep, double bodyStep,
                                        int chiefPoints);

// Ring kernels g_m and h_m (as defined above) for m = 0..maxOrder at a field
// point (rx, zx) and a source point (ry, zy) with normal (nr, nz).
void ringKernels (double rx, double zx, double ry, double zy, double nr, double nz, double k, int maxOrder, int refine,
                  std::complex<double>* g, std::complex<double>* h);

// The four Gauss nodes of a panel at which normal velocity is given.
inline constexpr int bemVelocityNodes = 4;
[[nodiscard]] std::array<double, bemVelocityNodes> bemNodeOffsets (const BemPanel&); // t along the panel

// One order's problem: normal velocity columns given at the velocity nodes of
// the panels of the listed parts ("active" parts), node-major:
// velocity[(panelIndexAmongActive * 4 + node) * columns + column].
struct BemProblem
{
    int m = 0;
    bool chief = true;             // append the CHIEF rows (least squares); false: the square system (LU)
    std::vector<int> activeParts;
    int columns = 0;
    std::vector<double> velocity;
    // Output: surface pressure per panel and column, [panel * columns + column],
    // for p in pascals per unit normal velocity (m/s).
    std::vector<std::complex<double>> pressure;
    // Optional field points (r, z) (at azimuth 0 of the cos(m theta) pattern):
    // pressure there, [point * columns + column]. The first problem's points
    // are used for every problem of one solve.
    std::vector<std::array<double, 2>> fieldPoints;
    std::vector<std::complex<double>> field;
};

class BemSolver
{
public:
    BemSolver (BemGeometry geometry, int refine = 1);
    // Solves every problem (each its own order) at wavenumber k; kernels are
    // computed once for all orders. useChief = false solves the square system.
    void solve (double k, double density, double soundSpeed, std::vector<BemProblem>& problems, bool useChief = true,
                int threads = 1) const;
    [[nodiscard]] const BemGeometry& geometry() const noexcept { return geo; }

private:
    BemGeometry geo;
    int refine;
};

} // namespace taikor::candidates::modal
