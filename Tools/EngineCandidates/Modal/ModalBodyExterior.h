#pragma once

// The exterior loading of the coupled band from the boundary-element solution
// of the whole drum body (ModalBem.h): per angular order m and frequency, the
// complex radiation impedance over the block's head modes, batter then rear,
//   Z_ij(w) = Lam_m int_heads p^(j)(r) phi_i(r) r dr,
// p^(j) the surface pressure for the outward normal velocity phi_j on its own
// head (zero elsewhere). Im Z / w is the added mass and Re Z the resistance,
// and both include the batter-rear terms that one shared exterior field
// carries: the two heads' volume velocities of a push-pull mode cancel in the
// far field instead of radiating as two independent baffled pistons.
//
// It covers the head modes below maxHz (the "columns": per order the first
// batterColumns(m) batter and rearColumns(m) rear modes, both ascending in
// frequency) and is tabulated on a global log frequency grid (gridPerOctave
// points per octave) from below half the lowest of them up to maxHz; each
// order's block (ModalBlock.cpp) overlays it on the per-head baffled King
// loading (ModalExterior.h) at the grid points in its range and interpolates
// the resistance to each mode's own frequency (cubic in log w). Pairs with a
// mode above maxHz, and all frequencies above it, keep the King loading: such
// a mode is subsonic with a near field that decays within 1/b1 < 1.5 cm of the
// head, inside the flat 52 mm rim annulus, so the finite body barely changes
// its loading (tests/modal body_exterior_convergence: maxHz, panel size and grid
// density against refined values).
//
// CAPSULES. The same solutions give the pressure at given field points (the
// capsules) per unit normal velocity of each column; the table keeps it minus
// a reference transfer supplied by the caller (the engine's own baffled
// Rayleigh / rear-monopole transfer at the same frequency), i.e. what the drum
// body adds at the capsule: diffraction at the rim corner and the barrel, and
// the rear head's field around the body (ModalMics.cpp uses it).

#include "ModalBem.h"
#include "ModalCavity.h"
#include "ModalMath.h"

#include <array>
#include <complex>
#include <functional>
#include <memory>
#include <vector>

namespace taikor::candidates::modal
{

struct BodyExteriorOptions
{
    double maxHz = 3900.0;            // top of the table
    int gridPerOctave = 12;
    double panelsPerWavelength = 10.0;// acoustic wavelength at maxHz, and the shortest head-mode wavelength
    int chiefPoints = 40;
    int threads = 8;
    // Field points (r, z), at azimuth 0 of each order's cos(m theta) pattern, and
    // the reference pressure per unit normal velocity of column c of order m at w.
    std::vector<std::array<double, 2>> fieldPoints;
    std::function<std::complex<double> (int m, int column, double omega, int point)> reference;
};

// The global grid: w_g = 2 pi 2^(g / gridPerOctave) rad/s for every integer g.
[[nodiscard]] double bodyGridOmega (int g, int gridPerOctave);

class BodyExterior
{
public:
    // batter / rear: the head modes of each order below maxHz (the columns), per m.
    BodyExterior (const AirSpec& air, double wallThickness, const std::vector<std::vector<const HeadMode*>>& batter,
                  const std::vector<std::vector<const HeadMode*>>& rear, const BodyExteriorOptions& options);

    [[nodiscard]] const std::vector<double>& grid() const noexcept { return omegas; }
    [[nodiscard]] int gridDensity() const noexcept { return gridPerOctave; }
    [[nodiscard]] double maxOmega() const noexcept { return omegas.empty() ? 0.0 : omegas.back(); }
    // Z (row-major nh x nh, nh = batter[m].size() + rear[m].size()) of order m at grid point g;
    // empty if the order has no head mode or the point is outside its range.
    [[nodiscard]] const std::vector<std::complex<double>>& impedance (int m, int g) const;
    // Added mass Im Z / w and resistance Re Z at any w within the table (cubic in log w).
    void addedMass (int m, double omega, std::vector<double>& out) const;
    void resistance (int m, double omega, std::vector<double>& out) const;
    // Field-point pressure minus the reference, per unit normal velocity of each
    // column, [point * nh + column], at any w within the table (cubic in log w).
    void fieldCorrection (int m, double omega, std::vector<std::complex<double>>& out) const;
    [[nodiscard]] int panels() const noexcept { return panelCount; }
    [[nodiscard]] int batterColumns (int m) const noexcept { return m >= 0 && static_cast<std::size_t> (m) < batterCount.size() ? batterCount[static_cast<std::size_t> (m)] : 0; }
    [[nodiscard]] int rearColumns (int m) const noexcept { return m >= 0 && static_cast<std::size_t> (m) < rearCount.size() ? rearCount[static_cast<std::size_t> (m)] : 0; }
    [[nodiscard]] double seconds() const noexcept { return buildSeconds; }
    // Largest |Z_ij - Z_ji| / sqrt(|Z_ii Z_jj|) before symmetrisation.
    [[nodiscard]] double maxAsymmetry() const noexcept { return asymmetry; }

private:
    std::vector<double> omegas;
    std::vector<std::vector<std::vector<std::complex<double>>>> table; // [m][g]
    std::vector<std::vector<std::vector<std::complex<double>>>> fields; // [m][g], [point * nh + column]
    std::vector<int> heads, batterCount, rearCount;                   // columns per m
    int firstIndex = 0, gridPerOctave = 12;                           // omegas[i] = bodyGridOmega (firstIndex + i)
    int panelCount = 0;
    double buildSeconds = 0.0, asymmetry = 0.0;
    void interpolate (int m, double omega, bool real, std::vector<double>& out) const;
    // Cubic Lagrange weights in log w over grid points [first, last] of order m's table (empty: none).
    bool weights (int m, double omega, int& first, std::vector<double>& w) const;
};

} // namespace taikor::candidates::modal
