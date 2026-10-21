#pragma once

// Shared pieces of the finite-difference heads: the grid a head lives on and
// the interface through which its geometric nonlinearity plugs into the
// scalar-auxiliary-variable (SAV) time step. The Berger and von Karman letters
// differ in exactly one object, the NonlinearPotential, and nothing else.

#include "CandidateEngine.h"

#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace taikor::candidates
{

// A clamped circular head on a Cartesian grid. Node (i, j) sits at
// (x0 + i h, y0 + j h) relative to the head centre, index j * nx + i. Nodes
// flagged interior are updated; every other node is held at w = 0, and the
// grid carries at least two clamped rings outside the rim so the 13-point
// biharmonic stencil never reads past the array.
struct MembraneGrid
{
    int nx = 0, ny = 0;
    double h = 0.0;
    double x0 = 0.0, y0 = 0.0;
    std::vector<std::uint8_t> interior;

    [[nodiscard]] int index (int i, int j) const noexcept { return j * nx + i; }
    [[nodiscard]] double x (int i) const noexcept { return x0 + i * h; }
    [[nodiscard]] double y (int j) const noexcept { return y0 + j * h; }
    [[nodiscard]] int size() const noexcept { return nx * ny; }
};

// The staircase disk every finite-difference head uses: node spacing h,
// interior where r < inclusionRadius. inclusionRadius is normally the head
// radius; a scheme may move it by a fraction of h to correct the staircase's
// systematic flattening of the mode frequencies, which is a numerical
// correction and must be documented as one where it is applied.
[[nodiscard]] inline MembraneGrid makeDiskGrid (double inclusionRadius, double h)
{
    MembraneGrid grid;
    grid.h = h;
    const int half = static_cast<int> (std::ceil (inclusionRadius / h)) + 2;
    grid.nx = grid.ny = 2 * half + 1;
    grid.x0 = grid.y0 = -half * h;
    grid.interior.assign (static_cast<std::size_t> (grid.size()), 0);
    for (int j = 0; j < grid.ny; ++j)
        for (int i = 0; i < grid.nx; ++i)
        {
            const double r = std::hypot (grid.x (i), grid.y (j));
            grid.interior[static_cast<std::size_t> (grid.index (i, j))] =
                r < inclusionRadius ? 1 : 0;
        }
    return grid;
}

// A conservative geometric nonlinearity of one head: a potential energy
// V(w) >= 0 of the discrete displacement field, stored in joules, and its
// gradient dV/dw_l at every node in newtons. The head's nodal equation is
//     sigma h^2 d_tt w_l = (linear terms) - dV/dw_l,
// and the SAV scheme (Bilbao, Ducceschi & Zama, "Explicit exactly
// energy-conserving methods for Hamiltonian systems", J. Comput. Phys. 472,
// 2023) needs exactly these two things at each step.
class NonlinearPotential
{
public:
    virtual ~NonlinearPotential() = default;
    [[nodiscard]] virtual std::string name() const = 0;
    virtual bool prepare (const MembraneGrid& grid, const HeadSpec& head,
                          double headRadius, std::string& error) = 0;
    // w holds grid.size() values, zero at clamped nodes. Returns V(w) and
    // writes dV/dw into gradient (grid.size() values, zero at clamped nodes).
    // Called once per time step with the current displacement.
    virtual double evaluate (const double* w, double* gradient) = 0;
};

// Berger (1955) average stretching with an immovable edge:
// V = xi / 4 (int |grad w|^2 dA)^2, xi = E h / (2 (1 - nu^2) pi a^2).
[[nodiscard]] std::unique_ptr<NonlinearPotential> makeBergerPotential();
// Foppl-von Karman membrane stretching with an immovable (tacked) edge:
// V = min_u 1/2 int N : eps dA, eps = sym grad u + 1/2 grad w (x) grad w,
// N = E h / (1 - nu^2) [(1 - nu) eps + nu tr(eps) I], u = 0 on the rim.
[[nodiscard]] std::unique_ptr<NonlinearPotential> makeVonKarmanPotential();

} // namespace taikor::candidates
