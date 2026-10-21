#pragma once

// Energy-consistent coupling of a head to the air planes on its two sides
// (Bilbao, JASA 131(1), 2012, eqs. 21-22 and 39-40; Bilbao & Webb, JAES
// 61(10), 2013, eqs. 21-22; Torin, PhD thesis 2016, eqs. 5.47-5.59):
//
//   air cell above:  Psi+^{n+1} = Psi^+ + (c^2 k^2 / Delta) J d_t. w
//   air cell below:  Psi-^{n+1} = Psi^- - (c^2 k^2 / Delta) J d_t. w
//   head node:       force h^2 rho I (d_t. Psi- - d_t. Psi+)
//
// where Psi^ is the leapfrog value with the face across the head dropped.
// I interpolates the air plane (cell centres of the columns inside the head,
// r_c < a) to each membrane node with quadratic B-spline weights (see
// Coupling.cpp for why not bilinear); J = (h^2 / Delta^2) I^T. That J is the
// scaled transpose of I is exactly what makes the air's energy gain equal the
// head's energy loss (Torin eq. 5.51), for any I. Near the rim a node's
// weights on columns outside the head are dropped and the rest renormalised
// to sum to one, so the air sees the head's whole volume velocity.
//
// Eliminating Psi^{n+1} gives the head update (alpha I + mu E) w^{n+1} = ...
// with E = I J, mu = rho c^2 k^2 / (2 sigma Delta), alpha = 1 + sigma0 k.
// It is solved in the air-plane space by the Woodbury identity:
//   (alpha I + mu I J)^-1 v = (v - nu I z) / alpha,  (I_a + nu G) z = J v,
// with G = J I (small, 25 non-zeros per row) and nu = mu / alpha ~ 6e-4 at
// 96 kHz, and (I_a + nu G)^-1 is the Neumann series sum_p (-nu G)^p to the
// order FdtdNumerics::neumannOrder, stored as one sparse matrix. Also
// J (A^-1 v) = z / alpha, which the air-plane update uses directly.

#include "HeadGrid.h"

#include <array>
#include <vector>

namespace taikor::candidates::fdtd
{

struct SparseMatrix
{
    int rows = 0, cols = 0;
    std::vector<int> start;    // rows + 1
    std::vector<int> column;
    std::vector<double> value;

    // y[r] = sum_j A[r, j] x[j] for r in [r0, r1).
    void multiply (const double* x, double* y, int r0, int r1) const noexcept
    {
        for (int r = r0; r < r1; ++r)
        {
            double sum = 0.0;
            for (int e = start[static_cast<std::size_t> (r)]; e < start[static_cast<std::size_t> (r) + 1]; ++e)
                sum += value[static_cast<std::size_t> (e)] * x[column[static_cast<std::size_t> (e)]];
            y[r] = sum;
        }
    }
    [[nodiscard]] int nonZeros() const noexcept { return static_cast<int> (value.size()); }
};

// A B for CSR matrices.
[[nodiscard]] SparseMatrix sparseProduct (const SparseMatrix& a, const SparseMatrix& b);
[[nodiscard]] SparseMatrix sparseTranspose (const SparseMatrix& a, double scale);

struct AirPlaneLayout
{
    double x0 = 0.0, y0 = 0.0, delta = 0.0; // cell centres at x0 + (i + 1/2) delta
    int nx = 0, ny = 0;
};

class HeadCoupling
{
public:
    // columns: the air columns (i, j) inside the head; their index in this
    // list is the air-plane index c used by every vector below.
    void build (const HeadGrid& head, const AirPlaneLayout& layout,
                const std::vector<std::array<int, 2>>& columns, double nu, int neumannOrder);

    int planeSize() const noexcept { return static_cast<int> (columnList.size()); }
    const std::vector<std::array<int, 2>>& columns() const noexcept { return columnList; }

    // The barrel's wall loss scales the head's flux into some plane cells
    // (FdtdEngine.h HeadState::sliceScale), so the head solves
    // (alpha I + mu I S J), S diagonal: (alpha I + mu I S J)^-1 v =
    // (v - nu I S z) / alpha with (I_a + nu G S) z = J v, and still
    // J (A^-1 v) = z / alpha. This rebuilds `inverse` as the Neumann series of
    // (I_a + nu G S)^-1 (S differs from 1 by ~1e-4 on the wall's ring).
    void setSliceScale (const std::vector<double>& scale, int neumannOrder);

    SparseMatrix interpolate;   // I: membrane node (full grid index) <- air plane
    SparseMatrix spread;        // J = (h^2 / Delta^2) I^T: air plane <- membrane node
    SparseMatrix gram;          // G = J I
    SparseMatrix inverse;       // (I_a + nu G)^-1, truncated Neumann series
    double transposeScale = 0.0; // Delta^2 / h^2, so I^T = transposeScale * J
    double nu = 0.0;

private:
    void buildInverse (const std::vector<double>* scale, int neumannOrder);
    std::vector<std::array<int, 2>> columnList;
};

} // namespace taikor::candidates::fdtd
