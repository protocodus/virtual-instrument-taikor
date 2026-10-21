#pragma once

// One stiff, lossy, clamped head on a Cartesian staircase disk: the linear
// part of the finite-difference scheme (Bilbao, Numerical Sound Synthesis,
// Wiley 2009, ch. 11-12; Bilbao & Webb, "Physical modeling of timpani drums in
// 3D on GPGPUs", JAES 61(10), 2013, eq. 14; Torin, PhD thesis, Edinburgh 2016,
// eqs. 2.165-2.168):
//
//   sigma d_tt w = T Lap w - D Lap^2 w - 2 sigma sigma0 d_t w
//                  + 2 sigma sigma1 Lap d_t w + (pressure jump) - dV/dw - (stick)
//
// discretised per node l in newtons (the nodal equation is multiplied by the
// node's area h^2):
//
//   sigma h^2 d_tt w = T h^2 L w - D h^2 L^T L w - 2 sigma sigma0 h^2 d_t. w
//                      + 2 sigma sigma1 h^2 L d_t- w + ...
//
// L is the 5-point Laplacian with w = 0 on every non-interior node; applied
// twice (the second time on the interior nodes only, reading the first
// clamped ring) it is exactly the 13-point biharmonic with two zero rings,
// i.e. the clamped plate operator, and L^T L makes the bending energy a sum
// of squares. The sigma1 term uses the backward difference d_t- so the update
// stays explicit; its energy analysis (a mass operator sigma (1 - sigma1 k
// (-L)), positive under the stability bound) is in energy() below.

#include "FdtdContracts.h"

#include <vector>

namespace taikor::candidates::fdtd
{

// Stability bound of the explicit stiff lossy membrane (Bilbao-Webb 2013
// eq. 23 with the sigma1 term; Torin 2016 eq. 2.168):
//   h >= sqrt(c^2 k^2 + 4 sigma1 k + sqrt((c^2 k^2 + 4 sigma1 k)^2 + 16 kappa^2 k^2)),
// c^2 = (T + extraTension) / sigma, kappa^2 = D / sigma. extraTension is a
// tension the scheme applies explicitly on top of T (the SAV step's
// nonlinear membrane force, FdtdDrumConfig.h headSpacing).
[[nodiscard]] double stableHeadSpacing (const HeadSpec& head, double k, double extraTension = 0.0) noexcept;

class HeadGrid
{
public:
    // spacing >= stableHeadSpacing; inclusionRadius ~ radius (see makeDiskGrid).
    void build (const HeadSpec& head, double radius, double k, double spacing,
                double inclusionRadius, bool lossless);

    const MembraneGrid& grid() const noexcept { return mesh; }
    const HeadSpec& spec() const noexcept { return head; }
    int size() const noexcept { return mesh.size(); }
    int rows() const noexcept { return mesh.ny; }
    int interiorCount() const noexcept { return interior; }
    double spacing() const noexcept { return h; }
    double timeStep() const noexcept { return k; }
    double radius() const noexcept { return headRadius; }
    double sigma0() const noexcept { return s0; }
    double sigma1() const noexcept { return s1; }
    // 1 + sigma0 k: the coefficient of w^{n+1} in the lossy update.
    double alpha() const noexcept { return 1.0 + s0 * k; }

    // Interior nodes of row j are [rowBegin, rowEnd); the interior plus the
    // first clamped ring of row j is [extBegin, extEnd).
    int rowBegin (int j) const noexcept { return interiorBegin[static_cast<std::size_t> (j)]; }
    int rowEnd (int j) const noexcept { return interiorEnd[static_cast<std::size_t> (j)]; }
    int extBegin (int j) const noexcept { return extendedBegin[static_cast<std::size_t> (j)]; }
    int extEnd (int j) const noexcept { return extendedEnd[static_cast<std::size_t> (j)]; }
    bool isInterior (int index) const noexcept { return mesh.interior[static_cast<std::size_t> (index)] != 0; }

    // lap = L w on rows [j0, j1) over the extended segments (w = 0 outside).
    void laplacian (const double* w, double* lap, int j0, int j1) const noexcept;

    // The explicit part of the update on interior rows [j0, j1):
    //   out = 2 w - (1 - sigma0 k) wPrev + k^2/sigma (T lap - D L lap)
    //         + 2 sigma1 k (lap - lapPrev),
    // i.e. alpha w^{n+1} = out + (coupling, nonlinear and contact terms).
    // lap must hold L w on the extended segments of rows [j0 - 1, j1 + 1);
    // lapPrev must hold L wPrev on the interior of rows [j0, j1).
    void explicitPart (const double* w, const double* wPrev, const double* lap,
                       const double* lapPrev, double* out, int j0, int j1) const noexcept;

    // Discrete energy at time n+1/2 from w^{n+1} (next) and w^n (now), in J:
    // kinetic (with the sigma1 mass correction) + tension + bending.
    [[nodiscard]] double energy (const double* next, const double* now) const noexcept;
    // Power dissipated by the hide loss at time n (J/s), from w^{n+1}, w^{n-1}.
    [[nodiscard]] double lossPower (const double* next, const double* prev) const noexcept;

    // Sums over interior nodes.
    [[nodiscard]] double dot (const double* a, const double* b) const noexcept;

private:
    MembraneGrid mesh;
    HeadSpec head;
    double headRadius = 0.0, k = 0.0, h = 0.0, s0 = 0.0, s1 = 0.0;
    int interior = 0;
    std::vector<int> interiorBegin, interiorEnd, extendedBegin, extendedEnd;
};

// In-vacuo stepping of one head (its own loss setting, no air, no
// nonlinearity) from rest at an initial displacement shape (w^0 = w^{-1} =
// shape), for eigenfrequency measurements. Returns, per probe, the trace
// s_n = <w^n, probe> for n = 0..steps-1.
[[nodiscard]] std::vector<std::vector<double>> simulateInVacuo (
    const HeadGrid& head, const std::vector<double>& shape,
    const std::vector<std::vector<double>>& probes, int steps);

// Frequency (Hz) of a signal dominated by one mode of a linear recurrence:
// cos(omega k) = sum s_n (s_{n+1} + s_{n-1}) / (2 sum s_n^2).
[[nodiscard]] double pronyFrequency (const std::vector<double>& s, double k) noexcept;

} // namespace taikor::candidates::fdtd
