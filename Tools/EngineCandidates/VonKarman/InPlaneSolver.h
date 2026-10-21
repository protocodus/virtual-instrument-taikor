#pragma once

// The quasi-static in-plane problem of a von Karman membrane with an immovable
// edge, discretised on the staggered (MAC-type) grid that the transverse
// displacement w lives on, and the multigrid-preconditioned conjugate-gradient
// solver used to minimise the stretching energy over the in-plane
// displacement every time step.
//
// Grid. The w nodes of the MembraneGrid are the cell centres of the in-plane
// grid: u (x displacement) lives on x-edges (i + 1/2, j), v on y-edges
// (i, j + 1/2), normal strains at the nodes and the engineering shear at the
// cell corners (i + 1/2, j + 1/2). Every array is stored as an n x n square
// with index k = j n + i (u at k is the edge east of node k, v at k the edge
// north of node k, the shear at k the corner north-east of node k), and a
// level vector concatenates [pad | u | pad | v | pad] so that every stencil
// shift is a plain offset and reads zeros outside the square.
//
// Immovable edge. u (resp. v) is free only if both nodes it joins are
// interior; every other edge displacement is held at zero. This places the
// in-plane Dirichlet condition on the same staircase rim as w's: an edge
// crossing the rim carries a slope of w but no in-plane displacement. Grid
// convergence of the stretching energy against exact continuum values with
// u(a) = 0 (tests/vonkarman) chose this rule over "either node interior" and
// "edge midpoint inside r < a", both of which converge an order slower for a
// field with a smooth edge.
//
// Precision. The conjugate-gradient iteration, its operator and the returned
// field are double; the multigrid preconditioner runs in single precision
// (it only has to approximate K^-1), which halves its memory traffic and
// doubles its SIMD width. The CG uses the flexible (Polak-Ribiere) beta, so
// the preconditioner's rounding-level asymmetry cannot stall it.

#include "Team.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace taikor::candidates::vk
{

// Plane-stress membrane constants, N/m: C11 = E h / (1 - nu^2), C12 = nu C11,
// C66 = (1 - nu) C11 / 2 (engineering shear strain).
struct PlaneStress
{
    double c11 = 0.0, c12 = 0.0, c66 = 0.0;
};

// One level of the in-plane hierarchy.
struct InPlaneLevel
{
    int n = 0;          // square side
    int N = 0;          // n * n entries per component
    int pad = 0;        // zero margin around each component block
    std::vector<std::uint8_t> freeU, freeV; // N each

    [[nodiscard]] long size() const noexcept { return 3L * pad + 2L * N; }
    [[nodiscard]] long base (int component) const noexcept
    {
        return component == 0 ? pad : 2L * pad + N;
    }

    // The finest level (and rediscretised coarse levels) apply K = D^T C D
    // matrix-free. A Galerkin coarse level stores its operator in two parts:
    // the constant stencil of the rows away from the rim, applied to every
    // row, and the explicit rows near the rim, which overwrite it.
    bool matrixFree = false;
    struct Tap
    {
        long offset = 0;          // column minus row, level-vector index
        double coefficient = 0.0;
    };
    std::vector<Tap> stencil[2];           // per row component
    std::vector<float> stencilF[2];        // single-precision copies for the preconditioner
    // Rim rows: corrections to the stencil in slot-major ELL storage (slot s
    // of rim row r at s * rows + r; unused slots have offset 0, coefficient 0).
    std::vector<long> irregularRow;        // level-vector index of each rim row
    int irregularWidth = 0;
    std::vector<long> irregularOffset;
    std::vector<double> irregularCoefficient;
    std::vector<float> irregularCoefficientF;

    std::vector<double> inverseDiagonal;   // level-vector sized, zero off free dofs
    std::vector<float> inverseDiagonalF;
    double lambdaMax = 0.0;                // upper bound on the spectrum of D^-1 K
};

// Marks the coarse dofs the prolongation from fine reaches (InPlaneOperator.cpp).
void markCoarseFree (const InPlaneLevel& fine, InPlaneLevel& coarse);

class InPlaneSolver
{
public:
    // Every value here is numerical, not physical, and tests/vonkarman backs
    // each one: the tolerance by the output's insensitivity to tightening it
    // ("tolerance"), the smoother degrees and the Galerkin coarse operators by
    // the measured convergence rate and cost ("solver"); the others only
    // bound work (the coarsest level is solved exactly whatever its size).
    struct Settings
    {
        double relativeTolerance = 3.0e-3; // ||b - K x||_2 <= tol ||b||_2
        int maxIterations = 200;           // PCG safety cap
        int chebyshevDegree = 3;           // smoother polynomial degree per sweep, finest level
        int coarseChebyshevDegree = 1;     // the same on the coarse levels
        int coarsestMaxFree = 160;         // exact dense solve below this many dofs
        bool galerkinCoarse = true;        // P^T K P coarse operators (else rediscretised)
    };

    // interior: n x n node flags (the MembraneGrid's), with at least two
    // clamped rings round the rim, as makeDiskGrid guarantees.
    bool build (int n, const std::vector<std::uint8_t>& interior,
                const PlaneStress& stiffness, const Settings& settings,
                std::string& error);

    [[nodiscard]] const InPlaneLevel& fine() const noexcept { return levels_.front(); }
    [[nodiscard]] int levelCount() const noexcept { return static_cast<int> (levels_.size()); }
    [[nodiscard]] const InPlaneLevel& level (int l) const noexcept { return levels_[static_cast<std::size_t> (l)]; }
    [[nodiscard]] const PlaneStress& stiffness() const noexcept { return c_; }
    [[nodiscard]] Settings& settings() noexcept { return settings_; }
    [[nodiscard]] const Settings& settings() const noexcept { return settings_; }

    // y = K x on level l (x zero off free dofs; y is written zero there too),
    // in double or single precision.
    void apply (int l, const double* x, double* y) const;
    void apply (int l, const float* x, float* y) const;

    // Solve K x = b on the finest level by PCG from the initial guess in x,
    // to ||b - K x|| <= tolerance ||b|| (settings().relativeTolerance when
    // tolerance <= 0). Returns the iteration count; relativeResidual receives
    // ||b - K x|| / ||b||.
    int solve (const double* b, double* x, double& relativeResidual, double tolerance = 0.0);
    // The same by a team of threads (Team.h; team.size > 1), every member
    // calling with the same arguments. The members split the finest two
    // levels and the conjugate-gradient passes by rows; the coarser levels,
    // a few percent of the work, run on member 0 (teamLevels). x must be
    // complete on entry (the caller syncs), and is complete on return.
    int solve (const double* b, double* x, double& relativeResidual, double tolerance, Team& team);

    // One single-precision multigrid V-cycle z = B r on level l.
    void vcycle (int l, const float* r, float* z);
    // The same by a team (Team.h), every member calling with the same r and
    // z: entry for entry the serial V-cycle's arithmetic, so the same z bit
    // for bit (tests/integration teamEvaluate). r must be complete on entry;
    // z is complete on return.
    void vcycle (int l, const float* r, float* z, Team& team);
    // The preconditioner on the finest level in double: z = B r.
    void precondition (const double* r, double* z);

    // Rayleigh quotient of D^-1 K after a power iteration (for tests).
    [[nodiscard]] double rayleighEstimate (int l, int iterations) const;

    // Relative residual of the initial guess of the last solve.
    [[nodiscard]] double lastInitialResidual() const noexcept { return lastInitialResidual_; }

    // Grid transfers between level coarse - 1 and coarse (public for tests).
    void prolong (int coarse, const double* xc, double* xf) const;   // xf = P xc (overwrites)
    void restrict_ (int coarse, const double* rf, double* rc) const; // rc = P^T rf (overwrites)
    void prolong (int coarse, const float* xc, float* xf) const;
    void restrict_ (int coarse, const float* rf, float* rc) const;

private:
    template <typename T> void applyT (int l, const T* x, T* y) const;
    // y = K x on the rows [j0, j1) of both components (member's scratch).
    template <typename T> void applyRows (int l, const T* x, T* y, int j0, int j1, int member) const;
    // Fine rows [j0, j1) of xf = P xc; coarse rows [j0, j1) of rc = P^T rf.
    template <typename T> void prolongRows (int coarse, const T* xc, T* xf, int j0, int j1) const;
    template <typename T> void restrictRows (int coarse, const T* rf, T* rc, int j0, int j1) const;
    void chebyshev (int l, const float* b, float* x, bool zeroInitial, Team& team);
    // Levels from this one down run on member 0 alone: level 2 of the
    // nagado batter is 45 x 45 (1/14 of the finest level's rows), where a
    // team's barriers would cost more than its split saves. Numerical.
    static constexpr int teamLevels = 2;
    template <typename T> void prolongT (int coarse, const T* xc, T* xf) const;
    template <typename T> void restrictT (int coarse, const T* rf, T* rc) const;
    void addWork (const InPlaneLevel& L);
    bool buildCoarse (int fineLevel, std::string& error);
    bool buildCoarsest (std::string& error);
    void estimateLambdaMax (int l);
    void chebyshev (int l, const float* b, float* x, bool zeroInitial);

    PlaneStress c_;
    Settings settings_;
    std::vector<InPlaneLevel> levels_;

    // Per-level single-precision work vectors (level-vector sized) and
    // matrix-free stress scratch in both precisions.
    struct Work
    {
        std::vector<float> b, x;    // this level's V-cycle right-hand side and correction
        std::vector<float> t, d, c; // residual / K d, Chebyshev direction, Chebyshev residual
        mutable std::vector<float> sxxF, syyF, sxyF;
        mutable std::vector<double> sxx, syy, sxy;
        mutable std::vector<float> rimF;   // rim-row accumulators
        mutable std::vector<double> rim;
    };
    std::vector<Work> work_;
    std::vector<double> pcgR_, pcgZ_, pcgP_, pcgQ_;
    std::vector<float> zOldF_;
    double lastInitialResidual_ = 0.0;

    // Coarsest level: dense Cholesky of the free dofs and the explicit inverse.
    std::vector<long> coarsestIndex_;
    std::vector<double> coarsestFactor_;
    std::vector<double> coarsestRhs_;
    std::vector<float> coarsestInverse_; // column-major
    std::vector<float> coarsestWork_;
};

} // namespace taikor::candidates::vk
