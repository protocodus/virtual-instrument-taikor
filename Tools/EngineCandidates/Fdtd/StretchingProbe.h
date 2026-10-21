#pragma once

// An optional second interface of a NonlinearPotential (FdtdContracts.h): what
// the finite-difference drum needs to know about a head's stretching beyond V
// and dV/dw, to size the head's grid and to report on the take. The shared
// contract stays as it is; the drum asks for this one with a dynamic_cast, and
// a potential without it is treated as unbounded (the drum then cannot size a
// grid for it and says so).
//
// Why the drum needs it. The SAV step applies the stretching force
// explicitly, -g psi with g = grad V / sqrt(2 V + eps), so for a perturbation
// d of the head it acts, with psi tracking sqrt(2 V), like the frozen
// stress-stiffness operator G(N): d^T G d = int N : (grad d (x) grad d) dA.
// The explicit scheme stays within its stability bound only if that extra
// stiffness is covered by the grid spacing (FdtdDrumConfig.h headSpacing,
// HeadGrid.h stableHeadSpacing); beyond it the grid's top modes can grow at
// almost no energy cost and psi loses track of sqrt(2 V) (tests/fdtd
// savTracking; tests/integration vonKarmanTracking). The measure that bound
// needs is the smallest N* with d^T G d <= N* d^T K_1 d for every d, K_1 the
// head's unit-tension membrane operator (d^T K_1 d = sum over edges of the
// squared differences): the grid is then stable if it is stable at the
// uniform tension T + N*.

#include <string>

namespace taikor::candidates::fdtd
{

class StretchingProbe
{
public:
    virtual ~StretchingProbe() = default;

    // N* (N/m) of the stress field of the last evaluate(), as above; 0 when
    // the head is flat. Called after evaluate() on the same thread.
    [[nodiscard]] virtual double peakStretchingTension() const = 0;

    // True when N* is bounded a priori by the energy the drum holds, so the
    // grid can be sized before the stroke (Berger: the uniform tension
    // xi int |grad w|^2 <= 2 sqrt(xi E), FdtdDrumConfig.h bergerTensionBound).
    // False when only a run can tell (von Karman: the stress concentrates
    // under the stick; the energy bounds only its mean square).
    [[nodiscard]] virtual bool boundedByEnergy() const = 0;

    // One line of the potential's own statistics for the drum's diagnostics
    // (solver work and accuracy), or empty.
    [[nodiscard]] virtual std::string report() const = 0;
};

// Another optional interface: a potential whose evaluate a team of the drum's
// threads can share (von Karman: the in-plane solve, VonKarman/Team.h). The
// drum gives the team one P1 work item per member (FdtdEngineStep.cpp); the
// items start together at the head of the queue, and the other threads take
// the air meanwhile.
class TeamPotential
{
public:
    virtual ~TeamPotential() = default;
    // Teams of `size` members from now on (1: the serial evaluate()); returns
    // the size in effect (the potential may cap it). Not concurrently with an
    // evaluate.
    virtual int setTeam (int size) = 0;
    // Member `member`'s share of evaluate (w, gradient); every member of the
    // team calls it with the same w and gradient, concurrently, and gets V.
    virtual double evaluateMember (const double* w, double* gradient, int member) = 0;
};

} // namespace taikor::candidates::fdtd
