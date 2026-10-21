#pragma once

// The finite-difference drum: two stiff lossy heads (HeadGrid) with their
// geometric nonlinearity through the scalar-auxiliary-variable (SAV) step,
// embedded in 3-D air (AirDomain) through energy-consistent coupling
// (Coupling.h), struck by a lumped bachi tip through an energy-conserving
// Hertz collision (Collision.h), heard through two pressure capsules.
//
// Per head, per step (nodal equation divided by sigma h^2, times k^2):
//   alpha w^{n+1} + mu E w^{n+1} + beta g (g . w^{n+1})
//     = rExpl + k^2 rho / sigma I (q- - q+) + mu E w^{n-1}
//       - gamma g psi^{n-1/2} + beta g (g . w^{n-1}) - gamma c f,
// with rExpl the explicit linear part (HeadGrid::explicitPart),
// q+- = (Psi^+- - Psi+-^{n-1}) / 2k the air planes' partial centred
// differences, g = grad V / sqrt(2 V + eps) and psi the SAV variable
// (S. Bilbao, M. Ducceschi & F. Zama, "Explicit exactly energy-conserving
// methods for Hamiltonian systems", J. Comput. Phys. 472, 2023, 111697,
// eqs. 26-29 and 34: psi^{n+1/2} = psi^{n-1/2} + g . (w^{n+1} - w^{n-1}) / 2),
// beta = k^2 / (4 sigma h^2), gamma = k^2 / (sigma h^2), c the footprint
// weights (sum 1) and f the contact force (batter only). The matrix is
// A + beta g g^T with A = alpha I + mu E: Sherman-Morrison for the rank-1
// SAV term, Woodbury in the air-plane space for A (Coupling.h), and the
// collision reduces to one scalar Newton solve (Collision.h).
//
// Threads: every step runs as four phases separated by barriers:
//   P1  (dynamic queue) air leapfrog (regular runs, irregular cells, PML),
//       each head's explicit part in row chunks, each head's nonlinear
//       potential (an expensive potential runs here, hidden behind the air;
//       a TeamPotential batter, von Karman, is shared by a team of the
//       threads, one item per member at the head of the queue), and each
//       head's air planes (Psi^, J w^{n-1}, the plane RHS);
//   P2a (static) b_r = J rExpl + G y and b_g = J g per plane cell; g . rExpl
//       partials; the barrel wall's boundary layer (WallLoss.h) on every cell
//       but the head planes' (theirs runs in P1 before the plane RHS and in P3
//       after the coupling term);
//   P2b (static) z = (I_a + nu G)^-1 b for b_r and b_g, and the plane dots;
//   P3  every thread forms the same scalars (Sherman-Morrison, Newton), then
//       w^{n+1} by rows, the air planes' coupling term, and (thread 0) the
//       microphones, the stick, psi and the contact statistics.

#include "AirDomain.h"
#include "Collision.h"
#include "Coupling.h"
#include "Decimator.h"
#include "FdtdContracts.h"
#include "FdtdDrumConfig.h"
#include "HeadCalibration.h"
#include "HeadGrid.h"
#include "StretchingProbe.h"
#include "ThreadPool.h"
#include "WallLoss.h"

#include <array>
#include <atomic>
#include <memory>
#include <string>
#include <vector>

namespace taikor::candidates::fdtd
{

// Options the verification suite changes; the render uses the defaults.
struct FdtdOptions
{
    bool lossless = false;        // sigma0 = sigma1 = 0 on both heads
    bool rigidBox = false;        // reflecting outer box instead of the PML
    bool bergerOnRear = true;     // rear head nonlinearity (always on in renders)
    bool withBatterNonlinearity = true;
    bool withAir = true;          // false: heads in vacuo (no air planes)
    bool exteriorAir = true;      // false: only the barrel's air (rigid box hugging the drum)
    int threads = 0;              // 0: candidateThreadCount()
    // Reuse known inclusion radii (skips the calibration runs); <= 0: calibrate.
    double batterInclusion = 0.0, rearInclusion = 0.0;
    bool selfChecks = true;       // strike() runs the lossless-energy, PML and doublet checks
    // Tests only: a head spacing above the stability bound (the bound is a
    // minimum), for near-ideal membranes whose bound would be tiny.
    double minimumHeadSpacing = 0.0;
    // Energy, J, the heads' grids are built for at prepare() (their spacing
    // admits the largest SAV tension that energy allows, FdtdDrumConfig.h
    // headSpacing). 0: the linear bound; strike() rebuilds the heads of a
    // drum at rest whose stroke brings more (the harness path).
    double strokeEnergy = 0.0;
    // A batter potential whose stretching tension the energy does not bound
    // (StretchingProbe.h; von Karman) has its grid sized at the stroke by
    // sizing runs (FdtdEngineSizing.cpp). batterDesignTension > 0 instead
    // fixes the extra tension its grid is built for (tests, and the sizing
    // runs themselves); sizeBatterByStress = false keeps the Berger bound.
    bool sizeBatterByStress = true;
    double batterDesignTension = 0.0;
    // Sizing runs only: explicit head spacings (m), overriding the bounds.
    double batterSpacing = 0.0, rearSpacing = 0.0;
    // Threads sharing the batter's evaluate when its potential is a
    // TeamPotential (StretchingProbe.h); 0: automaticTeamSize (threads).
    // Performance only: a take is reproducible for a given size, and sizes
    // differ only by the in-plane solve's rounding (tests/integration
    // teamEvaluate).
    int potentialTeam = 0;
    // The Kirchhoff viscothermal boundary layer on the barrel's interior side
    // wall (WallLoss.h); off in a lossless drum. trackWallLoss (tests) sums
    // the energy it dissipates (dissipatedWall()).
    bool wallLoss = true;
    bool trackWallLoss = false;
};

// Members of the batter potential's team for a pool of `threads`: the team
// should finish its evaluate about when the other threads finish the air.
// On the nagado Don the von Karman evaluate is about 1.7-2 ms of thread time
// per step against about 8 ms of P1 air work, and its team scales at 1.4x for
// two members and 1.9x for three (barriers inside the multigrid V-cycle); on
// the shared 10-core machine with 8 threads the attack's 20 ms cost 2.46,
// 2.20, 2.31 and 2.49 ms per step with 1-4 members, so a quarter of the pool
// (tests/integration performance measures it again). Numerical/performance,
// never physical.
[[nodiscard]] inline int automaticTeamSize (int threads) noexcept
{
    return threads < 4 ? 1 : threads / 4;
}

struct HeadState
{
    HeadGrid grid;
    std::unique_ptr<NonlinearPotential> potential;
    InclusionCalibration calibration;
    HeadCoupling coupling;
    bool nonlinear = true;
    bool struck = false;          // carries the collision (batter)
    // Full-grid vectors.
    std::vector<double> w[3];     // w^n in w[n mod 3]
    std::vector<double> lap[2];   // L w^n in lap[n mod 2] (interior and first ring)
    std::vector<double> rExpl, gradient, g, footprint, xc;
    std::vector<std::vector<double>> scratch; // per-thread L w for chunk halos
    // Air-plane vectors (size planeSize).
    std::vector<int> cellAbove, cellBelow;
    bool hasAbove = true, hasBelow = true; // air on each side of the head
    std::vector<double> jwPrev, y, gy, br, bg, zr, zg, zc, jc;
    // With the barrel's wall loss: each plane cell's wall entry (-1: none;
    // WallLoss::sliceEntry), the factor 1 / (1 + beta) its implicit thermal
    // term puts on the head's flux into it, and the coupling's per-cell scale
    // S = mean of the factors over the air sides, so the head solves
    // (alpha I + mu I S J) (Coupling.h). Empty without wall loss.
    std::vector<int> wallAbove, wallBelow;
    std::vector<double> factorAbove, factorBelow, sliceScale;
    // Scalars of the current step, written by the potential item.
    double V = 0.0, gg = 0.0, gc = 0.0, gwPrev = 0.0;
    // The potential's StretchingProbe, if it has one, and the tension the SAV
    // step applies this step: N* psi^{n-1/2} / sqrt(2 V + eps) (StretchingProbe.h).
    const StretchingProbe* probe = nullptr;
    // The potential's TeamPotential (batter only) and the members sharing
    // each evaluate: P1 items {0, head, member} for member < teamSize.
    TeamPotential* team = nullptr;
    int teamSize = 1;
    double savTension = 0.0;
    double peakSavTension = 0.0, peakSavTensionTime = 0.0; // over the take so far (s after the strike)
    double windowSavTension = 0.0;                         // peak over the sizing window (FdtdEngineSizing.cpp)
    // SAV consistency over the take, for a potential whose stretching the
    // energy does not bound (von Karman; the Berger letter is not monitored).
    // psi^{n-1/2} is compared with sqrt(V^n + V^{n-1} + eps), V the take's own
    // evaluates (their energy error is second order in the solver residual;
    // tests/integration vonKarmanTracking checks them against a 1e-9 solve),
    // and the misstated force is expressed as the spurious tension
    // N*(w^n) (psi / sqrt(2 V) - 1) (VonKarmanPotential.h, "SAV consistency").
    struct SavConsistency
    {
        double ratio = 1.0, spurious = 0.0;        // at the last step
        double worstAttack = 0.0, worstAttackTime = 0.0; // largest |spurious| within sizingWindow (s after the strike)
        double worstTail = 0.0, worstTailTime = 0.0;     // largest |spurious| after it
        double worstOffset = 0.0;                        // largest |psi^2 / 2 - (V^n + V^{n-1}) / 2|, J
    };
    bool monitorSav = false;
    double stretchTension = 0.0; // N* of this step's evaluate (StretchingProbe)
    double vPrevious = 0.0;      // V^{n-1}, for the consistency monitor
    SavConsistency consistency;
    double psi[2] = { 0.0, 0.0 }; // psi^{n+1/2} in psi[n & 1]
    // Constants.
    double mu = 0.0, nu = 0.0, alpha = 1.0, beta = 0.0, gamma = 0.0, loadScale = 0.0;
    double cXc = 0.0, footprintCells = 0.0;
    int footprintBegin = 0, footprintEnd = 0; // rows holding footprint nodes
    int chunks = 1;
};

// Scalars summed over threads after P2a / P2b.
enum Partial : int
{
    PartialGR, PartialCR, PartialGJY, PartialCJY, PartialGZR, PartialGZG, PartialCZR, PartialCZG, PartialGZC,
    PartialCount
};

// Numerical window, s after the strike: 20 ms holds the Don's ~5 ms contact
// and any re-contact of the attack. The contact force history is kept over
// it, and the batter's sizing runs measure their peak stretching tension
// over it (FdtdEngineSizing.cpp); the take itself reports its peak over its
// whole length against the design, so a later excursion would show.
constexpr double sizingWindow = 0.02;

struct ContactStatistics
{
    double firstContact = -1.0, lastContact = -1.0; // s after the strike
    double peakForce = 0.0, impulse = 0.0, peakPenetration = 0.0;
    int episodes = 0;
    bool inContact = false;
    double exitSpeed = 0.0; // stick velocity (+ out) at the end of the last contact
    std::vector<double> force; // contact force per internal step over the first 20 ms (diagnostics)
};

template <typename Real>
class FdtdEngine final : public Engine
{
public:
    FdtdEngine (std::unique_ptr<NonlinearPotential> batterPotential, FdtdOptions options = {});
    ~FdtdEngine() override;

    std::string description() const override;
    bool prepare (const DrumSpec& drum, const MicSpec& mics, double sampleRate, std::string& error) override;
    void strike (const StrokeSpec& stroke) override;
    void process (float* left, float* right, int samples) override;
    std::string diagnostics() const override;

    // Verification access: internal steps, energies at the last half step.
    void advance (int steps);
    long long stepsTaken() const noexcept { return step; }
    double energyHeads() const;
    double energyAir() const;   // with the wall loss's boundary storage
    double dissipatedWall() const; // energy the wall loss has dissipated (trackWallLoss), J
    const WallLoss<Real>& wallLoss() const noexcept { return wall; }
    double energyStick() const;
    double energySav() const;
    double dissipatedHeads() const noexcept { return dissipated; }
    double internalRate() const noexcept { return 1.0 / k; }
    const ContactStatistics& contact() const noexcept { return stats; }
    const HeadState& batterState() const noexcept { return heads[0]; }
    const HeadState& rearState() const noexcept { return heads[1]; }
    const AirDomain<Real>& air() const noexcept { return domain; }
    double stickPosition() const noexcept { return stickY[step % 3]; }
    double lastForce() const noexcept { return force; }
    double micPressure (int capsule) const noexcept { return lastMic[capsule]; }
    std::unique_ptr<NonlinearPotential> releaseBatterPotential() { return std::move (heads[0].potential); }
    // Test hook: head `index` at rest in the given displacement (w^n = w^{n-1}).
    void setHeadAtRest (int index, const std::vector<double>& shape);
    // The extra tension (N/m) head `index`'s grid admits (FdtdDrumConfig.h
    // headSpacing): the Berger energy bound, or the batter's sized design.
    double designTension (int index) const { return tensionBound (index); }
    // Seconds the batter's sizing runs took at the stroke (0 if none ran),
    // and the peak N* over sizingWindow of the run it accepted (the take's
    // own first sizingWindow, bit for bit).
    double sizingTime() const noexcept { return sizingSeconds; }
    double sizingPeak() const noexcept { return sizingAcceptedPeak; }
    // Non-empty when a rebuild at the stroke failed (a potential's prepare()
    // at the new grid): process() then outputs silence, so the render tool
    // fails loudly, and diagnostics() says why.
    const std::string& failure() const noexcept { return failed; }
    // Profiling (TAIKOR_FDTD_PROFILE at prepare()): seconds per thread in each
    // P1 item kind (0 potential, 1 air planes, 2 head explicit, 3 air regular,
    // 4 air irregular, 5 PML), then P2a, P2b, P3 and barrier waits (6-9).
    const std::vector<std::array<double, 10>>& profileSeconds() const noexcept { return profile; }

private:
    void runSteps (int steps, float* left, float* right);
    void stepPhases (int thread, long long n);
    void phase1Item (int thread, int item, long long n);
    void phase2a (int thread, long long n);
    void phase2b (int thread, long long n);
    void phase3 (int thread, long long n);
    void recentre (int thread, long long n);
    void buildHeads (const DrumSpec& drum, std::string& error, bool batterOnly = false);
    void sizeBatter (const StrokeSpec& stroke);
    bool needsBatterSizing() const;
    void buildAir (const DrumSpec& drum);
    void buildHeadCoupling (const DrumSpec& drum);
    double tensionBound (int index) const;
    void precompute();
    void runSelfChecks (const DrumSpec& drum, const MicSpec& mics, double sampleRate);

    FdtdOptions options;
    FdtdNumerics numerics;
    DrumSpec drumSpec;
    MicSpec micSpec;
    StrokeSpec strokeSpec;
    double k = 0.0, outputRate = 48000.0;
    AirDomain<Real> domain;
    WallLoss<Real> wall;
    SquareRootFit wallFit;
    bool wallActive = false;
    std::vector<double> wallLost;           // per thread, with trackWallLoss
    HeadState heads[2];
    HertzContact hertz;
    double stickMass = 0.0;
    double stickY[3] = { 0.0, 0.0, 0.0 };
    double force = 0.0;
    bool struckYet = false, prepared = false;
    long long step = 0;
    std::unique_ptr<ThreadPool> pool;
    std::unique_ptr<Barrier> barrier;
    std::atomic<int> queue { 0 };
    int itemCount = 0;
    std::vector<std::array<int, 3>> items; // {kind, a, b}
    std::vector<int> pmlSplit;              // P2a: plane bounds per thread, balanced by PML cells
    double interiorShift = 0.0;             // pending shift of the interior potential (see phase3)
    std::vector<std::array<int, 2>> airColumns; // air columns over the heads (both heads share them)
    AirPlaneLayout airLayout;
    double headEnergy = 0.0;                // J the head spacings admit (FdtdOptions::strokeEnergy)
    double rebuildSeconds = 0.0;            // head rebuild at the stroke
    bool strokeBeyondBound = false;         // a stroke on a moving drum brought more than headEnergy
    bool outputPreroll = false;             // process() first runs the FIR's delay ahead (see process)
    // The batter's sizing (FdtdEngineSizing.cpp): the design tension its grid
    // is built for (0: the Berger bound), the calibrated inclusion radius the
    // accepted sizing run found for that grid, and the runs' record.
    double batterDesign = 0.0, batterKnownInclusion = 0.0;
    double sizingSeconds = 0.0, sizingAcceptedPeak = 0.0;
    std::string sizingReport;
    std::string failed;                     // see failure()
    void fail (const std::string& what, const std::string& error);
    int threadCount = 1;                    // the pool's size, fixed at prepare()
    int batterTeam = 1;                     // members asked of a TeamPotential batter (<= threadCount)
    std::vector<std::array<double, PartialCount * 2>> partials; // [thread][head * PartialCount + p]
    Decimator decimator[2];
    double lastMic[2] = { 0.0, 0.0 };
    std::vector<float> pendingLeft, pendingRight;
    int pendingCount = 0;
    ContactStatistics stats;
    double dissipated = 0.0;
    // Diagnostics. Profiling (env TAIKOR_FDTD_PROFILE=1): seconds per thread
    // in each P1 item kind (0-5), P2a, P2b, P3 and barrier waits (6-9).
    bool profiling = false;
    std::vector<std::array<double, 10>> profile;
    double prepareSeconds = 0.0, processSeconds = 0.0;
    long long processedSteps = 0;
    std::string checkReport;
};

} // namespace taikor::candidates::fdtd
