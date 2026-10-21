#pragma once

// Foppl-von Karman membrane stretching of a head with an immovable (tacked)
// edge, as a NonlinearPotential for the finite-difference drum's SAV step.
//
// Continuum model (Foppl 1907; von Karman 1910; plane stress):
//   V(w) = min_u 1/2 int N : eps dA,
//   eps  = sym grad u + 1/2 grad w (x) grad w,
//   N    = C eps, C = E h / (1 - nu^2) [[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu)/2]]
//          on (eps_xx, eps_yy, gamma_xy = 2 eps_xy),
//   u = 0 on r = a (the byo-uchi head is nailed round the shell; Thomas &
//   Bilbao, JSV 315, 2008, show the in-plane edge condition changes the
//   nonlinear behaviour qualitatively, so the traction-free Airy form with
//   F = d_n F = 0 would model a different drum).
//
// In-plane inertia is neglected: the in-plane field is the quasi-static
// minimiser U* at every step, as in the Foppl-von Karman / SAV schemes this
// letter follows and as the candidate brief specifies. That holds while the
// stretching source varies slowly against the head's in-plane modes. For
// this batter as a plane-stress disk with u = 0 at r = a (c_L = sqrt(C /
// sigma) = 1961 m/s, c_T = sqrt((1 - nu) C / (2 sigma)) = 1160 m/s;
// displacement potentials J_n(k_L r) cos n theta and J_n(k_T r) sin n theta
// give the edge determinant alpha beta J_n'(alpha) J_n'(beta) - n^2
// J_n(alpha) J_n(beta) = 0, alpha = omega a / c_L, beta = omega a / c_T)
// the lowest in-plane modes are 1566 Hz (n = 1), 1814 Hz (n = 0 torsional,
// J_1(beta) = 0) and 2439 Hz (n = 2). The Don's stretching rises and falls
// within the 0.3-4 ms contact and carries the footprint's ~9.4 kHz patch
// oscillation, so part of it lies above those modes, where a real in-plane
// field would lag U* (and ring near the modes) rather than follow it. That
// part is small: over the first 30 ms of the production take the mean
// stress under the stick holds 0.49 % of its energy above 1566 Hz (0.47 %
// above 5 kHz, i.e. nearly all of it the chatter), the stress 0.1 m outboard
// 3e-5 and V 3.6e-4 (a probe of the drum recording (N_xx + N_yy) / 2 and V
// every step). The full system with in-plane inertia (S. Bilbao, O. Thomas,
// C. Touze & M. Ducceschi, "Conservative numerical methods for the full von
// Karman plate equations", Numer. Methods PDE 31, 2015, 1948-1970) is the
// upgrade path.
//
// The pre-tension T is the drum's linear part and is excluded here. Its cross
// term with the stretching strain is T int tr(eps) dA = T int div u dA +
// T/2 int |grad w|^2 dA; the first integral is T times the flux of u through
// the rim, zero for an immovable edge, and the second is exactly the linear
// membrane energy the drum already carries. In the discrete form the first
// term telescopes to the fixed edge values (zero) along every grid row and
// column, so the split is exact on the grid too.
//
// Discretisation (see InPlaneSolver.h for the staggered layout):
//   slopes a = D+x w on x-edges, b = D+y w on y-edges;
//   q_xx(node) = (a_e^2 + a_w^2)/4, q_yy(node) = (b_n^2 + b_s^2)/4,
//   q_xy(corner) = s_x s_y with s_x, s_y the averages of the two x- and two
//   y-edge slopes round the corner;
//   eps = D U / h + q, N = C eps, V_h = h^2/2 [sum_nodes N_xx eps_xx +
//   N_yy eps_yy + sum_corners N_xy gamma].
// U* minimises V_h over the free in-plane dofs (K U* = -h D^T C q, K SPD),
// and by the envelope theorem dV_h/dw = (dq/dw)^T h^2 N*: the returned
// gradient is the exact derivative of the returned discrete energy whenever
// U = U*, and its error is linear in the solver residual otherwise. The SAV
// scheme stays exactly energy-conserving for any gradient (Bilbao, Ducceschi
// & Zama, J. Comput. Phys. 472, 2023, remark after eq. 43); the residual
// tolerance is for the physics and tests/vonkarman shows the output does not
// move when it is tightened.
//
// Solver tolerance. The in-plane solve stops at ||b - K U|| <= tol ||b||
// with tol = 3e-3 (InPlaneSolver::Settings), from a warm start extrapolated
// linearly from the previous two solutions. On the test stroke (48 kHz grid
// stepped at 96 kHz, 50 ms) that keeps the head's volume velocity within
// about 1e-5 of a 1e-10 solve at about 2 PCG iterations per call
// (tests/vonkarman "tolerance"). An optional rule stops at
// tol (1 + T / s_q) instead (T the pre-tension, s_q the peak source stress),
// bounding the stress error by tol (T + s_q) rather than tol s_q; it halves
// the cost but its absolute stress error, about tol T, detunes the quiet
// ring-down by up to ~tol/2 in frequency, so it is off.
//
// Calling it from the drum. One instance per head; evaluate() is not
// re-entrant (it holds the warm-start history and scratch), distinct
// instances are independent, and it starts no threads, so it may run on its
// own thread concurrently with the drum's air update. Call it once per time
// step with the current w, as the SAV step expects; the warm start assumes
// consecutive calls come from consecutive steps (any other w just costs more
// iterations). prepare() takes about 0.1 s (0.10 s on the 48 kHz grid,
// 0.12 s on the 96 kHz one).
//
// Time step. The explicit SAV step treats the stretching force explicitly
// and stays stable whatever its size, but it is accurate only while the
// step resolves the local stretching tension. On the test stroke (in vacuo,
// with or without the hide loss) the mean stretching stress under the stick
// reaches 1.4-1.5 T when resolved; on a grid at the linear stability limit
// of its own step (48 kHz with h = 6.9 mm, or 96 kHz with h = 4.8 mm) the
// step under-resolves it, the SAV variable psi decouples from sqrt(2 V)
// about 1 ms into the contact and the stress runs to 3-4 T. Stepping the
// 48 kHz grid at 96 kHz tracks psi to 1.2 %, the 96 kHz grid at 192 kHz to
// 0.7 %, and the two agree on the peak energy to 4 % (tests/vonkarman
// "resolution"). What the step must resolve is the stretching tension
// itself: the explicit SAV force acts on the grid's top modes like an extra
// tension N* (peakStretchingTension below), so the grid is stable only at
// the spacing of the stability bound at T + N*. The drum therefore sizes
// this letter's batter grid by measuring N* on runs of the stroke
// (Fdtd/FdtdEngineSizing.cpp: the Don peaks at N* = 1.52 T, h 4.92 mm at
// 96 kHz), and the take then tracks psi to sqrt(2 V) within 0.5 % over the
// 60 ms attack (tests/integration vonKarmanTracking); after that the ratio
// grows as V decays, see the next paragraph.
//
// SAV consistency over the take. psi is advanced from the gradient
// (psi^{n+1/2} = psi^{n-1/2} + g . (w^{n+1} - w^{n-1}) / 2), never
// recomputed from V, so it keeps whatever offset from sqrt(2 V) the attack
// leaves (the explicit step's consistency error at the strongest
// stretching, and the solver residual in g): energy is conserved exactly
// with psi, not with V. On the Don the gap psi - sqrt(V^n + V^{n-1})
// settles at about 1.1e-4 sqrt(J) by 0.5 s and stays (a verification run of
// the production drum to 1 s against a 1e-10 solve), so psi / sqrt(2 V) -
// 1 grows without bound as V decays (3e-2 at 0.5 s, 0.4 at 1 s, 5e4 at
// 4.4 s, where 2 V nears eps), while the force it misstates, the spurious
// tension N* (psi / sqrt(2 V) - 1), does not: for a quartic V, N* /
// sqrt(2 V) depends only on the shape of w, so the spurious tension follows
// the shape of the residual motion rather than time. The drum measures it
// at every step against the take's own V (Fdtd/FdtdEngineReport.cpp, "SAV
// consistency"; tests/integration vonKarmanTracking finds that measure
// within 4e-7 T of a 1e-9 solve's). On the Don it peaks at 99 N/m
// (1.6e-2 T) at the contact's onset, which the explicit step lags by about
// a step, and after the attack's 20 ms stays below 1.92 N/m (3.1e-4 T, a
// brief excursion at 4.03 s; 0.50-0.69 N/m at every 0.5 s of the 4.4 s
// take), which moves a partial by under 0.3 cent (tests/integration
// vonKarmanWholeTake holds it below 1e-3 T). Late in the take that
// spurious tension is larger than the true, decayed stretching, but both
// are far below T. Power-balanced drift
// regulation (T. Risse, T. Helie & S. Bilbao, DAFx-25) would remove the
// offset, but the SAV step is shared with the Berger letter, and the two
// letters must differ only in the potential, so it is not used.
//
// Summing q_xx + q_yy over the nodes gives (1/2) sum over all edges of the
// squared slopes, the same discrete int |grad w|^2 the Berger potential and
// the 5-point membrane energy use, so the two nonlinear letters share one
// discrete geometry.

#include "Fdtd/StretchingProbe.h"
#include "FdtdContracts.h"
#include "InPlaneSolver.h"

#include <array>
#include <memory>
#include <vector>

namespace taikor::candidates
{

class VonKarmanPotential final : public NonlinearPotential, public fdtd::StretchingProbe, public fdtd::TeamPotential
{
public:
    struct Stats
    {
        long calls = 0;
        long iterations = 0;        // PCG iterations over all calls
        int maxIterations = 0;      // worst call
        int lastIterations = 0;
        double lastResidual = 0.0;  // relative residual of the last solve
        double maxResidual = 0.0;
        double initialResidualSum = 0.0; // relative residual of the warm-start guess, summed
        double initialResidualMax = 0.0;
        double toleranceSum = 0.0;       // effective relative tolerance, summed
    };

    // The name the drum prints in key.md's description and diagnostics: the
    // model and its references (CandidateEngine.h asks for them).
    [[nodiscard]] std::string name() const override
    {
        return "Foppl-von Karman stretching (A. Foppl 1907; T. von Karman 1910), quasi-static in-plane field with an "
               "immovable edge (O. Thomas & S. Bilbao, JSV 315, 2008), as in S. Bilbao, C. Webb, Z. Wang & "
               "M. Ducceschi, \"Real-time gong synthesis\", DAFx-23";
    }
    bool prepare (const MembraneGrid& grid, const HeadSpec& head,
                  double headRadius, std::string& error) override;
    double evaluate (const double* w, double* gradient) override;

    // One line for key.md: solver work and accuracy over the calls so far
    // (the drum reads it through StretchingProbe::report()).
    [[nodiscard]] std::string diagnostics() const;

    // StretchingProbe (Fdtd/StretchingProbe.h). The bound on the
    // stress-stiffness of the last evaluate(), from the quadratic form
    //   d^T G d = 2 h^2 [sum_nodes N_xx q_xx(d) + N_yy q_yy(d)
    //                    + sum_corners N_xy q_xy(d)]
    // (q(d) the discrete strain terms of d, applyStressOperator): with
    // q_xx = (a_e^2 + a_w^2)/4, q_yy likewise and |q_xy| = |s_x s_y| <=
    // (s_x^2 + s_y^2)/2, s^2 <= (a_1^2 + a_2^2)/2, every edge slope collects
    // the mean normal stress along it over its two nodes plus the mean |N_xy|
    // over its two corners, so
    //   N* = max over edges of [ (N_nn(end 1) + N_nn(end 2)) / 2
    //                            + (|N_xy(corner 1)| + |N_xy(corner 2)|) / 2 ]
    // bounds G by N* K_1 (d^T K_1 d = h^2 sum over edges of the squared
    // slopes). On the checkerboard mode, the one the stability bound is set
    // by, the shear terms vanish (its corner-averaged slopes are zero) and
    // G acts as the mean normal stress, so N* is conservative by at most the
    // shear.
    // N* bounds G only, not the whole Hessian of V. The rest, the in-plane
    // material stiffness dq[d]^T C_rel dq[d] (through dN/dw, C_rel the
    // in-plane-relaxed stiffness), is positive and does not vanish on the
    // checkerboard: there only dq_xy[d] vanishes, while dq_xx[d] = (a_e da_e
    // + a_w da_w) / 2 = -(-1)^(i+j) w_xx delta for d = delta (-1)^(i+j)
    // (w_xx the second difference of w over h^2), nonzero wherever w is
    // curved, so per node the term is of order C h^2 w_xx^2 delta^2 against
    // G's 4 (N_xx + N_yy) delta^2. It is small on the Don: at the take's
    // peak-N* state (1.30 ms; Lanczos, 250 iterations, full
    // reorthogonalisation) the largest eigenvalue of the stiffness the SAV
    // step applies, T K_1 + D B + H, is 0.942 of that of (T + design) K_1 +
    // D B the grid is sized for, against 0.940 with G in place of H and
    // 0.910 for the linear part alone; the material term adds 0.2 % and the
    // sized grid keeps a 5.8 % margin (FdtdEngineSizing.cpp). Not bounded by
    // the energy: the stress concentrates under the stick (tests "ratios":
    // 67-74 times Berger's uniform tension for a static dent), and the energy
    // bounds only its mean square.
    [[nodiscard]] double peakStretchingTension() const override;
    [[nodiscard]] bool boundedByEnergy() const override { return false; }
    [[nodiscard]] std::string report() const override { return diagnostics(); }

    // Test and diagnostic hooks; the drum uses only the base interface.
    [[nodiscard]] vk::InPlaneSolver::Settings& solverSettings() noexcept { return solver_.settings(); }
    [[nodiscard]] const vk::InPlaneSolver& solver() const noexcept { return solver_; }
    // The tension-relative stopping rule above (off by default; see there).
    void setTensionRelativeStopping (bool on) noexcept { tensionRelative_ = on; }
    [[nodiscard]] bool tensionRelativeStopping() const noexcept { return tensionRelative_; }
    // 0: start every solve from zero; 1: from the previous solution;
    // 2: from the linear extrapolation of the previous two (default).
    void setWarmStartOrder (int order) noexcept { warmOrder_ = order; }
    void resetWarmStart() noexcept { history_ = 0; }
    [[nodiscard]] const Stats& stats() const noexcept { return stats_; }
    void resetStats() noexcept { stats_ = {}; }
    // The in-plane field of the last evaluate (solver level-vector layout).
    [[nodiscard]] const std::vector<double>& inPlaneField() const noexcept { return u0_; }
    // Energy of w at the given in-plane field (level-vector layout), without
    // minimising: V_h(w, U) >= V_h(w) with equality at U = U*.
    [[nodiscard]] double energyAt (const double* w, const double* inPlane);
    // Stress resultants N (N/m) of the last evaluate: normal components at
    // the nodes, shear at the corners north-east of them (grid.size() each).
    [[nodiscard]] const double* stressXX() const noexcept { return nxx_.data() + pad_; }
    [[nodiscard]] const double* stressYY() const noexcept { return nyy_.data() + pad_; }
    [[nodiscard]] const double* stressXY() const noexcept { return nxy_.data() + pad_; }
    // out = h^2 (dq/dw [d])^T N with N frozen at the last evaluate: the
    // discrete -div(N grad d) h^2, the "stress-stiffness" operator. Linear
    // and symmetric in d; G(N(w)) w = dV/dw.
    void applyStressOperator (const double* d, double* out);

    // TeamPotential (Fdtd/StretchingProbe.h): the evaluate by a team of
    // threads (Team.h). Every member calls evaluateMember with the same w and
    // gradient and gets the same V; the gradient is complete when any member
    // returns. A team of one is evaluate(), bit for bit; a larger team rounds
    // its sums in a different order (reproducibly for a given size).
    int setTeam (int size) override;
    double evaluateMember (const double* w, double* gradient, int member) override;
    [[nodiscard]] int teamSize() const noexcept { return teamSize_; }

private:
    double evaluateTeam (const double* w, double* gradient, vk::Team& team);
    void computeStrainSources (const double* w);
    double energyAndGradient (const double* inPlane, double* gradient);
    // Slopes and corner-averaged slopes of a node field (padded copy in wp).
    void slopes (double* wp, double* ax, double* by, double* sx, double* sy) const;
    void slopesRows (const double* wp, double* ax, double* by, int k0, int k1) const;
    void cornerSlopesRows (const double* ax, const double* by, double* sx, double* sy, int k0, int k1) const;
    // computeStrainSources in stages over node ranges (see the .cpp).
    void sourcesMask (const double* w, int k0, int k1);
    double sourcesStrain (int k0, int k1); // returns the range's peak source trace
    void sourcesRhs (int k0, int k1);
    double energyRows (const double* inPlane, int k0, int k1); // sum of N : eps over the range
    // gradient = h^2 (dq/dw)^T N for the stored N and the given slopes.
    void stressDivergence (const double* ax, const double* by, const double* sx,
                           const double* sy, double* gradient);
    void divergenceFlux (const double* ax, const double* by, const double* sx, const double* sy, int k0, int k1);
    void divergenceGradient (double* gradient, int k0, int k1);

    MembraneGrid grid_;
    vk::PlaneStress c_;
    vk::InPlaneSolver solver_;
    int n_ = 0, N_ = 0, pad_ = 0;
    double h_ = 0.0;

    // Node-centred work arrays, each N_ + 2 pad_ long, index pad_ + k.
    std::vector<double> w_, ax_, by_, sx_, sy_, qxx_, qyy_, qxy_;
    std::vector<double> nxx_, nyy_, nxy_, fx_, fy_;
    std::vector<double> rhs_, u0_, u1_;
    std::vector<double> dw_, dax_, dby_, dsx_, dsy_; // applyStressOperator scratch
    int history_ = 0;
    int warmOrder_ = 2;
    bool tensionRelative_ = false;
    double tension_ = 0.0;         // HeadSpec::tension, N/m
    double sourceStress_ = 0.0;    // s_q of the last computeStrainSources
    mutable double peakTension_ = -1.0; // N* of the last evaluate (< 0: not computed yet)
    Stats stats_;
    int teamSize_ = 1;
    std::unique_ptr<fdtd::Barrier> teamBarrier_;
    vk::TeamShared teamShared_;
    std::array<vk::Team, vk::maxTeam> teams_ {}; // each member's Team
};

} // namespace taikor::candidates
