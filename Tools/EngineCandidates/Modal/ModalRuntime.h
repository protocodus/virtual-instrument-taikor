#pragma once

// Time stepping of the coupled-mode drum with the bachi.
//
// MODES. Each coupled mode (mass-normalised coordinate y, frequency w, decay
// alpha) is stepped by the two-point Hamiltonian scheme of van Walstijn &
// Bridges, "Simulation of distributed contact in string instruments: a modal
// expansion approach", EUSIPCO 2016 (their eqs. 18-21): the trapezoidal rule in
// first-order form, y^{n+1} - y^n = q^{n+1} + q^n, q^{n+1} - q^n =
// -a mu y - b delta y + xi F, with scaled momentum q = (dt/2) p, xi = dt^2/2 and
// the "pre-warped" coefficients
//   a* = (1 - 2 R O + R^2)/(1 + 2 R O + R^2),  b* = 2 (1 - R^2)/(1 + 2 R O + R^2),
//   R = exp(-alpha dt), O = cos(w dt),
// which reproduce each mode's (w, alpha) exactly. The numerical energy
// (q^2 + a y^2)/xi falls by (b/xi) s^2 per step and is otherwise conserved.
//
// CONTACT. One scalar gap eta = wbar - z between the footprint-averaged batter
// displacement wbar = sum c_k y_k and the stick tip z, with the one-sided Hertz
// potential Phi(eta) = K/(alpha_H + 1) [eta]_+^{alpha_H + 1}. The force is its
// discrete gradient over the step, F = (Phi(eta^{n+1}) - Phi(eta^n)) /
// (eta^{n+1} - eta^n) (Chatziioannou & van Walstijn, "Energy conserving schemes
// for the simulation of musical instrument contact dynamics", JSV 339:262,
// 2015; van Walstijn & Bridges 2016 eq. 17), so stick + bank + Phi conserve
// energy exactly in the lossless case. Because every mode's update is linear in
// F, eta^{n+1} - eta^n = U - W F with scalars U, W > 0, and the implicit step
// reduces to one monotone scalar equation solved by safeguarded Newton. The
// contact stays armed for the whole take: re-contact is automatic.
//
// BERGER. Each head's tension rises by dT = E h/(2 (1 - nu^2) pi a^2) S with
// S = int |grad w|^2 dS, a quadratic form in the modal coordinates, so the
// Berger energy beta S^2 (beta = E h/(8 (1 - nu^2) pi a^2)) is the square of
// psi = sqrt(2 beta) S. It is stepped with the explicit scalar-auxiliary-
// variable method (Bilbao, Ducceschi & Zama, "Explicit exactly energy-
// conserving methods for Hamiltonian systems", J. Comput. Phys. 472:111697,
// 2023): psi^{n+1} = psi^n + g . (y^{n+1} - y^n) with g = grad psi at the
// predicted midpoint y + q, force -(psi^{n+1} + psi^n)/2 g on the modes. The
// modified energy (modes + psi^2/2 + stick + Phi) is conserved exactly; the two
// heads' psi and the contact force form a 2x2 linear system and one scalar
// equation per step (Sherman-Morrison on the diagonal modal update).
//
// OVERSAMPLING. Samples in which contact is possible are stepped Nos times at
// dt/Nos (the exact modal coefficients are simply recomputed for the substep).
// The state is carried between the two rates through each mode's analytic
// amplitude A = y + j (y mu_r - q)/mu_i, mu = (zeta - 1)/(zeta + 1), zeta =
// exp((-alpha + j w) dt), which is the rate-independent (y, dy/dt) of the mode.
//
// OUTPUT. p_mic = sum_k Re{H_k A_k} - sum_l F_l W_l: the microphone transfer
// at the complex pole (ModalMics.cpp) applied to each mode's analytic
// amplitude, exact for the free decay, plus the contact force's own share of
// the head acceleration and the retardation of the forced motion across the
// head, as one kernel per capsule on the history of the substep forces F_l
// (ModalFeedthrough.cpp). The sum is causal and its time integral over a
// stroke vanishes (tests/modal output_*).

#include "CandidateEngine.h"

#include <array>
#include <complex>
#include <string>
#include <vector>

namespace taikor::candidates::modal
{

struct ModalModel;

struct BankMode
{
    double omega = 0.0, alpha = 0.0;
    double contact = 0.0;                          // c_k: footprint mean of the batter shape
    std::array<std::complex<double>, 2> mic {};    // H_k per capsule, Pa per metre of A
};

// One head's Berger quadratic form S = sum_groups [sum_i d_i q_i^2 - sum_r (l_r . q)^2]
// with q_i = sum_(k, coef) in rows[i] of coef y_k.
struct BergerGroup
{
    std::vector<double> diagonal;
    std::vector<std::vector<double>> lowRank;
    std::vector<std::vector<std::pair<int, double>>> rows;
};
struct BergerHead
{
    double beta = 0.0;                             // E h / (8 (1 - nu^2) pi a^2)
    double tensionPerS = 0.0;                      // E h / (2 (1 - nu^2) pi a^2), for diagnostics
    std::vector<BergerGroup> groups;
};

// Modes whose damping couples them (ModalBuild.h CoupledMode::damping): the
// decay-rate matrix A (p x p, row-major, diagonal = each mode's alpha) is
// stepped densely, which keeps the non-proportional damping of near-degenerate
// modes (a hide-damped head mode against a lightly damped barrel air mode does
// not hybridise when their decay rates differ by more than their coupling).
struct DampingCluster
{
    std::vector<int> members;
    std::vector<double> decay;
};

struct Bank
{
    std::vector<BankMode> modes;
    std::vector<DampingCluster> clusters;
    std::array<BergerHead, 2> berger;              // batter, rear
    double tipMass = 0.0, hertzStiffness = 0.0, hertzExponent = 1.5, impactSpeed = 0.0;
    // Contact-force term of the capsule pressure (ModalFeedthrough.cpp): per
    // capsule, the kernel K integrated over lag bins [l d, (l + 1) d) of width
    // d = feedthroughBin (the contact substep), in Pa per N; p -= sum_l F_l W_l.
    std::array<std::vector<double>, 2> feedthrough;
    double feedthroughBin = 0.0;
    // Build statistics.
    int prunedModes = 0, bergerEntries = 0, bergerDenseEntries = 0, clusteredModes = 0, clippedClusters = 0;
    double feedthroughSeconds = 0.0;
};

struct RuntimeOptions
{
    int oversample = 16;            // contact substeps per sample (contact_impulse_convergence)
    double pruneResidue = 1.0e-7;   // relative strike-to-capsule residue below which a mode is dropped
    double pruneContact = 1.0e-5;   // ... provided its footprint projection is also below this
    double bergerTolerance = 2.0e-4;// sparsification of the Berger head rows (modal_berger_sparsity)
    double bergerRankTolerance = 1.0e-6; // low-rank truncation of the Berger form (modal_berger_sparsity)
    bool berger = true;
    bool feedthrough = true;        // tests only: false drops the contact-force term of the pressure
    int forceRate = 0;              // tests only: 1 = never oversample, 2 = always oversample
    int threads = 8;
};

// Builds the bank for one stroke: angular factors of the capsules relative to
// the strike azimuth, the contact projection of every mode, residue pruning and
// the sparsified Berger rows; with sampleRate > 0 also the contact-force
// feedthrough kernel at the substep rate sampleRate * oversample.
[[nodiscard]] Bank makeBank (const ModalModel&, const DrumSpec&, const MicSpec&, const StrokeSpec&, const RuntimeOptions&,
                             double sampleRate = 0.0);

// The feedthrough kernel (ModalFeedthrough.cpp); original[i] is the model index
// of bank mode i.
void buildFeedthrough (Bank&, const ModalModel&, const DrumSpec&, const MicSpec&, const StrokeSpec&,
                       const std::vector<int>& original, double bin, int threads);

struct ContactReport
{
    int episodes = 0, baseRateContactSteps = 0, oversampledSamples = 0;
    double firstContact = -1.0, lastContact = -1.0; // s after the strike
    double peakForce = 0.0, impulse = 0.0, exitVelocity = 0.0, maxPenetration = 0.0;
    double energyStart = 0.0, energyEnd = 0.0, dissipated = 0.0; // over the first contact window
    bool energyWindowClosed = false;
    double peakTensionRise[2] = { 0.0, 0.0 };       // N/m, batter and rear
    double maxPsiDrift = 0.0;                       // max |psi - sqrt(2 beta) S(y)| / peak psi (every 64 samples)
};

class Stepper
{
public:
    Stepper (Bank bank, double sampleRate, const RuntimeOptions& options);
    void strike();
    void process (float* left, float* right, int samples);
    // Double-precision variant for tests (either pointer may be null).
    void processDouble (double* left, double* right, int samples);
    [[nodiscard]] const ContactReport& report() const noexcept { return contactReport; }
    [[nodiscard]] int modeCount() const noexcept { return n; }
    [[nodiscard]] double totalEnergy() const;       // numerical energy at the current rate
    [[nodiscard]] double stickVelocity() const;
    [[nodiscard]] double seconds() const noexcept { return renderSeconds; }
    // Tests: sets mode k's displacement and velocity (at the current rate).
    void setModeState (int k, double displacement, double velocity);
    // Tests: the modal displacements, and the time integral of the pressure of
    // the remaining free decay, -sum_k Re{H_k A_k / s_k}, per capsule.
    [[nodiscard]] const std::vector<double>& displacements() const noexcept { return y; }
    [[nodiscard]] std::array<double, 2> freeDecayIntegral() const;

private:
    struct Rate
    {
        double dt = 0.0, xi = 0.0, xiStick = 0.0;
        std::vector<double> a, d, bOverXi, muR, muI;
        std::vector<double> gc;                     // G c (G = (I + a + B)^-1, block-diagonal)
        std::vector<std::vector<double>> clusterG, clusterBOverXi; // per cluster, p x p
        double cec = 0.0;                           // xi c . G c
    };
    void applyG (const Rate&, const double* in, double* out) const;
    void setRate (Rate&, double dt);
    double step (const Rate&);                      // returns the contact force of the step
    void berger (const std::vector<double>& x, double& sB, double& sR);
    void convert (const Rate& from, const Rate& to);
    void pushForce (double force);
    void subtractFeedthrough (double& left, double& right) const;

    Bank data;
    RuntimeOptions options;
    int n = 0;
    Rate base, over;
    const Rate* current = nullptr;
    std::vector<double> y, q, c, xt, u, gB, gR, h1L, h2L, h1R, h2R, v1, v2, w, ggB, ggR, sv;
    std::vector<char> clustered;                    // mode belongs to a damping cluster
    std::vector<double> headQ, headDs;             // scratch for the Berger rows
    double z = 0.0, qs = 0.0, psi[2] = { 0.0, 0.0 };
    double eta = 0.0, etaRate = 0.0;
    bool struck = false, lastSampleContact = false;
    double time = 0.0;                              // seconds since strike
    ContactReport contactReport;
    bool inEnergyWindow = false, energyStarted = false;
    bool inContact = false;
    double stepDissipation = 0.0, stepS[2] = { 0.0, 0.0 };
    long long sampleCount = 0;
    double peakPsi[2] = { 0.0, 0.0 };
    double renderSeconds = 0.0;
    // Contact force per substep over the feedthrough kernel's span (ring buffer).
    bool useFeedthrough = false;
    std::vector<double> forceHistory;
    std::size_t historyHead = 0, forceAge = 0;
};

} // namespace taikor::candidates::modal
