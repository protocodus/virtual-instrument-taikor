#pragma once

// Precompute of the coupled-mode drum: per angular order m, one real symmetric
// generalised eigenproblem over the batter's and the rear head's in-vacuo
// modes (ModalMath.h) and the barrel's rigid-cylinder air modes (ModalCavity.h),
// with the exterior added mass of the whole drum body (ModalBodyExterior.h,
// both heads in one field; per-head baffled King loading above its table,
// ModalExterior.h) iterated to self-consistency at every coupled mode's own
// frequency (Christian, Davis,
// Tubis, Anderson, Mills & Rossing, "Effects of air loading on timpani
// membrane vibrations", JASA 76(5):1336, 1984; Tiwari & Gupta, "Effects of air
// loading on the acoustics of an Indian musical drum", JASA 141(4):2611, 2017),
// followed by first-order damping projection and the microphone transfer of
// every coupled mode at its complex pole. See ModalBuild.cpp for the method.

#include "CandidateEngine.h"
#include "ModalCavity.h"
#include "ModalExterior.h"
#include "ModalMath.h"

#include <array>
#include <complex>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace taikor::candidates::modal
{

// Numerical parameters. Each is justified by a convergence test in
// tests/modal (named beside it); none carries physics.
struct BuildOptions
{
    double headMaxHz = 21600.0;    // batter modes resolved to 0.45 fs (the spec's band)
    // Head modes solved in the coupled blocks; build_coupled_band: the capsule
    // signal below 3 kHz is unchanged (0.1 dB) by 4.5 / 6 kHz, the barrel's
    // ring above is not (a stated limitation, ModalDrum.cpp).
    double coupledMaxHz = 3000.0;
    double airMaxHz = 4500.0;      // explicit barrel air modes (cavity_convergence)
    int gridPerOctave = 12;        // exterior-loading frequency grid (build_grid_convergence)
    double gradientTolerance = 1.0e-10; // Berger/sigma1 form truncation (heads_gradient_form)
    double crossGammaFactor = 3.0; // cross-modes to crossGammaFactor * b1max + 40/a (cavity_convergence)
    bool exterior = true;          // off only in tests that isolate the interior
    bool interior = true;
    // The drum body's exterior (ModalBodyExterior.h) over the coupled band; off: per-head baffles.
    bool bodyExterior = true;
    double bodyMaxHz = 1500.0;     // top of the body table (modes and frequencies); 0: 1.3 coupledMaxHz (body_exterior_convergence)
    double bodyPanelsPerWavelength = 10.0; // (body_exterior_convergence: 10 vs 16)
    int bodyChiefPoints = 40;      // (bem_drum_body)
    double dampingBandwidth = 8.0; // couple damping of modes closer than this x (alpha_k + alpha_l) (build_damping_coupling)
    double dampingFloor = 1.0e-3;  // ... when |A_kl| >= this x sqrt(alpha_k alpha_l)
    int threads = 8;
};

struct CoupledMode
{
    int m = 0;
    bool high = false;             // above coupledMaxHz: a lone loaded batter mode
    double omega = 0.0, alpha = 0.0;
    // Mass-normalised participation on the in-vacuo modes of order m
    // (indices into ModalModel::batter[m] / rear[m]); q_i = sum_k coef x_k.
    std::vector<int> batterIndex, rearIndex;
    std::vector<double> batterCoef, rearCoef;
    double shareBatter = 0.0, shareRear = 0.0, shareAir = 0.0; // kinetic-energy shares
    double alphaHide = 0.0, alphaRadiation = 0.0, alphaWall = 0.0;
    // Off-diagonal decay-rate couplings (other mode, A_kl in 1/s) to close modes
    // of the same block; indices are block-local in solveCoupledBlock and global
    // (into ModalModel::modes) after buildModel.
    std::vector<std::pair<int, double>> damping;
    // Microphone transfers at the complex pole s = -alpha + j omega, without the
    // batter's angular factor cos(m (theta_mic - theta_strike)), which is applied
    // at strike time: p_mic = Re{ (micBatter cos(..) + micRear) A(t) }, A = the
    // analytic modal displacement.
    std::array<std::complex<double>, 2> micBatter {}, micRear {};
    // What the drum body adds to those at the capsules (diffraction at the rim and
    // barrel, the rear head's field round the body), for the modes of the body
    // table's band, also without the angular factor (ModalMics.cpp). Reported in
    // the diagnostics, not heard: see ModalMics.cpp.
    std::array<std::complex<double>, 2> micBody {};
};

struct ModalModel
{
    double radius = 0.0;
    HeadMaterial batterMaterial, rearMaterial;
    std::vector<std::vector<HeadMode>> batter, rear;   // per m
    std::vector<GradientForm> batterGradient, rearGradient;
    std::vector<CoupledMode> modes;
    // Summary for diagnostics.
    int airModes = 0, coupledBlocks = 0, highModes = 0, batterCount = 0, rearCount = 0, dampingLinks = 0;
    int bodyPanels = 0, bodyGridPoints = 0;
    double bodyMaxHz = 0.0, bodySeconds = 0.0;
    double bodyAsymmetry = 0.0;    // largest |Z_ic - Z_ci| / sqrt|Z_ii Z_cc| before symmetrising (ModalBodyExterior.cpp)
    double seconds = 0.0;
    std::string log;
};

[[nodiscard]] std::unique_ptr<ModalModel> buildModel (const DrumSpec&, const MicSpec&, const BuildOptions&,
                                                      std::string& error);

// One coupled block without damping or microphones, for the tests:
// frequencies (rad/s) and energy shares, ascending.
struct BlockSummary
{
    std::vector<double> omega, alpha, shareBatter, shareRear, shareAir;
};
[[nodiscard]] BlockSummary solveBlockForTest (const DrumSpec&, const BuildOptions&, int m);

// The baffled Rayleigh integral of a radial profile u(r) cos(m theta) at complex
// frequency s, seen from a capsule at radius rhoMic (at azimuth 0 of the profile)
// and height zMic: int u(r) cos(m theta) e^{-s R/c} / (2 pi R) dS (ModalMics.cpp).
// maxWavenumber bounds the profile's radial wavenumber; refine multiplies both
// node counts (convergence tests).
[[nodiscard]] std::complex<double> rayleighIntegral (int m, double radius, const std::function<double (double)>& profile,
                                                     double maxWavenumber, std::complex<double> s, double soundSpeed,
                                                     double rhoMic, double zMic, int refine);
// The quadrature nodes of rayleighIntegral for |s|/c = kAbs: visit(weight, R) for
// every node, so that the integral is sum weight e^{-s R/c} (the same nodes give
// the delay density of the transfer, ModalFeedthrough.cpp).
void rayleighNodes (int m, double radius, const std::function<double (double)>& profile, double maxWavenumber,
                    double kAbs, double rhoMic, double zMic, int refine,
                    const std::function<void (double weight, double distance)>& visit);
// A coupled mode's batter profile u(r) = sum x_i phi_i(r) (the returned function
// refers to mode and batter, which must outlive it) and its largest b1.
[[nodiscard]] std::function<double (double)> batterProfile (const CoupledMode& mode,
                                                            const std::vector<std::vector<HeadMode>>& batter,
                                                            double radius, double& maxWavenumber);
// A coupled mode's rear volume displacement sum x_j V_j (m = 0 only) and the
// straight distance from the rear head's centre to a capsule (ModalMics.cpp).
[[nodiscard]] double rearVolume (const CoupledMode& mode, const std::vector<std::vector<HeadMode>>& rear);
[[nodiscard]] double rearDistance (double length, const std::array<double, 3>& capsule);

} // namespace taikor::candidates::modal
