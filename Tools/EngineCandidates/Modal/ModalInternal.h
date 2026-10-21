#pragma once

// Shared state of the precompute (ModalBuild.cpp, ModalBlock.cpp, ModalMics.cpp).

#include "ModalBodyExterior.h"
#include "ModalBuild.h"

namespace taikor::candidates::modal
{

struct HeadLoss
{
    double sigma0 = 0.0, sigma1 = 0.0, arealDensity = 0.0;
};

struct DrumContext
{
    BuildOptions options;
    AirSpec air;
    HeadMaterial batterMaterial, rearMaterial;
    HeadLoss batterLoss, rearLoss;
    std::vector<std::vector<HeadMode>> batter, rear;          // per m
    std::vector<GradientForm> batterGradient, rearGradient;   // per m (empty forms beyond the list)
    std::vector<std::vector<CrossMode>> cross;                // per m
    std::unique_ptr<ExteriorTable> exterior;     // per-head baffled King loading
    std::unique_ptr<BodyExterior> body;          // the drum body's exterior (coupled band)
    double wall = 0.0;
    double omegaCoupled = 0.0, omegaAir = 0.0;
};

// Builds modes, gradient forms, cross-mode catalogue and the exterior tables
// (with the capsules as the body table's field points when mics is given).
bool prepareContext (const DrumSpec&, const BuildOptions&, DrumContext&, std::string& error, const MicSpec* mics = nullptr);

// The drum body's correction to a coupled mode's capsule transfers
// (CoupledMode::micBody), from the body table (ModalBodyExterior.h).
void computeBodyTransfers (const DrumContext&, CoupledMode&);

// Solves the coupled block of order m (head modes below omegaCoupled and air
// modes below omegaAir): frequencies by the self-consistent exterior loading,
// mass-normalised participations, energy shares and damping. The modes are
// appended to out.
void solveCoupledBlock (const DrumContext&, int m, std::vector<CoupledMode>& out);

// The batter modes of order m above omegaCoupled, each on its own with its self
// exterior loading and evanescent interior loading at its own frequency.
void solveHighModes (const DrumContext&, int m, std::vector<CoupledMode>& out);

// Microphone transfers of one mode (ModalMics.cpp).
void computeMicTransfers (const DrumContext&, const MicSpec&, CoupledMode&);

// Runs body(i) for i in [0, count) on the given number of threads, taking the
// items in order (so the largest items can be placed first).
void parallelFor (int count, int threads, const std::function<void (int)>& body);

} // namespace taikor::candidates::modal
