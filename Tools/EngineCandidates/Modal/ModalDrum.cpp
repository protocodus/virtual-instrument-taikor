// The complete coupled-mode drum (engine id "modal"): both clamped stiff hides
// resolved mode by mode to 0.45 fs, coupled per angular order through the
// barrel's rigid-cylinder air modes with Kirchhoff's side-wall boundary-layer
// loss (ModalCavity.h), loaded by and radiating
// into the air around the whole drum body (ModalBodyExterior.h; baffled
// half-spaces above its table, ModalExterior.h), heard through near-field
// Rayleigh microphones at each mode's complex pole plus the contact force's
// feedthrough (ModalMics.cpp, ModalFeedthrough.cpp), struck by a re-contacting
// Hertz bachi against the whole bank with exact modal Berger tension on both
// heads (ModalRuntime.h). No statistical noise layer.

#include "CandidateEngine.h"
#include "ModalBuild.h"
#include "ModalRuntime.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <thread>

namespace taikor::candidates
{
namespace
{
using namespace modal;

int threadCount()
{
    if (const char* env = std::getenv ("TAIKOR_CANDIDATE_THREADS"))
    {
        const int n = std::atoi (env);
        if (n > 0)
            return n;
    }
    const unsigned hw = std::thread::hardware_concurrency();
    return hw > 0 ? static_cast<int> (std::min (hw, 8u)) : 8;
}

class ModalDrum final : public Engine
{
public:
    [[nodiscard]] std::string description() const override
    {
        return "Complete coupled-mode drum: clamped stiff-membrane modes of both heads (Leissa 1969; Suzuki, Yamaguchi "
               "& Izumi 2009), the batter resolved to 0.45 fs; coupled per angular order with the rear head and the "
               "barrel's rigid-cylinder air modes (Dowell, Gorman & Smith 1977), heads below 3 kHz and barrel air "
               "below 4.5 kHz. Exterior loading and radiation of the whole drum body, both heads in one field, by an "
               "axisymmetric boundary-element solution with CHIEF (Seybert, Soenarko, Rizzo & Shippy 1986; Soenarko "
               "1993; Schenck 1968) below 1.5 kHz, per-head baffled King loading above (King 1934), self-consistent at "
               "every mode's frequency (Christian et al. 1984; Tiwari & Gupta 2017). Barrel loss: Kirchhoff's "
               "viscothermal boundary layer on the rigid side wall (Kirchhoff 1868; Pierce, Acoustics, ch. 10; "
               "mu 1.81e-5 Pa s, Pr 0.71, gamma 1.4), the plane-wave law (1/a) sqrt(w nu/2) (1 + (gamma-1)/sqrt Pr) "
               "(12.4 / 14.3 / 16.0 / 17.5 dB/s at 3 / 4 / 5 / 6 kHz) applied to each coupled mode's whole wall field "
               "at its own frequency; its viscous part follows the velocity along the wall, so radial barrel modes "
               "keep only the thermal part (0.32 of the plane-wave rate). Non-proportional damping of close modes. "
               "Rayleigh near-field microphones at the complex poles plus the contact force's feedthrough and "
               "retardation (causal output). Re-contacting energy-conserving Hertz bachi on the whole bank (van "
               "Walstijn & Bridges 2016; Chatziioannou & van Walstijn 2015); exact modal Berger by SAV (Bilbao, "
               "Ducceschi & Zama 2023). No statistical noise layer. "
               "IN THIS LETTER'S ATTACK ONLY: the contact force's chatter near 9.6 kHz (the Hertz stiffness against "
               "the hide under the footprint; force ripple above 5 kHz up to 85 N of the 114 N peak) is carried "
               "because the batter is resolved to 21.6 kHz, and reaches the capsules through the force's "
               "feedthrough: the 8-11.5 kHz band peaks about 8 dB below the take's peak over the first 0.2-4.5 ms; "
               "the finite-difference letters' grids cannot carry it above about 5-9 kHz. "
               "NOT REPRESENTED: head-air coupling above 3 kHz and barrel air above 4.5 kHz, so the rigid barrel's "
               "high-frequency reverberation is short: at 3-4.5 kHz by 1.5 dB at 20-50 ms and 4-7 dB from 50 ms to "
               "0.6 s, at 4.5-8 kHz by 3 dB at 20-50 ms and 18-24 dB from 50 ms to 0.6 s (measured against a "
               "4.5 / 6 kHz band that is itself not converged, so lower bounds); the rear head above 3 kHz (nothing "
               "would drive it without head-air coupling); the drum body's exterior above 1.5 kHz (baffled heads "
               "there: against a body table to 2.6 kHz the decay at 1-2 kHz differs by at most 5.5 %); bulk "
               "atmospheric absorption (ISO 9613-1), which needs a humidity the specs do not carry; the Berger "
               "force's own pressure term (about "
               "2 % of the first 10 ms); and the body's diffraction at the capsules, which hear a baffled Rayleigh "
               "batter and a compact-monopole rear (the body's field, computed and listed in the diagnostics, would "
               "lower the (0,2) and (1,2) partials by 3-4 dB and the (0,1) pair by 1-2 dB).";
    }

    bool prepare (const DrumSpec& drumIn, const MicSpec& micsIn, double sampleRate, std::string& error) override
    {
        drum = drumIn;
        mics = micsIn;
        fs = sampleRate;
        build.headMaxHz = 0.45 * sampleRate;
        build.threads = threadCount();
        runtime.threads = build.threads;
        model = buildModel (drum, mics, build, error);
        stepper.reset();
        return model != nullptr;
    }

    void strike (const StrokeSpec& strokeIn) override
    {
        stroke = strokeIn;
        const auto t0 = std::chrono::steady_clock::now();
        Bank bank = makeBank (*model, drum, mics, stroke, runtime, fs);
        feedthroughSeconds = bank.feedthroughSeconds;
        feedthroughSpan = static_cast<double> (bank.feedthrough[0].size()) * bank.feedthroughBin;
        feedthroughStatic[0] = feedthroughStatic[1] = 0.0;
        for (int j = 0; j < 2; ++j)
            for (double w : bank.feedthrough[static_cast<std::size_t> (j)])
                feedthroughStatic[j] += w;
        bankModes = static_cast<int> (bank.modes.size());
        pruned = bank.prunedModes;
        bergerEntries = bank.bergerEntries;
        bergerDense = bank.bergerDenseEntries;
        clusters = static_cast<int> (bank.clusters.size());
        clippedClusters = bank.clippedClusters;
        clusteredModes = bank.clusteredModes;
        largestCluster = 0;
        for (const auto& cl : bank.clusters)
            largestCluster = std::max (largestCluster, static_cast<int> (cl.members.size()));
        stepper = std::make_unique<Stepper> (std::move (bank), fs, runtime);
        stepper->strike();
        bankSeconds = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
    }

    void process (float* left, float* right, int samples) override
    {
        if (! stepper)
        {
            std::fill (left, left + samples, 0.0f);
            std::fill (right, right + samples, 0.0f);
            return;
        }
        stepper->process (left, right, samples);
    }

    [[nodiscard]] std::string diagnostics() const override;

private:
    DrumSpec drum;
    MicSpec mics;
    StrokeSpec stroke;
    double fs = 48000.0;
    BuildOptions build;
    RuntimeOptions runtime;
    std::unique_ptr<ModalModel> model;
    std::unique_ptr<Stepper> stepper;
    int bankModes = 0, pruned = 0, bergerEntries = 0, bergerDense = 0, clusters = 0, clusteredModes = 0, largestCluster = 0;
    int clippedClusters = 0;
    double bankSeconds = 0.0, feedthroughSeconds = 0.0, feedthroughSpan = 0.0, feedthroughStatic[2] = { 0.0, 0.0 };
};

std::string ModalDrum::diagnostics() const
{
    std::ostringstream out;
    out.precision (4);
    if (! model)
        return "not prepared";
    out << "Modes: " << model->batterCount << " batter (cos partners, in vacuo to " << build.headMaxHz / 1000.0
        << " kHz), " << model->rearCount << " rear (to " << build.coupledMaxHz / 1000.0 << " kHz), "
        << model->airModes << " barrel air modes (to " << build.airMaxHz / 1000.0 << " kHz); "
        << model->coupledBlocks << " coupled angular blocks (heads below " << build.coupledMaxHz / 1000.0
        << " kHz), " << model->highModes << " lone loaded batter modes above; " << model->modes.size()
        << " coupled modes, " << bankModes << " kept for this stroke (" << pruned << " pruned by residue). "
        << "Non-proportional damping: " << clusters << " clusters of close modes (" << clusteredModes
        << " modes, largest " << largestCluster << ") stepped with their dense decay matrix (" << clippedClusters
        << " needed their couplings shrunk to stay passive).\n";
    if (model->bodyPanels > 0)
        out << "Exterior: drum-body boundary elements (" << model->bodyPanels << " panels, " << model->bodyGridPoints
            << " frequencies to " << model->bodyMaxHz << " Hz, " << model->bodySeconds
            << " s) for the head modes below that frequency, both heads in one field; per-head baffled King loading "
               "for the rest.\n";
    {
        // The barrel's side-wall loss (ModalCavity.h) and what it does to the
        // air-dominated coupled modes, band by band.
        const AirSpec air { drum.radius, drum.depth, drum.airDensity, drum.soundSpeed };
        out << "Barrel side wall: Kirchhoff viscothermal boundary layer (mu 1.81e-5 Pa s, Pr 0.71, gamma 1.4); "
               "plane-wave law ";
        for (double f : { 3000.0, 4000.0, 5000.0, 6000.0 })
            out << (f > 3000.0 ? " / " : "") << 8.685889638 * kirchhoffPlaneWaveRate (air, 2.0 * pi * f);
        out << " dB/s at 3 / 4 / 5 / 6 kHz, applied to every coupled mode's whole wall field at its own frequency "
               "(radial barrel modes keep only the thermal part, 0.32 of that law). Air-dominated coupled modes "
               "(air share > 0.5): ";
        const double edges[] = { 1000.0, 2000.0, 3000.0, build.airMaxHz };
        for (int b = 0; b < 3; ++b)
        {
            std::vector<double> decay, wall;
            for (const auto& m : model->modes)
                if (! m.high && m.shareAir > 0.5 && m.omega >= 2.0 * pi * edges[b] && m.omega < 2.0 * pi * edges[b + 1])
                {
                    decay.push_back (8.685889638 * m.alpha);
                    wall.push_back (m.alphaWall / m.alpha);
                }
            const auto median = [] (std::vector<double>& v)
            {
                if (v.empty())
                    return 0.0;
                std::nth_element (v.begin(), v.begin() + static_cast<std::ptrdiff_t> (v.size() / 2), v.end());
                return v[v.size() / 2];
            };
            out << (b > 0 ? "; " : "") << edges[b] / 1000.0 << "-" << edges[b + 1] / 1000.0 << " kHz: "
                << decay.size() << " modes, median decay " << median (decay) << " dB/s, median wall share "
                << median (wall);
        }
        out << ". Bulk atmospheric absorption is not modelled.\n";
    }
    out << "\n";

    std::vector<const CoupledMode*> order;
    for (const auto& m : model->modes)
        if (m.shareAir < 0.5)
            order.push_back (&m);
    std::sort (order.begin(), order.end(), [] (auto x, auto y) { return x->omega < y->omega; });
    out << "Lowest head-dominated coupled modes (energy shares batter/rear/air; each mode's own projected decay, "
           "parts hide/radiation/wall in 1/s; T60 = 6.91/alpha; last column: the left-capsule level of the whole "
           "drum body's boundary-element field against the transfer heard (baffled batter plus compact rear "
           "monopole), not applied; large for rear modes of order m >= 1, which that transfer barely hears):\n\n"
        << "| m | f (Hz) | T60 (s) | batter | rear | air | hide | radiation | wall | body at capsule (dB) |\n"
        << "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n";
    const double thetaS = std::atan2 (stroke.strikeY, stroke.strikeX);
    const double thetaMic = std::atan2 (mics.capsule[0][1], mics.capsule[0][0]);
    for (std::size_t i = 0; i < std::min<std::size_t> (order.size(), 16); ++i)
    {
        const auto& m = *order[i];
        out << "| " << m.m << " | " << m.omega / (2.0 * pi) << " | " << 6.907755 / m.alpha << " | " << m.shareBatter
            << " | " << m.shareRear << " | " << m.shareAir << " | " << m.alphaHide << " | " << m.alphaRadiation
            << " | " << m.alphaWall << " | ";
        const std::complex<double> heard = m.micBatter[0] * std::cos (m.m * (thetaMic - thetaS)) + m.micRear[0];
        const std::complex<double> body = m.micBody[0] * std::cos (m.m * (thetaMic - thetaS));
        if (body != 0.0 && std::abs (heard) > 0.0)
            out << 20.0 * std::log10 (std::abs (heard + body) / std::abs (heard));
        else
            out << "-";
        out << " |\n";
    }
    if (stepper)
    {
        const auto& r = stepper->report();
        const double mv = stroke.tipMass * stroke.impactSpeed;
        out << "\nContact (tip mass " << stroke.tipMass << " kg at " << stroke.impactSpeed << " m/s): "
            << r.episodes << " episode(s) from " << r.firstContact * 1000.0 << " to " << r.lastContact * 1000.0
            << " ms after the strike; peak force " << r.peakForce << " N; impulse " << r.impulse << " N s = "
            << r.impulse / mv << " m v0; stick leaves at " << r.exitVelocity << " m/s; max penetration "
            << r.maxPenetration * 1e6 << " um; " << r.oversampledSamples << " samples oversampled x"
            << runtime.oversample << ", " << r.baseRateContactSteps << " contact steps at the base rate.\n"
            << "Energy over the first contact window (stick + bank + Hertz + Berger SAV): start "
            << r.energyStart << " J, end " << r.energyEnd << " J + dissipated " << r.dissipated
            << " J; balance error " << (r.energyEnd + r.dissipated - r.energyStart) / r.energyStart << ".\n"
            << "Berger: peak tension rise " << r.peakTensionRise[0] << " N/m on the batter ("
            << 100.0 * r.peakTensionRise[0] / drum.batter.tension << " %), " << r.peakTensionRise[1]
            << " N/m on the rear; SAV drift |psi - sqrt(2 beta) S| <= " << r.maxPsiDrift << " of peak; "
            << bergerEntries << " of " << bergerDense << " head-row entries kept.\n";
        out << "Output: pole-residue sum plus the contact-force feedthrough and retardation kernel ("
            << feedthroughSpan * 1000.0 << " ms at the contact substep; integral " << feedthroughStatic[0] << " / "
            << feedthroughStatic[1] << " Pa/N at the left / right capsule): the pressure is causal, and its time "
            << "integral over the stroke vanishes but for the Berger force's own term (about 3 % of the integral "
            << "without the kernel).\n";
        out << "Cost: " << model->log << "; bank " << bankSeconds << " s (feedthrough kernel " << feedthroughSeconds
            << " s); render " << stepper->seconds() << " s (" << bankModes << " modes).\n";
    }
    return out.str();
}
} // namespace

std::unique_ptr<Engine> makeModalEngine() { return std::make_unique<ModalDrum>(); }

} // namespace taikor::candidates
