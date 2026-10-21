// Convergence of the self-consistent exterior loading with the frequency grid,
// and sanity of the assembled nagado model (passivity, ordering of the (0,1)
// pair against the literature's push-pull / push-push description).

#include "ModalBuild.h"
#include "ModalFixtures.h"
#include "ModalRuntime.h"
#include "ModalTestUtil.h"

#include <algorithm>
#include <complex>

using namespace taikor::candidates;
using namespace taikor::candidates::modal;

MODAL_TEST (build_grid_convergence)
{
    const DrumSpec d = modaltest::nagadoDrum();
    BuildOptions coarse;
    coarse.headMaxHz = coarse.coupledMaxHz;
    BuildOptions fine = coarse;
    fine.gridPerOctave = 2 * coarse.gridPerOctave;
    for (int m : { 0, 1, 6 })
    {
        const auto a = solveBlockForTest (d, coarse, m), b = solveBlockForTest (d, fine, m);
        double worstF = 0.0, worstA = 0.0;
        for (std::size_t k = 0; k < std::min (a.omega.size(), b.omega.size()); ++k)
        {
            if (a.shareAir[k] > 0.5)
                continue; // the exterior loads only head coordinates
            worstF = std::max (worstF, std::abs (a.omega[k] - b.omega[k]) / b.omega[k]);
            worstA = std::max (worstA, std::abs (a.alpha[k] - b.alpha[k]) / b.alpha[k]);
        }
        std::printf ("  m %d: %zu coupled modes\n", m, a.omega.size());
        modaltest::check (a.omega.size() == b.omega.size(), "same mode count", double (a.omega.size()), double (b.omega.size()), 0.0);
        // 3e-5 is 0.05 cent.
        modaltest::below ("grid 12 -> 24 per octave: max relative change of head-dominated frequencies", worstF, 3e-5);
        modaltest::below ("grid 12 -> 24 per octave: max relative change of their decay rates", worstA, 2e-3);
    }
}

MODAL_TEST (build_nagado_sanity)
{
    const DrumSpec d = modaltest::nagadoDrum();
    BuildOptions o;
    o.headMaxHz = 8000.0;
    std::string error;
    const auto model = buildModel (d, modaltest::nagadoMics(), o, error);
    modaltest::check (model != nullptr, "model built", model ? 1.0 : 0.0, 1.0, 0.0);
    double minAlpha = 1e300, maxShareError = 0.0;
    for (const auto& m : model->modes)
    {
        minAlpha = std::min (minAlpha, m.alpha);
        maxShareError = std::max (maxShareError, std::abs (m.shareBatter + m.shareRear + m.shareAir - 1.0));
    }
    modaltest::check (minAlpha > 0.0, "every coupled mode is passive (alpha > 0)", minAlpha, 0.0, 0.0);
    modaltest::below ("energy shares sum to one", maxShareError, 1e-10);
    // The (0,1) pair: the lower member moves both heads the same way (the air
    // column translates: push-pull in Suzuki & Hwang's terms, below the uncoupled
    // frequency), the upper compresses the barrel (push-push, above it).
    std::vector<const CoupledMode*> zero;
    for (const auto& m : model->modes)
        if (m.m == 0 && m.shareAir < 0.5)
            zero.push_back (&m);
    std::sort (zero.begin(), zero.end(), [] (auto x, auto y) { return x->omega < y->omega; });
    const double f0 = model->batter[0][0].omega / (2.0 * pi);
    std::printf ("  (0,1) pair %.2f / %.2f Hz around the in-vacuo batter (0,1) %.2f Hz\n", zero[0]->omega / (2.0 * pi),
                 zero[1]->omega / (2.0 * pi), f0);
    // Outward-normal coordinates: push-pull = opposite outward signs of the two heads' (0,1).
    const double pp = zero[0]->batterCoef[0] * zero[0]->rearCoef[0];
    const double ppush = zero[1]->batterCoef[0] * zero[1]->rearCoef[0];
    modaltest::check (pp < 0.0 && ppush > 0.0, "lower member push-pull, upper push-push (outward signs)", pp, ppush, 0.0);
    modaltest::check (zero[1]->omega > 2.0 * pi * f0 * 0.99, "push-push above the uncoupled (0,1)", zero[1]->omega / (2.0 * pi), f0, 0.0);
}

namespace
{
// Band energies (dB) of a capsule signal in frames of [t0, t1) s, by FFT masks
// over Hann frames of 20 ms (hop 10 ms): bands[b] to bands[b + 1] Hz.
std::vector<std::vector<double>> bandEnvelopes (const std::vector<double>& x, const std::vector<double>& bands,
                                                const std::vector<std::pair<double, double>>& windows)
{
    const int frame = 960, hop = 480;
    const std::size_t nb = bands.size() - 1;
    std::vector<std::vector<double>> out (windows.size(), std::vector<double> (nb, 0.0));
    for (int start = 0; start + frame <= static_cast<int> (x.size()); start += hop)
    {
        const double centre = (start + 0.5 * frame) / 48000.0;
        for (std::size_t w = 0; w < windows.size(); ++w)
        {
            if (centre < windows[w].first || centre >= windows[w].second)
                continue;
            for (int kbin = 1; kbin < frame / 2; ++kbin)
            {
                const double f = kbin * 48000.0 / frame;
                for (std::size_t b = 0; b < nb; ++b)
                {
                    if (f < bands[b] || f >= bands[b + 1])
                        continue;
                    std::complex<double> acc = 0.0;
                    for (int i = 0; i < frame; ++i)
                        acc += (0.5 - 0.5 * std::cos (2.0 * pi * i / frame)) * x[static_cast<std::size_t> (start + i)]
                             * std::polar (1.0, -2.0 * pi * kbin * i / frame);
                    out[w][b] += std::norm (acc);
                }
            }
        }
    }
    for (auto& row : out)
        for (auto& v : row)
            v = 10.0 * std::log10 (v + 1e-300);
    return out;
}
} // namespace

MODAL_TEST (build_coupled_band)
{
    // Batter modes above coupledMaxHz are solved on their own (self exterior
    // loading, evanescent interior loading), and barrel air modes exist only up
    // to airMaxHz. Moving the boundary from 3 / 4.5 kHz (the default) to
    // 4.5 / 6 kHz measures what that leaves out: (a) the lone batter modes'
    // frequencies against the coupled ones at 3-4.2 kHz; (b) the rendered
    // capsule signal in bands, clustered modes included.
    const DrumSpec d = modaltest::nagadoDrum();
    BuildOptions lone, coupled;
    lone.headMaxHz = coupled.headMaxHz = 8000.0;
    coupled.coupledMaxHz = 4500.0;
    coupled.airMaxHz = 6000.0;
    std::string error;
    const auto a = buildModel (d, modaltest::nagadoMics(), lone, error);
    const auto b = buildModel (d, modaltest::nagadoMics(), coupled, error);
    double worstF = 0.0;
    int count = 0;
    for (const auto& mode : a->modes)
    {
        if (! mode.high || mode.omega > 2.0 * pi * 4200.0)
            continue;
        const int index = mode.batterIndex[0];
        // The coupled mode carrying most of that in-vacuo mode.
        const CoupledMode* match = nullptr;
        double best = 0.0;
        for (const auto& other : b->modes)
            if (other.m == mode.m && ! other.high)
                for (std::size_t i = 0; i < other.batterIndex.size(); ++i)
                    if (other.batterIndex[i] == index && std::abs (other.batterCoef[i]) > best)
                    {
                        best = std::abs (other.batterCoef[i]);
                        match = &other;
                    }
        if (match == nullptr || match->shareBatter < 0.5)
            continue;
        worstF = std::max (worstF, std::abs (mode.omega - match->omega) / match->omega);
        ++count;
    }
    std::printf ("  %d batter modes at 3-4.2 kHz matched\n", count);
    modaltest::check (count > 100, "enough modes compared", double (count), 100.0, 0.0);
    modaltest::below ("lone vs coupled treatment: max relative frequency difference", worstF, 2e-3);

    // The side-wall loss sees the wall field of the explicit air modes only
    // (ModalCavity.h wallLoss); the wider band adds the next ones' wall pressure.
    // The 30 lowest head-dominated modes, matched by order and frequency.
    {
        std::vector<const CoupledMode*> low;
        for (const auto& mode : a->modes)
            if (! mode.high && mode.shareAir < 0.5)
                low.push_back (&mode);
        std::sort (low.begin(), low.end(), [] (auto x, auto y) { return x->omega < y->omega; });
        low.resize (std::min<std::size_t> (low.size(), 30));
        double worstWall = 0.0, worstAlpha = 0.0;
        for (const CoupledMode* mode : low)
        {
            const CoupledMode* match = nullptr;
            for (const auto& other : b->modes)
                if (other.m == mode->m && ! other.high
                    && (match == nullptr || std::abs (other.omega - mode->omega) < std::abs (match->omega - mode->omega)))
                    match = &other;
            worstWall = std::max (worstWall, std::abs (mode->alphaWall / match->alphaWall - 1.0));
            worstAlpha = std::max (worstAlpha, std::abs (mode->alphaWall - match->alphaWall) / mode->alpha);
        }
        modaltest::below ("30 lowest modes: max relative change of the wall loss with the wider air band", worstWall, 0.05);
        modaltest::below ("30 lowest modes: max change of the wall loss with the wider air band, relative to the decay",
                          worstAlpha, 1e-3);
    }

    std::vector<double> signal[2];
    for (int i = 0; i < 2; ++i)
    {
        Stepper stepper (makeBank (i == 0 ? *a : *b, d, modaltest::nagadoMics(), modaltest::nagadoStroke(), RuntimeOptions {}, 48000.0),
                         48000.0, RuntimeOptions {});
        stepper.strike();
        signal[i].resize (28800);
        stepper.processDouble (signal[i].data(), nullptr, 28800);
    }
    const std::vector<double> bands { 40.0, 250.0, 1000.0, 2000.0, 3000.0, 4500.0, 8000.0 };
    const std::vector<std::pair<double, double>> windows { { 0.0, 0.02 }, { 0.02, 0.05 }, { 0.05, 0.1 }, { 0.1, 0.2 },
                                                           { 0.2, 0.4 }, { 0.4, 0.6 } };
    const auto ea = bandEnvelopes (signal[0], bands, windows), eb = bandEnvelopes (signal[1], bands, windows);
    double below3k = 0.0;
    for (std::size_t band = 0; band + 1 < bands.size(); ++band)
    {
        std::printf ("  %4.0f-%4.0f Hz, default minus 4.5 / 6 kHz band (dB) in 0-20 / 20-50 / 50-100 / 100-200 / 200-400 / 400-600 ms:",
                     bands[band], bands[band + 1]);
        for (std::size_t w = 0; w < windows.size(); ++w)
        {
            const double diff = ea[w][band] - eb[w][band];
            std::printf (" %+.1f", diff);
            if (bands[band + 1] <= 3000.0)
                below3k = std::max (below3k, std::abs (diff));
        }
        std::printf ("\n");
    }
    modaltest::below ("below 3 kHz: max band-envelope change (dB)", below3k, 0.5);
}
