// (3) Rigid heads recover the rigid-cylinder air modes; (4) the Suzuki & Hwang
// 2008 two-head coupling fixture; convergence of the explicit air band and of
// the cross-mode truncation.

#include "ModalBuild.h"
#include "ModalFixtures.h"
#include "ModalTestUtil.h"

#include <algorithm>

using namespace taikor::candidates;
using namespace taikor::candidates::modal;

MODAL_TEST (cavity_rigid_heads)
{
    DrumSpec d = modaltest::nagadoDrum();
    // Heads a million times heavier: they cannot follow the air, so the
    // air-dominated coupled modes must sit on c sqrt(gamma^2 + (p pi/L)^2).
    d.batter.arealDensity *= 1e6;
    d.rear.arealDensity *= 1e6;
    BuildOptions o;
    o.exterior = false;
    o.headMaxHz = o.coupledMaxHz = 3.0;   // heads now resonate below 1 Hz
    o.airMaxHz = 1500.0;
    for (int m : { 0, 1, 3 })
    {
        const auto block = solveBlockForTest (d, o, m);
        std::vector<double> airFreq;
        for (std::size_t k = 0; k < block.omega.size(); ++k)
            if (block.shareAir[k] > 0.5)
                airFreq.push_back (block.omega[k]);
        const auto cross = crossModes (m, d.radius, 2.0 * pi * o.airMaxHz / d.soundSpeed + 40.0 / d.radius);
        AirSpec spec { d.radius, d.depth, d.airDensity, d.soundSpeed };
        const auto air = airModes (m, cross, spec, 2.0 * pi * o.airMaxHz);
        std::vector<double> exact;
        for (const auto& a : air)
            exact.push_back (a.omega);
        std::sort (exact.begin(), exact.end());
        std::printf ("  m %d: %zu air-dominated modes, %zu rigid-cylinder modes\n", m, airFreq.size(), exact.size());
        double worst = 0.0;
        for (std::size_t k = 0; k < std::min (airFreq.size(), exact.size()); ++k)
            worst = std::max (worst, std::abs (airFreq[k] - exact[k]) / exact[k]);
        modaltest::check (airFreq.size() == exact.size(), "same number of air modes", double (airFreq.size()), double (exact.size()), 0.0);
        modaltest::below ("max relative deviation from rigid-cylinder frequencies", worst, 1e-5);
    }
}

MODAL_TEST (cavity_suzuki_hwang)
{
    // Suzuki & Hwang neglect the hide's bending stiffness ("the effect of bending
    // stiffness will be neglected"), check against FEM with E = 1e6 Pa, and use
    // rho = 1.21 kg/m^3; the fixture follows them.
    DrumSpec d = modaltest::suzukiHwangDrum();
    d.airDensity = 1.21;
    for (auto* h : { &d.batter, &d.rear })
        h->bendingStiffness = 1.0e6 * 0.002 * 0.002 * 0.002 / (12.0 * (1.0 - 0.09));
    BuildOptions o;
    o.exterior = false;
    o.coupledMaxHz = 3000.0;
    o.headMaxHz = 3000.0;
    o.airMaxHz = 4500.0;
    const auto block = solveBlockForTest (d, o, 0);
    // The two lowest head-dominated modes are the (0,1) push-pull and push-push pair.
    std::vector<double> pair;
    for (std::size_t k = 0; k < block.omega.size() && pair.size() < 2; ++k)
        if (block.shareAir[k] < 0.5)
            pair.push_back (block.omega[k] / (2.0 * pi));
    std::printf ("  interior only: %.2f / %.2f Hz (Suzuki-Hwang bulged barrel 110.5 / 120.0; shipping cylinder solver 110.1 / 122.6)\n",
                 pair[0], pair[1]);
    modaltest::nearRel ("push-pull vs the shipping cylinder solver 110.1 Hz", pair[0], 110.1, 0.015);
    modaltest::nearRel ("push-push vs the shipping cylinder solver 122.6 Hz", pair[1], 122.6, 0.015);
    modaltest::nearRel ("push-pull vs Suzuki-Hwang 110.5 Hz", pair[0], 110.5, 0.03);
    modaltest::nearRel ("push-push vs Suzuki-Hwang 120.0 Hz (bulged barrel: larger volume)", pair[1], 120.0, 0.03);
    o.exterior = true;
    const auto loaded = solveBlockForTest (d, o, 0);
    std::vector<double> lp;
    for (std::size_t k = 0; k < loaded.omega.size() && lp.size() < 2; ++k)
        if (loaded.shareAir[k] < 0.5)
            lp.push_back (loaded.omega[k] / (2.0 * pi));
    std::printf ("  with the baffled exterior too: %.2f / %.2f Hz\n", lp[0], lp[1]);
}

MODAL_TEST (cavity_convergence)
{
    const DrumSpec d = modaltest::nagadoDrum();
    BuildOptions o;
    o.exterior = false;
    o.coupledMaxHz = 3000.0;
    o.headMaxHz = 3000.0;
    // Every coupled mode, matched across runs by nearest frequency (head/air
    // classification of strongly mixed pairs is not stable, their frequencies are).
    const auto headModes = [&] (const BuildOptions& opt, int m) { return solveBlockForTest (d, opt, m).omega; };
    const auto compare = [] (const std::vector<double>& a, const std::vector<double>& b, double below, double above)
    {
        double worst = 0.0;
        for (double x : a)
        {
            if (x < 2.0 * pi * below || x >= 2.0 * pi * above)
                continue;
            double best = 1e300;
            for (double y : b)
                best = std::min (best, std::abs (x - y) / y);
            worst = std::max (worst, best);
        }
        return worst;
    };
    for (int m : { 0, 2 })
    {
        BuildOptions lo = o, hi = o, top = o;
        lo.airMaxHz = 4500.0;
        hi.airMaxHz = 6000.0;
        top.airMaxHz = 9000.0;
        const auto a = headModes (lo, m), b = headModes (hi, m), t = headModes (top, m);
        std::printf ("  m %d: %zu coupled modes\n", m, a.size());
        for (double band : { 1000.0, 2000.0, 3000.0 })
            std::printf ("    to %.0f Hz: max change 4.5 -> 9 kHz air band %.2e, 6 -> 9 kHz %.2e\n", band,
                         compare (a, t, 0.0, band), compare (b, t, 0.0, band));
        modaltest::below ("air band 4.5 -> 9 kHz: max relative change of coupled modes < 1 kHz", compare (a, t, 0.0, 1000.0), 1e-4);
        modaltest::below ("air band 4.5 -> 9 kHz: max relative change of coupled modes < 3 kHz", compare (a, t, 0.0, 3000.0), 5e-4);
        BuildOptions wide = lo;
        wide.crossGammaFactor = 6.0;
        modaltest::below ("cross-mode truncation 3 b1 -> 6 b1: max relative change", compare (a, headModes (wide, m), 0.0, 3000.0), 1e-7);
    }
}
