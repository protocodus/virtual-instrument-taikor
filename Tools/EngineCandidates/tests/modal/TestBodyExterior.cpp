// The drum body's exterior in the coupled blocks (Modal/ModalBodyExterior.h):
// its numerical parameters against refined values on a reduced model (heads
// coupled to 2 kHz), and its effect on the lowest modes against the per-head
// baffled King loading it replaces.

#include "ModalBem.h"
#include "ModalBuild.h"
#include "ModalFixtures.h"
#include "ModalTestUtil.h"

#include <algorithm>
#include <cmath>
#include <complex>

using namespace taikor::candidates;
using namespace taikor::candidates::modal;

namespace
{
BuildOptions reduced()
{
    BuildOptions o;
    o.headMaxHz = 4000.0;
    o.coupledMaxHz = 2000.0;
    o.airMaxHz = 3000.0;
    return o;
}

// Largest relative frequency and decay-rate change of the coupled modes in
// [lo, hi) Hz between two builds of the same structure (modes compared by index).
struct Change
{
    double frequency = 0.0, decay = 0.0;
    int modes = 0;
};
Change compare (const ModalModel& a, const ModalModel& b, double lo, double hi)
{
    Change c;
    if (a.modes.size() != b.modes.size())
    {
        c.frequency = c.decay = 1e300;
        return c;
    }
    for (std::size_t k = 0; k < a.modes.size(); ++k)
    {
        const auto& x = a.modes[k];
        const auto& y = b.modes[k];
        const double f = x.omega / (2.0 * pi);
        if (x.high || f < lo || f >= hi)
            continue;
        c.frequency = std::max (c.frequency, std::abs (x.omega - y.omega) / y.omega);
        c.decay = std::max (c.decay, std::abs (x.alpha - y.alpha) / y.alpha);
        ++c.modes;
    }
    return c;
}
} // namespace

MODAL_TEST (body_exterior_convergence)
{
    const DrumSpec d = modaltest::nagadoDrum();
    std::string error;
    BuildOptions full = reduced();
    full.bodyMaxHz = 0.0;                       // the table over the blocks' whole range (2.6 kHz)
    BuildOptions top = reduced();               // the default top, 1.5 kHz
    BuildOptions mesh = reduced();
    mesh.bodyPanelsPerWavelength = 16.0;
    const auto ref = buildModel (d, modaltest::nagadoMics(), full, error);
    const auto base = buildModel (d, modaltest::nagadoMics(), top, error);
    const auto fine = buildModel (d, modaltest::nagadoMics(), mesh, error);
    std::printf ("  body table: %d panels to %.0f Hz (%.1f s); %d panels to %.0f Hz (%.1f s); %d panels (%.1f s)\n",
                 ref->bodyPanels, ref->bodyMaxHz, ref->bodySeconds, base->bodyPanels, base->bodyMaxHz, base->bodySeconds,
                 fine->bodyPanels, fine->bodySeconds);
    // Reciprocity: the impedance table is symmetric in the continuum; the
    // collocation leaves an asymmetry of the order of its error, which
    // ModalBodyExterior.cpp removes. It must shrink as the mesh is refined.
    std::printf ("  impedance asymmetry before symmetrising: %.2e (10 panels/wavelength), %.2e (16); %.2e (table to 2.6 kHz)\n",
                 base->bodyAsymmetry, fine->bodyAsymmetry, ref->bodyAsymmetry);
    modaltest::check (fine->bodyAsymmetry < base->bodyAsymmetry, "impedance asymmetry shrinks from 10 to 16 panels per wavelength",
                      fine->bodyAsymmetry, base->bodyAsymmetry, 0.0);
    for (const auto& [lo, hi] : { std::pair<double, double> { 0.0, 1000.0 }, std::pair<double, double> { 1000.0, 2000.0 },
                                  std::pair<double, double> { 2000.0, 3000.0 } })
    {
        const Change t = compare (*base, *ref, lo, hi), p = compare (*base, *fine, lo, hi);
        std::printf ("  %4.0f-%4.0f Hz (%d modes): top 1.5 kHz vs 2.6 kHz: df/f %.1e, da/a %.4f; 10 vs 16 panels/wavelength: df/f %.1e, da/a %.4f\n",
                     lo, hi, t.modes, t.frequency, t.decay, p.frequency, p.decay);
        if (hi <= 1000.0)
        {
            modaltest::below ("below 1 kHz, table top 1.5 vs 2.6 kHz: max relative decay change", t.decay, 1e-2);
            modaltest::below ("below 1 kHz, 10 vs 16 panels per wavelength: max relative decay change", p.decay, 1e-2);
        }
        modaltest::below ("table top 1.5 vs 2.6 kHz: max relative frequency change", t.frequency, 2e-4);
        modaltest::below ("10 vs 16 panels per wavelength: max relative frequency change", p.frequency, 2e-4);
        modaltest::below ("table top 1.5 vs 2.6 kHz: max relative decay change", t.decay, 0.08);
    }
}

MODAL_TEST (body_exterior_lowest_modes)
{
    // The (0,1) pair against the per-head baffled loading: one exterior field
    // lets the push-pull member (both heads moving the same way, opposite
    // outward volumes) radiate as a dipole, far less than two baffled pistons,
    // while the push-push member (a monopole of the whole body) radiates into
    // 4 pi rather than two 2 pi half-spaces.
    const DrumSpec d = modaltest::nagadoDrum();
    std::string error;
    BuildOptions body = reduced(), king = reduced();
    king.bodyExterior = false;
    const auto a = buildModel (d, modaltest::nagadoMics(), body, error);
    const auto b = buildModel (d, modaltest::nagadoMics(), king, error);
    const auto pair = [] (const ModalModel& model)
    {
        std::vector<const CoupledMode*> zero;
        for (const auto& m : model.modes)
            if (m.m == 0 && m.shareAir < 0.5)
                zero.push_back (&m);
        std::sort (zero.begin(), zero.end(), [] (auto x, auto y) { return x->omega < y->omega; });
        return std::pair { zero[0], zero[1] };
    };
    const auto [pa, qa] = pair (*a);
    const auto [pb, qb] = pair (*b);
    std::printf ("  push-pull %.2f Hz: radiation %.3f 1/s (baffled heads %.3f), T60 %.2f s (%.2f s)\n", pa->omega / (2.0 * pi),
                 pa->alphaRadiation, pb->alphaRadiation, 6.907755 / pa->alpha, 6.907755 / pb->alpha);
    std::printf ("  push-push %.2f Hz: radiation %.3f 1/s (baffled heads %.3f), T60 %.2f s (%.2f s)\n", qa->omega / (2.0 * pi),
                 qa->alphaRadiation, qb->alphaRadiation, 6.907755 / qa->alpha, 6.907755 / qb->alpha);
    modaltest::below ("push-pull radiation, body / baffled heads", pa->alphaRadiation / pb->alphaRadiation, 0.5);
    modaltest::check (qa->alphaRadiation > 2.0 * pa->alphaRadiation, "push-push radiates more than push-pull",
                      qa->alphaRadiation, pa->alphaRadiation, 0.0);
}

MODAL_TEST (body_capsule_field)
{
    // The body's addition to the capsule transfer (CoupledMode::micBody, a
    // diagnostic) against a direct boundary-element solve at the mode's own
    // frequency with the mode's full velocity on both heads: checks the table's
    // field points, their interpolation and the reference subtraction.
    const DrumSpec d = modaltest::nagadoDrum();
    const MicSpec mics = modaltest::nagadoMics();
    std::string error;
    const auto model = buildModel (d, mics, reduced(), error);
    std::vector<const CoupledMode*> list;
    for (const auto& cm : model->modes)
        if (cm.m <= 2 && cm.shareAir < 0.5 && cm.omega < 2.0 * pi * 300.0)
            list.push_back (&cm);
    std::sort (list.begin(), list.end(), [] (auto x, auto y) { return x->omega < y->omega; });
    BemSolver solver (drumGeometry (d.radius, d.shellWallThickness, d.depth, 0.004, 0.02, 40));
    double worst = 0.0;
    for (const auto* cm : list)
    {
        BemProblem pr;
        pr.m = cm->m;
        pr.activeParts = { 0, 1 };
        pr.columns = 1;
        for (const auto& p : solver.geometry().panels)
        {
            if (p.part > 1)
                continue;
            for (double t : bemNodeOffsets (p))
            {
                const double r = std::min (p.r (t), d.radius);
                double u = 0.0;
                const auto& index = p.part == 0 ? cm->batterIndex : cm->rearIndex;
                const auto& coef = p.part == 0 ? cm->batterCoef : cm->rearCoef;
                const auto& heads = p.part == 0 ? model->batter : model->rear;
                for (std::size_t i = 0; i < index.size(); ++i)
                    u += coef[i] * evaluateShape (heads[static_cast<std::size_t> (cm->m)][static_cast<std::size_t> (index[i])], d.radius, r).value;
                pr.velocity.push_back (u);
            }
        }
        pr.fieldPoints.push_back ({ std::hypot (mics.capsule[0][0], mics.capsule[0][1]), mics.capsule[0][2] });
        std::vector<BemProblem> problems { pr };
        solver.solve (cm->omega / d.soundSpeed, d.airDensity, d.soundSpeed, problems, true, 8);
        const std::complex<double> s (0.0, cm->omega);
        const std::complex<double> direct = s * problems[0].field[0];
        const std::complex<double> table = cm->micBatter[0] + cm->micRear[0] + cm->micBody[0];
        const double error = std::abs (table - direct) / std::abs (direct);
        std::printf ("  m %d, %.1f Hz: capsule |H| heard %.4g, with the body %.4g, direct solve %.4g (difference %.1e)\n", cm->m,
                     cm->omega / (2.0 * pi), std::abs (cm->micBatter[0] + cm->micRear[0]), std::abs (table), std::abs (direct), error);
        if (std::abs (direct) > 1e-3 * std::abs (cm->micBatter[0]))
            worst = std::max (worst, error);
    }
    modaltest::below ("body-corrected capsule transfer vs direct solve, lowest modes (relative)", worst, 0.05);
}
