// The per-stroke bank: what the strike position and the capsules select from
// the precomputed coupled modes (ModalBuild.h).
//
// - Contact projection c_k = sum_i x_ki <phi_i>_disk: the exact mean of each
//   batter shape over the uniform footprint disk (ModalMath.h diskAverage), so
//   the same disk average senses the head and spreads the force.
// - Capsule transfer H_k = micBatter cos(m (theta_mic - theta_strike)) + micRear.
// - Pruning by residue only: a mode is dropped when both its strike-to-capsule
//   residue |H_k| |c_k| / w_k (the pressure amplitude of its impulse response)
//   is below pruneResidue of the largest one and its footprint projection is
//   below pruneContact of the largest one, i.e. it is neither heard nor driven
//   (numerical tolerances; tests/modal modal_pruning).
// - Berger rows: each head mode's coordinate q_i = sum_k x_ki y_k, keeping the
//   entries whose size sqrt(d_i) |x_ki| |c_k| / w_k (the gradient amplitude the
//   entry carries per unit strike impulse) is at least bergerTolerance of the
//   head's largest (numerical; tests/modal modal_berger_sparsity). Dropping an
//   entry changes the form, not its structure: the sparsified S is still
//   Gamma's quadratic form of some head displacement, >= 0, so the SAV scheme
//   stays exactly energy-conserving.

#include "ModalBuild.h"
#include "ModalEigen.h"
#include "ModalRuntime.h"

#include <functional>

#include <algorithm>
#include <cmath>

namespace taikor::candidates::modal
{

Bank makeBank (const ModalModel& model, const DrumSpec& drum, const MicSpec& mics, const StrokeSpec& stroke,
               const RuntimeOptions& o, double sampleRate)
{
    Bank bank;
    const double a = model.radius;
    const double thetaS = std::atan2 (stroke.strikeY, stroke.strikeX);
    const double rs = std::hypot (stroke.strikeX, stroke.strikeY);
    const double b = stroke.footprintRadius;
    std::vector<std::vector<double>> disk (model.batter.size());
    for (std::size_t m = 0; m < model.batter.size(); ++m)
        for (const auto& mode : model.batter[m])
            disk[m].push_back (diskAverage (mode, a, rs, b));
    double thetaMic[2];
    for (int j = 0; j < 2; ++j)
        thetaMic[j] = std::atan2 (mics.capsule[static_cast<std::size_t> (j)][1], mics.capsule[static_cast<std::size_t> (j)][0]);

    const std::size_t total = model.modes.size();
    std::vector<BankMode> all (total);
    double maxResidue = 0.0, maxContact = 0.0;
    std::vector<double> residue (total);
    for (std::size_t k = 0; k < total; ++k)
    {
        const CoupledMode& cm = model.modes[k];
        BankMode& bm = all[k];
        bm.omega = cm.omega;
        bm.alpha = cm.alpha;
        const auto um = static_cast<std::size_t> (cm.m);
        for (std::size_t i = 0; i < cm.batterIndex.size(); ++i)
            bm.contact += cm.batterCoef[i] * disk[um][static_cast<std::size_t> (cm.batterIndex[i])];
        double heard = 0.0;
        for (int j = 0; j < 2; ++j)
        {
            const double angular = std::cos (cm.m * (thetaMic[j] - thetaS));
            bm.mic[static_cast<std::size_t> (j)] = cm.micBatter[static_cast<std::size_t> (j)] * angular
                                                 + cm.micRear[static_cast<std::size_t> (j)];
            heard = std::max (heard, std::abs (bm.mic[static_cast<std::size_t> (j)]));
        }
        residue[k] = heard * std::abs (bm.contact) / std::max (cm.omega, 1e-9);
        maxResidue = std::max (maxResidue, residue[k]);
        maxContact = std::max (maxContact, std::abs (bm.contact));
    }
    std::vector<int> index (total, -1);
    for (std::size_t k = 0; k < total; ++k)
    {
        const bool keep = residue[k] >= o.pruneResidue * maxResidue || std::abs (all[k].contact) >= o.pruneContact * maxContact;
        if (! keep)
        {
            ++bank.prunedModes;
            continue;
        }
        index[k] = static_cast<int> (bank.modes.size());
        bank.modes.push_back (all[k]);
    }

    // Damping clusters: connected components of the kept couplings between kept
    // modes. If dropping the weak couplings leaves a decay matrix that is not
    // positive semi-definite, its couplings are shrunk by the least common factor
    // that restores it (the diagonals, each mode's own decay, are kept), so the
    // damping can never inject energy.
    {
        std::vector<int> parent (bank.modes.size());
        for (std::size_t i = 0; i < parent.size(); ++i)
            parent[i] = static_cast<int> (i);
        const std::function<int (int)> find = [&] (int x) { return parent[static_cast<std::size_t> (x)] == x ? x : parent[static_cast<std::size_t> (x)] = find (parent[static_cast<std::size_t> (x)]); };
        for (std::size_t k = 0; k < total; ++k)
            if (index[k] >= 0)
                for (const auto& [other, value] : model.modes[k].damping)
                    if (index[static_cast<std::size_t> (other)] >= 0)
                        parent[static_cast<std::size_t> (find (index[k]))] = find (index[static_cast<std::size_t> (other)]);
        std::vector<std::vector<int>> groups (bank.modes.size());
        for (std::size_t i = 0; i < parent.size(); ++i)
            groups[static_cast<std::size_t> (find (static_cast<int> (i)))].push_back (static_cast<int> (i));
        std::vector<int> original (bank.modes.size());
        for (std::size_t k = 0; k < total; ++k)
            if (index[k] >= 0)
                original[static_cast<std::size_t> (index[k])] = static_cast<int> (k);
        for (auto& members : groups)
        {
            if (members.size() < 2)
                continue;
            DampingCluster cl;
            const std::size_t p = members.size();
            cl.members = members;
            Matrix a (static_cast<int> (p));
            for (std::size_t i = 0; i < p; ++i)
            {
                const CoupledMode& cm = model.modes[static_cast<std::size_t> (original[static_cast<std::size_t> (members[i])])];
                a (static_cast<int> (i), static_cast<int> (i)) = cm.alpha;
                for (const auto& [other, value] : cm.damping)
                {
                    const int bankIndex = index[static_cast<std::size_t> (other)];
                    const auto at = std::find (members.begin(), members.end(), bankIndex);
                    if (bankIndex >= 0 && at != members.end())
                        a (static_cast<int> (i), static_cast<int> (at - members.begin())) = value;
                }
            }
            // Keep every diagonal (each mode's own projected decay) and shrink the
            // couplings by the least factor that makes the matrix PSD, if needed.
            std::vector<double> values;
            Matrix vectors;
            if (symmetricEigen (a, values, vectors) && values.front() < 0.0)
            {
                ++bank.clippedClusters;
                double lo = 0.0, hi = 1.0;
                const auto scaled = [&] (double f)
                {
                    Matrix s = a;
                    for (int i = 0; i < s.n; ++i)
                        for (int j = 0; j < s.n; ++j)
                            if (i != j)
                                s (i, j) *= f;
                    return s;
                };
                for (int iteration = 0; iteration < 50; ++iteration)
                {
                    const double mid = 0.5 * (lo + hi);
                    std::vector<double> v;
                    Matrix vv;
                    symmetricEigen (scaled (mid), v, vv);
                    (v.front() >= 0.0 ? lo : hi) = mid;
                }
                a = scaled (lo);
            }
            cl.decay = a.a;
            bank.clusteredModes += static_cast<int> (p);
            bank.clusters.push_back (std::move (cl));
        }
    }

    // Berger rows per head and angular order.
    const HeadSpec* specs[2] = { &drum.batter, &drum.rear };
    const std::vector<std::vector<HeadMode>>* lists[2] = { &model.batter, &model.rear };
    const std::vector<GradientForm>* forms[2] = { &model.batterGradient, &model.rearGradient };
    for (int h = 0; h < 2; ++h)
    {
        const HeadSpec& s = *specs[h];
        const double stretch = s.youngsModulus * s.thickness / ((1.0 - s.poissonRatio * s.poissonRatio) * pi * a * a);
        BergerHead& head = bank.berger[static_cast<std::size_t> (h)];
        head.beta = stretch / 8.0;
        head.tensionPerS = stretch / 2.0;
        const auto& list = *lists[h];
        head.groups.resize (list.size());
        for (std::size_t m = 0; m < list.size(); ++m)
        {
            auto& g = head.groups[m];
            const auto& form = (*forms[h])[m];
            g.diagonal = form.diagonal;
            // Any prefix of a pivoted Cholesky leaves a PSD residual, so dropping the
            // trailing vectors keeps Gamma_approx >= Gamma; the runtime keeps the
            // pivots above bergerRankTolerance (tests/modal modal_berger_sparsity).
            for (std::size_t r = 0; r < form.lowRank.size(); ++r)
                if (form.pivots[r] > o.bergerRankTolerance)
                    g.lowRank.push_back (form.lowRank[r]);
            g.rows.assign (list[m].size(), {});
        }
        // Entry size: the head mode's gradient amplitude sqrt(d_i) |x_ki| carried by
        // coupled mode k per unit strike impulse (|c_k| / w_k, its impulse-response
        // amplitude), relative to the largest entry of this head.
        double largest = 0.0;
        const auto entrySize = [&] (std::size_t k, std::size_t i, const GradientForm& form, int headIndex)
        {
            const CoupledMode& cm = model.modes[k];
            const auto& coef = h == 0 ? cm.batterCoef : cm.rearCoef;
            return std::abs (coef[i]) * std::sqrt (form.diagonal[static_cast<std::size_t> (headIndex)])
                 * std::abs (all[k].contact) / std::max (cm.omega, 1e-9);
        };
        for (int pass = 0; pass < 2; ++pass)
            for (std::size_t k = 0; k < total; ++k)
            {
                if (index[k] < 0)
                    continue;
                const CoupledMode& cm = model.modes[k];
                const auto& idx = h == 0 ? cm.batterIndex : cm.rearIndex;
                const auto& coef = h == 0 ? cm.batterCoef : cm.rearCoef;
                if (idx.empty())
                    continue;
                const auto& form = (*forms[h])[static_cast<std::size_t> (cm.m)];
                auto& g = head.groups[static_cast<std::size_t> (cm.m)];
                for (std::size_t i = 0; i < idx.size(); ++i)
                {
                    const double size = entrySize (k, i, form, idx[i]);
                    if (pass == 0)
                    {
                        largest = std::max (largest, size);
                        continue;
                    }
                    ++bank.bergerDenseEntries;
                    if (size < o.bergerTolerance * largest)
                        continue;
                    g.rows[static_cast<std::size_t> (idx[i])].emplace_back (index[k], coef[i]);
                    ++bank.bergerEntries;
                }
            }
    }
    bank.tipMass = stroke.tipMass;
    bank.hertzStiffness = stroke.hertzStiffness;
    bank.hertzExponent = stroke.hertzExponent;
    bank.impactSpeed = stroke.impactSpeed;
    if (sampleRate > 0.0 && o.feedthrough)
    {
        std::vector<int> original (bank.modes.size());
        for (std::size_t k = 0; k < total; ++k)
            if (index[k] >= 0)
                original[static_cast<std::size_t> (index[k])] = static_cast<int> (k);
        buildFeedthrough (bank, model, drum, mics, stroke, original, 1.0 / (sampleRate * std::max (1, o.oversample)),
                          o.threads);
    }
    return bank;
}

} // namespace taikor::candidates::modal
