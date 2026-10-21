#include "ModalBodyExterior.h"
#include "ModalInternal.h"

#include <algorithm>
#include <chrono>
#include <cmath>

namespace taikor::candidates::modal
{

double bodyGridOmega (int g, int gridPerOctave) { return 2.0 * pi * std::pow (2.0, static_cast<double> (g) / gridPerOctave); }

BodyExterior::BodyExterior (const AirSpec& air, double wallThickness, const std::vector<std::vector<const HeadMode*>>& batter,
                            const std::vector<std::vector<const HeadMode*>>& rear, const BodyExteriorOptions& o)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int orders = static_cast<int> (std::max (batter.size(), rear.size()));
    const auto um = [] (int m) { return static_cast<std::size_t> (m); };
    heads.assign (um (orders), 0);
    batterCount.assign (um (orders), 0);
    rearCount.assign (um (orders), 0);
    for (int m = 0; m < orders; ++m)
    {
        batterCount[um (m)] = um (m) < batter.size() ? static_cast<int> (batter[um (m)].size()) : 0;
        rearCount[um (m)] = um (m) < rear.size() ? static_cast<int> (rear[um (m)].size()) : 0;
    }
    double b1Max = 0.0, lowest = 1e300;
    std::vector<double> lowestOf (um (orders), 1e300);
    for (int m = 0; m < orders; ++m)
    {
        for (const auto* list : { &batter, &rear })
            if (um (m) < list->size())
                for (const HeadMode* h : (*list)[um (m)])
                {
                    b1Max = std::max (b1Max, h->b1);
                    lowestOf[um (m)] = std::min (lowestOf[um (m)], h->omega);
                    ++heads[um (m)];
                }
        lowest = std::min (lowest, lowestOf[um (m)]);
    }
    // The global grid: w_g = 2 pi 2^(g/G) rad/s, from below half the lowest head mode to maxHz.
    const int gFirst = static_cast<int> (std::floor (o.gridPerOctave * std::log2 (0.45 * lowest / (2.0 * pi))));
    const int gLast = static_cast<int> (std::ceil (o.gridPerOctave * std::log2 (o.maxHz)));
    for (int g = gFirst; g <= gLast; ++g)
        omegas.push_back (bodyGridOmega (g, o.gridPerOctave));
    firstIndex = gFirst;
    gridPerOctave = o.gridPerOctave;

    // Panels: panelsPerWavelength over the acoustic wavelength at the top of the
    // table on the body, and over the shortest head-mode wavelength 2 pi / b1 on the heads.
    const double lambda = air.soundSpeed / o.maxHz;
    const double bodyStep = lambda / o.panelsPerWavelength;
    const double headStep = std::min (bodyStep, 2.0 * pi / std::max (b1Max, 1e-9) / o.panelsPerWavelength);
    BemSolver solver (drumGeometry (air.radius, wallThickness, air.length, headStep, bodyStep, o.chiefPoints));
    const auto& panels = solver.geometry().panels;
    panelCount = static_cast<int> (panels.size());

    // Per order: velocity at the head panels' nodes (batter modes on part 0,
    // rear modes on part 1) and the moments U_i[panel] = Lam_m sum_q w_q r_q phi_i(r_q).
    std::vector<double> gx, gw;
    special::gaussLegendre (bemVelocityNodes, gx, gw);
    std::vector<std::vector<double>> velocity (um (orders));
    std::vector<std::vector<double>> moments (um (orders)); // [panel * nh + i]
    for (int m = 0; m < orders; ++m)
    {
        const int nh = heads[um (m)];
        if (nh == 0)
            continue;
        const std::size_t nb = um (m) < batter.size() ? batter[um (m)].size() : 0;
        moments[um (m)].assign (panels.size() * static_cast<std::size_t> (nh), 0.0);
        for (std::size_t j = 0; j < panels.size(); ++j)
        {
            const auto& p = panels[j];
            if (p.part > 1)
                continue;
            const auto& list = p.part == 0 ? batter : rear;
            const std::size_t offset = p.part == 0 ? 0 : nb;
            for (int q = 0; q < bemVelocityNodes; ++q)
            {
                const double t = 0.5 * p.length * (gx[static_cast<std::size_t> (q)] + 1.0);
                const double r = p.r (t);
                const double w = 0.5 * p.length * gw[static_cast<std::size_t> (q)];
                std::vector<double> row (static_cast<std::size_t> (nh), 0.0);
                if (um (m) < list.size())
                    for (std::size_t i = 0; i < list[um (m)].size(); ++i)
                    {
                        const double value = evaluateShape (*list[um (m)][i], air.radius, std::min (r, air.radius)).value;
                        row[offset + i] = value;
                        moments[um (m)][j * static_cast<std::size_t> (nh) + offset + i] += angularNorm (m) * w * r * value;
                    }
                velocity[um (m)].insert (velocity[um (m)].end(), row.begin(), row.end());
            }
        }
    }

    // The body's lowest interior Dirichlet eigenvalue per order, k^2 = (j_m1 / b)^2 + (pi / L)^2:
    // below 0.9 of it the plain equation is regular and CHIEF rows are not needed.
    const double b = air.radius + wallThickness;
    const auto jZeros = special::besselJZeros (orders, 4.0 * orders + 10.0);
    std::vector<double> irregular (um (orders), 0.0);
    for (int m = 0; m < orders; ++m)
        irregular[um (m)] = std::hypot (jZeros[um (m)].front() / b, pi / air.length);

    table.assign (um (orders), std::vector<std::vector<std::complex<double>>> (omegas.size()));
    fields.assign (um (orders), std::vector<std::vector<std::complex<double>>> (omegas.size()));
    for (std::size_t g = 0; g < omegas.size(); ++g)
    {
        std::vector<BemProblem> problems;
        std::vector<int> order;
        for (int m = 0; m < orders; ++m)
            if (heads[um (m)] > 0 && omegas[g] >= 0.45 * lowestOf[um (m)])
            {
                BemProblem pr;
                pr.m = m;
                pr.chief = omegas[g] / air.soundSpeed > 0.9 * irregular[um (m)];
                pr.activeParts = { 0, 1 };
                pr.columns = heads[um (m)];
                pr.velocity = velocity[um (m)];
                pr.fieldPoints = o.fieldPoints;
                problems.push_back (std::move (pr));
                order.push_back (m);
            }
        if (problems.empty())
            continue;
        solver.solve (omegas[g] / air.soundSpeed, air.density, air.soundSpeed, problems, true, o.threads);
        if (! o.fieldPoints.empty())
            parallelFor (static_cast<int> (problems.size()), o.threads, [&] (int qi)
            {
                const auto q = static_cast<std::size_t> (qi);
                const int m = order[q];
                const int nh = heads[um (m)];
                auto& f = fields[um (m)][g];
                f = problems[q].field; // [point * nh + column]
                if (o.reference)
                    for (std::size_t point = 0; point < o.fieldPoints.size(); ++point)
                        for (int c = 0; c < nh; ++c)
                            f[point * static_cast<std::size_t> (nh) + static_cast<std::size_t> (c)]
                                -= o.reference (m, c, omegas[g], static_cast<int> (point));
            });
        for (std::size_t q = 0; q < problems.size(); ++q)
        {
            const int m = order[q];
            const int nh = heads[um (m)];
            auto& z = table[um (m)][g];
            z.assign (static_cast<std::size_t> (nh) * static_cast<std::size_t> (nh), 0.0);
            const auto& p = problems[q].pressure;
            const auto& u = moments[um (m)];
            for (std::size_t j = 0; j < panels.size(); ++j)
            {
                if (panels[j].part > 1)
                    continue;
                for (int i = 0; i < nh; ++i)
                {
                    const double ui = u[j * static_cast<std::size_t> (nh) + static_cast<std::size_t> (i)];
                    if (ui == 0.0)
                        continue;
                    for (int c = 0; c < nh; ++c)
                        z[static_cast<std::size_t> (i * nh + c)] += ui * p[j * static_cast<std::size_t> (nh) + static_cast<std::size_t> (c)];
                }
            }
            // Reciprocity makes Z symmetric; the discretisation leaves an asymmetry
            // of the order of its error (tests/modal body_exterior_convergence: it
            // shrinks with the mesh), removed here.
            for (int i = 0; i < nh; ++i)
                for (int c = i + 1; c < nh; ++c)
                {
                    const auto ic = static_cast<std::size_t> (i * nh + c), ci = static_cast<std::size_t> (c * nh + i);
                    const std::complex<double> mean = 0.5 * (z[ic] + z[ci]);
                    asymmetry = std::max (asymmetry, std::abs (z[ic] - z[ci]) / std::sqrt (std::abs (z[static_cast<std::size_t> (i * nh + i)] * z[static_cast<std::size_t> (c * nh + c)])));
                    z[ic] = z[ci] = mean;
                }
        }
    }
    buildSeconds = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
}

const std::vector<std::complex<double>>& BodyExterior::impedance (int m, int g) const
{
    static const std::vector<std::complex<double>> empty;
    if (m < 0 || static_cast<std::size_t> (m) >= table.size())
        return empty;
    const int index = g - firstIndex;
    if (index < 0 || static_cast<std::size_t> (index) >= omegas.size())
        return empty;
    return table[static_cast<std::size_t> (m)][static_cast<std::size_t> (index)];
}

bool BodyExterior::weights (int m, double omega, int& first, std::vector<double>& w) const
{
    w.clear();
    if (m < 0 || static_cast<std::size_t> (m) >= table.size() || omegas.size() < 2)
        return false;
    const double x = gridPerOctave * std::log2 (omega / (2.0 * pi)) - firstIndex;
    if (x < 0.0 || x > static_cast<double> (omegas.size() - 1))
        return false;
    const int lo = std::clamp (static_cast<int> (std::floor (x)), 0, static_cast<int> (omegas.size()) - 2);
    const auto& row = table[static_cast<std::size_t> (m)];
    if (row[static_cast<std::size_t> (lo)].empty() || row[static_cast<std::size_t> (lo + 1)].empty())
        return false;
    // Cubic Lagrange interpolation in log w through the four surrounding grid
    // points where the table has them (error fourth order in the grid step,
    // tests/modal build_grid_convergence), linear at the table's ends.
    first = lo - 1;
    int last = lo + 2;
    if (first < 0 || row[static_cast<std::size_t> (first)].empty() || last >= static_cast<int> (omegas.size())
        || row[static_cast<std::size_t> (last)].empty())
    {
        first = lo;
        last = lo + 1;
    }
    for (int p = first; p <= last; ++p)
    {
        double weight = 1.0;
        for (int q = first; q <= last; ++q)
            if (q != p)
                weight *= (x - q) / static_cast<double> (p - q);
        w.push_back (weight);
    }
    return true;
}

void BodyExterior::interpolate (int m, double omega, bool real, std::vector<double>& out) const
{
    out.clear();
    int first = 0;
    std::vector<double> w;
    if (! weights (m, omega, first, w))
        return;
    const auto& row = table[static_cast<std::size_t> (m)];
    out.assign (row[static_cast<std::size_t> (first)].size(), 0.0);
    for (std::size_t p = 0; p < w.size(); ++p)
    {
        const auto& z = row[static_cast<std::size_t> (first) + p];
        const double wg = omegas[static_cast<std::size_t> (first) + p];
        for (std::size_t i = 0; i < z.size(); ++i)
            out[i] += w[p] * (real ? z[i].real() : z[i].imag() / wg);
    }
}

void BodyExterior::fieldCorrection (int m, double omega, std::vector<std::complex<double>>& out) const
{
    out.clear();
    int first = 0;
    std::vector<double> w;
    if (! weights (m, omega, first, w) || fields[static_cast<std::size_t> (m)][static_cast<std::size_t> (first)].empty())
        return;
    const auto& row = fields[static_cast<std::size_t> (m)];
    out.assign (row[static_cast<std::size_t> (first)].size(), 0.0);
    for (std::size_t p = 0; p < w.size(); ++p)
    {
        const auto& f = row[static_cast<std::size_t> (first) + p];
        for (std::size_t i = 0; i < f.size(); ++i)
            out[i] += w[p] * f[i];
    }
}

void BodyExterior::addedMass (int m, double omega, std::vector<double>& out) const { interpolate (m, omega, false, out); }
void BodyExterior::resistance (int m, double omega, std::vector<double>& out) const { interpolate (m, omega, true, out); }

} // namespace taikor::candidates::modal
