// Berger's average-stretching nonlinearity for a clamped head with an
// immovable edge (H. M. Berger, "A new approach to the analysis of large
// deflections of plates", J. Appl. Mech. 22, 1955, 465-472), in its exact
// form: the membrane force of a uniform in-plane stretch,
//   N = E h / (2 (1 - nu^2) pi a^2) * int |grad w|^2 dA,
// stores the potential
//   V = xi / 4 (int |grad w|^2 dA)^2,  xi = E h / (2 (1 - nu^2) pi a^2)
// (Bilbao & Webb, JAES 61(10), 2013, eqs. 14-16; Torin, PhD thesis 2016,
// sec. 3.3). The shipping engine scales this by a calibrated 0.10; here the
// coefficient is the exact one (shared decision 5).
//
// Discretely, int |grad w|^2 dA = sum over grid edges of (w_a - w_b)^2 (each
// edge stands for an area h^2 and a slope (w_a - w_b)/h), with w = 0 on the
// clamped nodes, so V is an exact function of the nodal displacements and its
// gradient is exact:
//   dV/dw_l = xi G * sum_{neighbours n} (w_l - w_n) = -xi G h^2 (L w)_l,
// where G is the edge sum. The SAV step (FdtdDrum.cpp) needs V and dV/dw only.

#include "FdtdContracts.h"
#include "StretchingProbe.h"

#include <cmath>
#include <vector>

namespace taikor::candidates
{
namespace
{
class BergerPotential final : public NonlinearPotential, public fdtd::StretchingProbe
{
public:
    std::string name() const override { return "Berger (exact, immovable edge)"; }

    bool prepare (const MembraneGrid& grid, const HeadSpec& head, double headRadius, std::string& error) override
    {
        mesh = grid;
        const double nu = head.poissonRatio;
        xi = head.youngsModulus * head.thickness
             / (2.0 * (1.0 - nu * nu) * 3.14159265358979323846 * headRadius * headRadius);
        if (! (xi > 0.0) || grid.nx < 5)
        {
            error = "Berger potential: invalid head or grid";
            return false;
        }
        // Rows that hold interior nodes, and the span to visit in each (one
        // node of margin so every edge touching the interior is seen).
        spans.clear();
        for (int j = 1; j < grid.ny - 1; ++j)
        {
            int first = -1, last = -1;
            for (int i = 0; i < grid.nx; ++i)
                if (grid.interior[static_cast<std::size_t> (grid.index (i, j))])
                {
                    if (first < 0)
                        first = i;
                    last = i;
                }
            if (first >= 0)
                spans.push_back ({ j, first, last + 1 });
        }
        return true;
    }

    double evaluate (const double* w, double* gradient) override
    {
        const int nx = mesh.nx;
        // G: every x-edge and y-edge with at least one interior end. Summing
        // the edges leaving each interior node to the right and upwards, plus
        // the edges entering it from a clamped node on the left or below,
        // counts each such edge once.
        double sum = 0.0;
        for (const auto& s : spans)
        {
            const int base = s.row * nx;
            for (int i = s.begin; i < s.end; ++i)
            {
                const int l = base + i;
                const double right = w[l + 1] - w[l];
                const double up = w[l + nx] - w[l];
                sum += right * right + up * up;
                if (! mesh.interior[static_cast<std::size_t> (l - 1)])
                    sum += w[l] * w[l];
                if (! mesh.interior[static_cast<std::size_t> (l - nx)])
                    sum += w[l] * w[l];
            }
        }
        const double value = 0.25 * xi * sum * sum;
        const double scale = xi * sum;
        tension = scale;
        for (int l = 0; l < mesh.size(); ++l)
            gradient[l] = 0.0;
        for (const auto& s : spans)
        {
            const int base = s.row * nx;
            for (int i = s.begin; i < s.end; ++i)
            {
                const int l = base + i;
                gradient[l] = scale * (4.0 * w[l] - w[l - 1] - w[l + 1] - w[l - nx] - w[l + nx]);
            }
        }
        return value;
    }

    // The Berger tension is uniform, so its stress-stiffness is exactly
    // tension * K_1 and N* is the tension itself (StretchingProbe.h).
    double peakStretchingTension() const override { return tension; }
    bool boundedByEnergy() const override { return true; }
    std::string report() const override { return {}; }

private:
    struct Span
    {
        int row, begin, end;
    };
    MembraneGrid mesh;
    double xi = 0.0;
    double tension = 0.0; // xi int |grad w|^2 of the last evaluate, N/m
    std::vector<Span> spans;
};
} // namespace

std::unique_ptr<NonlinearPotential> makeBergerPotential()
{
    return std::make_unique<BergerPotential>();
}

} // namespace taikor::candidates
