// The dense generalised symmetric eigen-solver: residuals and M-orthonormality
// on a random pencil with clustered eigenvalues.

#include "ModalEigen.h"
#include "ModalTestUtil.h"

#include <algorithm>
#include <random>

using namespace taikor::candidates::modal;

MODAL_TEST (eigen_generalized_random)
{
    std::mt19937 rng (7);
    std::normal_distribution<double> g;
    const int n = 180;
    Matrix k (n), m (n), r (n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            r (i, j) = g (rng);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j)
        {
            double s = 0.0;
            for (int p = 0; p < n; ++p)
                s += r (i, p) * r (j, p);
            m (i, j) = m (j, i) = 0.01 * s / n + (i == j ? 1.0 + i : 0.0);
            const double kk = g (rng) * 1e3;
            k (i, j) = k (j, i) = (i == j ? 1e6 * (1 + (i / 3)) : kk); // clustered diagonal
        }
    std::vector<double> values;
    Matrix x;
    const bool ok = generalizedEigen (k, m, values, x);
    modaltest::check (ok, "converged", ok ? 1.0 : 0.0, 1.0, 0.0);
    double residual = 0.0, ortho = 0.0;
    for (int c = 0; c < n; ++c)
    {
        double scale = 0.0;
        for (int i = 0; i < n; ++i)
        {
            double kx = 0.0, mx = 0.0;
            for (int j = 0; j < n; ++j)
            {
                kx += k (i, j) * x (j, c);
                mx += m (i, j) * x (j, c);
            }
            residual = std::max (residual, std::abs (kx - values[static_cast<std::size_t> (c)] * mx));
            scale = std::max (scale, std::abs (kx));
        }
        residual /= std::max (scale, 1.0);
        for (int c2 = c; c2 < std::min (n, c + 4); ++c2)
        {
            double dot = 0.0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    dot += x (i, c) * m (i, j) * x (j, c2);
            ortho = std::max (ortho, std::abs (dot - (c == c2 ? 1.0 : 0.0)));
        }
    }
    modaltest::below ("max relative residual |Kx - lambda Mx|", residual, 1e-10);
    modaltest::below ("max |x^T M x - I| (near-diagonal band)", ortho, 1e-10);
    bool sorted = std::is_sorted (values.begin(), values.end());
    modaltest::check (sorted, "eigenvalues ascending", sorted, 1.0, 0.0);
}
