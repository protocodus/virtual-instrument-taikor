#include "ModalEigen.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace taikor::candidates::modal
{
namespace
{
// Householder reduction of the symmetric matrix held in w to tridiagonal form
// d (diagonal) and e (sub-diagonal, e[0] = 0), with the orthogonal
// transformation V accumulated as its transpose in w (row k = column k of V).
void tridiagonalize (Matrix& w, std::vector<double>& d, std::vector<double>& e)
{
    // Works on the transpose w = V^T so the inner loops run along rows.
    const int n = w.n;
    const auto v = [&w, n] (int r, int c) -> double& { return w.a[static_cast<std::size_t> (c) * static_cast<std::size_t> (n) + static_cast<std::size_t> (r)]; };
    for (int j = 0; j < n; ++j)
        d[static_cast<std::size_t> (j)] = v (n - 1, j);
    for (int i = n - 1; i > 0; --i)
    {
        double scale = 0.0, h = 0.0;
        for (int k = 0; k < i; ++k)
            scale += std::abs (d[static_cast<std::size_t> (k)]);
        if (scale == 0.0)
        {
            e[static_cast<std::size_t> (i)] = d[static_cast<std::size_t> (i - 1)];
            for (int j = 0; j < i; ++j)
            {
                d[static_cast<std::size_t> (j)] = v (i - 1, j);
                v (i, j) = 0.0;
                v (j, i) = 0.0;
            }
        }
        else
        {
            for (int k = 0; k < i; ++k)
            {
                d[static_cast<std::size_t> (k)] /= scale;
                h += d[static_cast<std::size_t> (k)] * d[static_cast<std::size_t> (k)];
            }
            double f = d[static_cast<std::size_t> (i - 1)];
            double g = std::sqrt (h);
            if (f > 0.0)
                g = -g;
            e[static_cast<std::size_t> (i)] = scale * g;
            h -= f * g;
            d[static_cast<std::size_t> (i - 1)] = f - g;
            for (int j = 0; j < i; ++j)
                e[static_cast<std::size_t> (j)] = 0.0;
            for (int j = 0; j < i; ++j)
            {
                f = d[static_cast<std::size_t> (j)];
                v (j, i) = f;
                g = e[static_cast<std::size_t> (j)] + v (j, j) * f;
                for (int k = j + 1; k <= i - 1; ++k)
                {
                    g += v (k, j) * d[static_cast<std::size_t> (k)];
                    e[static_cast<std::size_t> (k)] += v (k, j) * f;
                }
                e[static_cast<std::size_t> (j)] = g;
            }
            f = 0.0;
            for (int j = 0; j < i; ++j)
            {
                e[static_cast<std::size_t> (j)] /= h;
                f += e[static_cast<std::size_t> (j)] * d[static_cast<std::size_t> (j)];
            }
            const double hh = f / (h + h);
            for (int j = 0; j < i; ++j)
                e[static_cast<std::size_t> (j)] -= hh * d[static_cast<std::size_t> (j)];
            for (int j = 0; j < i; ++j)
            {
                f = d[static_cast<std::size_t> (j)];
                g = e[static_cast<std::size_t> (j)];
                for (int k = j; k <= i - 1; ++k)
                    v (k, j) -= (f * e[static_cast<std::size_t> (k)] + g * d[static_cast<std::size_t> (k)]);
                d[static_cast<std::size_t> (j)] = v (i - 1, j);
                v (i, j) = 0.0;
            }
        }
        d[static_cast<std::size_t> (i)] = h;
    }
    for (int i = 0; i < n - 1; ++i)
    {
        v (n - 1, i) = v (i, i);
        v (i, i) = 1.0;
        const double h = d[static_cast<std::size_t> (i + 1)];
        if (h != 0.0)
        {
            for (int k = 0; k <= i; ++k)
                d[static_cast<std::size_t> (k)] = v (k, i + 1) / h;
            for (int j = 0; j <= i; ++j)
            {
                double g = 0.0;
                for (int k = 0; k <= i; ++k)
                    g += v (k, i + 1) * v (k, j);
                for (int k = 0; k <= i; ++k)
                    v (k, j) -= g * d[static_cast<std::size_t> (k)];
            }
        }
        for (int k = 0; k <= i; ++k)
            v (k, i + 1) = 0.0;
    }
    for (int j = 0; j < n; ++j)
    {
        d[static_cast<std::size_t> (j)] = v (n - 1, j);
        v (n - 1, j) = 0.0;
    }
    v (n - 1, n - 1) = 1.0;
    e[0] = 0.0;
}

// Implicit QL on the tridiagonal (d, e). vt holds the transformation
// transposed (row k = column k of the eigenvector matrix) so the plane
// rotations touch two contiguous rows.
bool implicitQL (Matrix& vt, std::vector<double>& d, std::vector<double>& e)
{
    const int n = vt.n;
    for (int i = 1; i < n; ++i)
        e[static_cast<std::size_t> (i - 1)] = e[static_cast<std::size_t> (i)];
    e[static_cast<std::size_t> (n - 1)] = 0.0;
    double f = 0.0, tst1 = 0.0;
    const double eps = std::ldexp (1.0, -52);
    for (int l = 0; l < n; ++l)
    {
        const auto L = static_cast<std::size_t> (l);
        tst1 = std::max (tst1, std::abs (d[L]) + std::abs (e[L]));
        int m = l;
        while (m < n)
        {
            if (std::abs (e[static_cast<std::size_t> (m)]) <= eps * tst1)
                break;
            ++m;
        }
        if (m > l)
        {
            int iteration = 0;
            do
            {
                if (++iteration > 60)
                    return false;
                double g = d[L];
                double p = (d[L + 1] - g) / (2.0 * e[L]);
                double r = std::hypot (p, 1.0);
                if (p < 0.0)
                    r = -r;
                d[L] = e[L] / (p + r);
                d[L + 1] = e[L] * (p + r);
                const double dl1 = d[L + 1];
                double h = g - d[L];
                for (int i = l + 2; i < n; ++i)
                    d[static_cast<std::size_t> (i)] -= h;
                f += h;
                p = d[static_cast<std::size_t> (m)];
                double c = 1.0, c2 = c, c3 = c;
                const double el1 = e[L + 1];
                double s = 0.0, s2 = 0.0;
                for (int i = m - 1; i >= l; --i)
                {
                    const auto I = static_cast<std::size_t> (i);
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[I];
                    h = c * p;
                    r = std::hypot (p, e[I]);
                    e[I + 1] = s * r;
                    s = e[I] / r;
                    c = p / r;
                    p = c * d[I] - s * g;
                    d[I + 1] = h + s * (c * g + s * d[I]);
                    double* rowI = &vt.a[I * static_cast<std::size_t> (n)];
                    double* rowJ = rowI + n;
                    for (int k = 0; k < n; ++k)
                    {
                        const double hk = rowJ[k];
                        rowJ[k] = s * rowI[k] + c * hk;
                        rowI[k] = c * rowI[k] - s * hk;
                    }
                }
                p = -s * s2 * c3 * el1 * e[L] / dl1;
                e[L] = s * p;
                d[L] = c * p;
            } while (std::abs (e[L]) > eps * tst1);
        }
        d[L] += f;
        e[L] = 0.0;
    }
    return true;
}
} // namespace

bool symmetricEigen (const Matrix& symmetric, std::vector<double>& values, Matrix& vectors)
{
    const int n = symmetric.n;
    values.assign (static_cast<std::size_t> (n), 0.0);
    vectors = Matrix (n);
    if (n == 0)
        return true;
    if (n == 1)
    {
        values[0] = symmetric (0, 0);
        vectors (0, 0) = 1.0;
        return true;
    }
    Matrix vt = symmetric; // symmetric, so it is its own transpose
    std::vector<double> e (static_cast<std::size_t> (n), 0.0);
    tridiagonalize (vt, values, e);
    if (! implicitQL (vt, values, e))
        return false;
    std::vector<int> order (static_cast<std::size_t> (n));
    std::iota (order.begin(), order.end(), 0);
    std::sort (order.begin(), order.end(),
               [&] (int x, int y) { return values[static_cast<std::size_t> (x)] < values[static_cast<std::size_t> (y)]; });
    std::vector<double> sorted (static_cast<std::size_t> (n));
    for (int k = 0; k < n; ++k)
    {
        const int from = order[static_cast<std::size_t> (k)];
        sorted[static_cast<std::size_t> (k)] = values[static_cast<std::size_t> (from)];
        for (int i = 0; i < n; ++i)
            vectors (i, k) = vt (from, i);
    }
    values = std::move (sorted);
    return true;
}

bool generalizedEigen (const Matrix& stiffness, const Matrix& mass, std::vector<double>& values, Matrix& vectors)
{
    const int n = mass.n;
    // Cholesky M = L L^T (lower triangle of l).
    Matrix l (n);
    for (int j = 0; j < n; ++j)
    {
        double diag = mass (j, j);
        for (int k = 0; k < j; ++k)
            diag -= l (j, k) * l (j, k);
        if (! (diag > 0.0))
            return false;
        const double root = std::sqrt (diag);
        l (j, j) = root;
        for (int i = j + 1; i < n; ++i)
        {
            double sum = mass (i, j);
            for (int k = 0; k < j; ++k)
                sum -= l (i, k) * l (j, k);
            l (i, j) = sum / root;
        }
    }
    // X = L^{-1} K row by row (X(i,:) = (K(i,:) - sum_k L(i,k) X(k,:)) / L(i,i)),
    // then A = L^{-1} X^T = L^{-1} K L^{-T} the same way; rows are contiguous.
    const auto forwardRows = [&] (Matrix& x)
    {
        for (int i = 0; i < n; ++i)
        {
            double* row = &x.a[static_cast<std::size_t> (i) * static_cast<std::size_t> (n)];
            for (int k = 0; k < i; ++k)
            {
                const double lik = l (i, k);
                const double* rk = &x.a[static_cast<std::size_t> (k) * static_cast<std::size_t> (n)];
                for (int c = 0; c < n; ++c)
                    row[c] -= lik * rk[c];
            }
            const double inv = 1.0 / l (i, i);
            for (int c = 0; c < n; ++c)
                row[c] *= inv;
        }
    };
    Matrix b = stiffness;
    forwardRows (b);
    Matrix a (n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            a (i, j) = b (j, i);
    forwardRows (a);
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
        {
            const double mean = 0.5 * (a (i, j) + a (j, i));
            a (i, j) = mean;
            a (j, i) = mean;
        }
    Matrix y;
    if (! symmetricEigen (a, values, y))
        return false;
    // X = L^{-T} Y: solve L^T X = Y from the last row up, row by row.
    vectors = y;
    for (int i = n - 1; i >= 0; --i)
    {
        double* row = &vectors.a[static_cast<std::size_t> (i) * static_cast<std::size_t> (n)];
        for (int k = i + 1; k < n; ++k)
        {
            const double lki = l (k, i);
            const double* rk = &vectors.a[static_cast<std::size_t> (k) * static_cast<std::size_t> (n)];
            for (int c = 0; c < n; ++c)
                row[c] -= lki * rk[c];
        }
        const double inv = 1.0 / l (i, i);
        for (int c = 0; c < n; ++c)
            row[c] *= inv;
    }
    return true;
}

} // namespace taikor::candidates::modal
