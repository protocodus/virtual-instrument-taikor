#pragma once

// Dense real symmetric eigen-solvers for the per-angular-order coupled blocks
// (at most a few hundred unknowns each): Householder tridiagonalisation with
// accumulated transformations followed by the implicit QL algorithm with
// Wilkinson shifts (the EISPACK tred2/tql2 pair, as in the public-domain JAMA
// EigenvalueDecomposition; Golub & Van Loan, Matrix Computations, 4th ed.,
// secs. 8.3.1 and 8.3.4). The generalised problem K x = lambda M x with M
// symmetric positive definite is reduced by the Cholesky factor M = L L^T.

#include <vector>

namespace taikor::candidates::modal
{

// Row-major square matrix.
struct Matrix
{
    int n = 0;
    std::vector<double> a;
    Matrix() = default;
    explicit Matrix (int size) : n (size), a (static_cast<std::size_t> (size) * static_cast<std::size_t> (size), 0.0) {}
    double& operator() (int i, int j) { return a[static_cast<std::size_t> (i) * static_cast<std::size_t> (n) + static_cast<std::size_t> (j)]; }
    double operator() (int i, int j) const { return a[static_cast<std::size_t> (i) * static_cast<std::size_t> (n) + static_cast<std::size_t> (j)]; }
};

// Eigenvalues ascending in values; vectors(i, k) is component i of vector k,
// orthonormal. Returns false if the QL iteration fails to converge.
bool symmetricEigen (const Matrix& symmetric, std::vector<double>& values, Matrix& vectors);

// K x = lambda M x. Eigenvectors are M-orthonormal (x_k^T M x_l = delta_kl).
// Returns false if M is not positive definite or the iteration fails.
bool generalizedEigen (const Matrix& stiffness, const Matrix& mass, std::vector<double>& values, Matrix& vectors);

} // namespace taikor::candidates::modal
