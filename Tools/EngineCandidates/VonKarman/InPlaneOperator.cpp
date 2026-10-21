// The in-plane stiffness operator K = D^T C D on the staggered grid (finest
// level, matrix-free), its Galerkin coarse versions (constant stencil plus
// rim corrections), and the grid-transfer operators between levels. See
// InPlaneSolver.h for the layout.

#include "InPlaneSolver.h"

#include <algorithm>
#include <cstring>
#include <type_traits>

namespace taikor::candidates::vk
{

namespace
{
// Prolongation weights of one coarse-to-fine transfer along one axis. The
// coarse cell I covers the fine cells 2I - 2 and 2I - 1 (the offset keeps
// every coarse index non-negative), so the coarse face I sits on the fine
// face 2I - 1.
//  - Normal axis of the component (x for u, y for v): a fine face on a coarse
//    face takes it with weight 1; a fine face inside a coarse cell averages
//    the cell's two faces (linear interpolation).
//  - Tangential axis: the fine row 2J - 2 lies 1/2 from the coarse row J and
//    3/2 from J - 1, so linear interpolation gives 3/4 and 1/4; the fine row
//    2J - 1 takes 3/4 of J and 1/4 of J + 1.
// This is the standard bilinear MAC-grid prolongation; with the Galerkin
// coarse operator P^T K P the choice affects only the convergence rate.
struct AxisWeights
{
    int index[2];
    double weight[2];
    int count;
};

inline AxisWeights normalWeights (int f) noexcept
{
    if ((f & 1) != 0)
        return { { (f + 1) / 2, 0 }, { 1.0, 0.0 }, 1 };
    return { { f / 2, f / 2 + 1 }, { 0.5, 0.5 }, 2 };
}

inline AxisWeights tangentialWeights (int f) noexcept
{
    if ((f & 1) != 0)
        return { { (f + 1) / 2, (f + 1) / 2 + 1 }, { 0.75, 0.25 }, 2 };
    return { { f / 2 + 1, f / 2 }, { 0.75, 0.25 }, 2 };
}
} // namespace

template <typename T>
void InPlaneSolver::applyT (int l, const T* x, T* y) const
{
    applyRows (l, x, y, 0, levels_[static_cast<std::size_t> (l)].n, 0);
}

template <typename T>
void InPlaneSolver::applyRows (int l, const T* x, T* y, int j0, int j1, int member) const
{
    const InPlaneLevel& L = levels_[static_cast<std::size_t> (l)];
    const Work& W = work_[static_cast<std::size_t> (l)];
    T* yu = y + L.base (0);
    T* yv = y + L.base (1);
    if (L.matrixFree)
    {
        // Energy (h^2/2) sum N : eps with eps = D x / h; its gradient in x is
        // h^2 (1/h) D^T C D x / h = D^T C D x, independent of h. One pass over
        // the rows with rolling one-row stress buffers (in L1):
        //   y_u = sxx[k] - sxx[k+1] + sxy[k-n] - sxy[k]   (sxx row j, sxy rows j-1, j)
        //   y_v = syy[k] - syy[k+n] + sxy[k-1] - sxy[k]   (syy rows j, j+1, sxy row j)
        // with sxx, syy at the nodes (u row j, v rows j-1, j) and sxy at the
        // corners (u rows j, j+1, v row j).
        const int n = L.n;
        T* buffer;
        if constexpr (std::is_same_v<T, float>)
            buffer = W.sxxF.data() + static_cast<std::size_t> (member) * 6 * static_cast<std::size_t> (n);
        else
            buffer = W.sxx.data() + static_cast<std::size_t> (member) * 6 * static_cast<std::size_t> (n);
        T* nodeXX[2] = { buffer, buffer + n };
        T* nodeYY[2] = { buffer + 2 * n, buffer + 3 * n };
        T* corner[2] = { buffer + 4 * n, buffer + 5 * n };
        const T* u = x + L.base (0);
        const T* v = x + L.base (1);
        const T c11 = static_cast<T> (c_.c11), c12 = static_cast<T> (c_.c12), c66 = static_cast<T> (c_.c66);
        auto nodeRow = [&] (int j, T* sxx, T* syy)
        {
            const T* __restrict ur = u + j * n;
            const T* __restrict vr = v + j * n;
            for (int i = 0; i < n; ++i)
            {
                const T exx = ur[i] - ur[i - 1];
                const T eyy = vr[i] - vr[i - n];
                sxx[i] = c11 * exx + c12 * eyy;
                syy[i] = c12 * exx + c11 * eyy;
            }
        };
        auto cornerRow = [&] (int j, T* sxy)
        {
            const T* __restrict ur = u + j * n;
            const T* __restrict vr = v + j * n;
            for (int i = 0; i < n; ++i)
                sxy[i] = c66 * (ur[i + n] - ur[i] + vr[i + 1] - vr[i]);
        };
        // Row -1 of the corner array and rows < 0 of the nodes read the pads.
        // A row range starts from its first row's neighbours the same way.
        cornerRow (j0 - 1, corner[(j0 & 1) ^ 1]);
        nodeRow (j0, nodeXX[j0 & 1], nodeYY[j0 & 1]);
        const std::uint8_t* fu = L.freeU.data();
        const std::uint8_t* fv = L.freeV.data();
        for (int j = j0; j < j1; ++j)
        {
            const int cur = j & 1, nxt = cur ^ 1;
            cornerRow (j, corner[cur]);
            if (j + 1 < n)
                nodeRow (j + 1, nodeXX[nxt], nodeYY[nxt]);
            else
                std::fill (nodeYY[nxt], nodeYY[nxt] + n, T (0));
            const T* __restrict sxx = nodeXX[cur];
            const T* __restrict syy = nodeYY[cur];
            const T* __restrict syyUp = nodeYY[nxt];
            const T* __restrict sxy = corner[cur];
            const T* __restrict sxyDown = corner[nxt]; // row j - 1
            const std::uint8_t* mu = fu + j * n;
            const std::uint8_t* mv = fv + j * n;
            T* __restrict yur = yu + j * n;
            T* __restrict yvr = yv + j * n;
            // i = n - 1 (the last column) is always clamped: no free dof reads i + 1.
            for (int i = 0; i + 1 < n; ++i)
            {
                const T gu = sxx[i] - sxx[i + 1] + sxyDown[i] - sxy[i];
                const T gv = syy[i] - syyUp[i] + (i > 0 ? sxy[i - 1] : T (0)) - sxy[i];
                yur[i] = mu[i] != 0 ? gu : T (0);
                yvr[i] = mv[i] != 0 ? gv : T (0);
            }
            yur[n - 1] = T (0);
            yvr[n - 1] = T (0);
        }
        return;
    }
    // Constant stencil on every row: four 16-byte vectors of rows at a time
    // with the accumulators in registers, the stencil's taps streaming x.
    typedef T Vec __attribute__ ((vector_size (16)));
    constexpr int lanes = static_cast<int> (16 / sizeof (T));
    constexpr int block = 4 * lanes;
    auto load = [] (const T* p) { Vec v; std::memcpy (&v, p, sizeof (Vec)); return v; };
    const int kBegin = j0 * L.n, kEnd = j1 * L.n;
    for (int c = 0; c < 2; ++c)
    {
        const std::uint8_t* mask = c == 0 ? L.freeU.data() : L.freeV.data();
        const T* xs = x + L.base (c);
        T* ys = y + L.base (c);
        const auto& taps = L.stencil[c];
        const std::size_t count = std::min<std::size_t> (taps.size(), 64);
        T coef[64];
        long off[64];
        for (std::size_t t = 0; t < count; ++t)
        {
            if constexpr (std::is_same_v<T, float>)
                coef[t] = L.stencilF[c][t];
            else
                coef[t] = taps[t].coefficient;
            off[t] = taps[t].offset;
        }
        int k0 = kBegin;
        for (; k0 + block <= kEnd; k0 += block)
        {
            Vec a0 {}, a1 {}, a2 {}, a3 {};
            for (std::size_t t = 0; t < count; ++t)
            {
                const Vec ct = Vec {} + coef[t];
                const T* xt = xs + k0 + off[t];
                a0 += ct * load (xt);
                a1 += ct * load (xt + lanes);
                a2 += ct * load (xt + 2 * lanes);
                a3 += ct * load (xt + 3 * lanes);
            }
            T out[block];
            std::memcpy (out, &a0, sizeof (Vec));
            std::memcpy (out + lanes, &a1, sizeof (Vec));
            std::memcpy (out + 2 * lanes, &a2, sizeof (Vec));
            std::memcpy (out + 3 * lanes, &a3, sizeof (Vec));
            for (int k = 0; k < block; ++k)
                ys[k0 + k] = mask[k0 + k] != 0 ? out[k] : T (0);
        }
        for (int k = k0; k < kEnd; ++k)
        {
            T a = T (0);
            for (std::size_t t = 0; t < count; ++t)
                a += coef[t] * xs[k + off[t]];
            ys[k] = mask[k] != 0 ? a : T (0);
        }
    }
    // Rows near the rim: add their differences from the stencil, slot by
    // slot (independent gathers, no dependency chain along a row).
    // The rim rows are listed in ascending order (u block, then v block), so
    // a row range's are two runs of the list.
    const std::size_t rows = L.irregularRow.size();
    if (rows == 0)
        return;
    const Work& scratch = work_[static_cast<std::size_t> (l)];
    T* acc;
    if constexpr (std::is_same_v<T, float>)
        acc = scratch.rimF.data();
    else
        acc = scratch.rim.data();
    const long* rowIndex = L.irregularRow.data();
    for (int c = 0; c < 2; ++c)
    {
        const auto r0 = static_cast<std::size_t> (std::lower_bound (rowIndex, rowIndex + rows, L.base (c) + kBegin) - rowIndex);
        const auto r1 = static_cast<std::size_t> (std::lower_bound (rowIndex, rowIndex + rows, L.base (c) + kEnd) - rowIndex);
        for (std::size_t r = r0; r < r1; ++r)
            acc[r] = T (0);
        for (int slot = 0; slot < L.irregularWidth; ++slot)
        {
            const long* off = L.irregularOffset.data() + static_cast<std::size_t> (slot) * rows;
            const T* coef;
            if constexpr (std::is_same_v<T, float>)
                coef = L.irregularCoefficientF.data() + static_cast<std::size_t> (slot) * rows;
            else
                coef = L.irregularCoefficient.data() + static_cast<std::size_t> (slot) * rows;
            for (std::size_t r = r0; r < r1; ++r)
                acc[r] += coef[r] * x[rowIndex[r] + off[r]];
        }
        for (std::size_t r = r0; r < r1; ++r)
            y[rowIndex[r]] += acc[r];
    }
}

void InPlaneSolver::apply (int l, const double* x, double* y) const { applyT (l, x, y); }
void InPlaneSolver::apply (int l, const float* x, float* y) const { applyT (l, x, y); }

template <typename T>
void InPlaneSolver::prolongT (int coarse, const T* xc, T* xf) const
{
    prolongRows (coarse, xc, xf, 0, levels_[static_cast<std::size_t> (coarse - 1)].n);
}

template <typename T>
void InPlaneSolver::prolongRows (int coarse, const T* xc, T* xf, int j0, int j1) const
{
    // xf = mask_f (P xc), row by row: each fine row blends two coarse rows
    // with the tangential (u) or normal (v) weights into R, then maps R along
    // the row in pairs of fine faces i = 2m - 1 (odd) and 2m (even):
    //   u (normal axis x): odd on coarse face m, even between faces m and m+1;
    //   v (tangential x):  odd 3/4 R[m] + 1/4 R[m+1], even 3/4 R[m+1] + 1/4 R[m].
    // Coarse vectors are zero off their free dofs; i = 0 is always clamped.
    const InPlaneLevel& F = levels_[static_cast<std::size_t> (coarse - 1)];
    const InPlaneLevel& C = levels_[static_cast<std::size_t> (coarse)];
    const int nf = F.n, nc = C.n;
    T R[1024];
    for (int comp = 0; comp < 2; ++comp)
    {
        const std::uint8_t* freeF = comp == 0 ? F.freeU.data() : F.freeV.data();
        const T* src = xc + C.base (comp);
        T* dst = xf + F.base (comp);
        for (int j = j0; j < j1; ++j)
        {
            const AxisWeights ay = comp == 0 ? tangentialWeights (j) : normalWeights (j);
            const T* r0 = src + ay.index[0] * nc;
            const T* r1 = src + ay.index[ay.count - 1] * nc;
            const T w0 = static_cast<T> (ay.weight[0]), w1 = static_cast<T> (ay.count > 1 ? ay.weight[1] : 0.0);
            for (int m = 0; m < nc; ++m)
                R[m] = w0 * r0[m] + w1 * r1[m];
            const std::uint8_t* mask = freeF + j * nf;
            T* out = dst + j * nf;
            out[0] = T (0);
            const int pairs = (nf - 1) / 2;
            if (comp == 0)
                for (int m = 1; m <= pairs; ++m)
                {
                    const T odd = R[m];
                    const T even = T (0.5) * (R[m] + R[m + 1]);
                    out[2 * m - 1] = mask[2 * m - 1] != 0 ? odd : T (0);
                    out[2 * m] = mask[2 * m] != 0 ? even : T (0);
                }
            else
                for (int m = 1; m <= pairs; ++m)
                {
                    const T odd = T (0.75) * R[m] + T (0.25) * R[m + 1];
                    const T even = T (0.75) * R[m + 1] + T (0.25) * R[m];
                    out[2 * m - 1] = mask[2 * m - 1] != 0 ? odd : T (0);
                    out[2 * m] = mask[2 * m] != 0 ? even : T (0);
                }
            if ((nf - 1) % 2 != 0)
                out[nf - 1] = T (0); // nf even: the last face is clamped
        }
    }
}

template <typename T>
void InPlaneSolver::restrictT (int coarse, const T* rf, T* rc) const
{
    restrictRows (coarse, rf, rc, 0, levels_[static_cast<std::size_t> (coarse)].n);
}

template <typename T>
void InPlaneSolver::restrictRows (int coarse, const T* rf, T* rc, int j0, int j1) const
{
    // rc = P^T rf as a gather (fine vectors are zero off their free dofs):
    //   u(I, J) = sum_{a = -2..0} (1/2, 1, 1/2)_a sum_{b = -3..0} (1/4, 3/4, 3/4, 1/4)_b u_f(2I + a, 2J + b),
    //   v(I, J) = the same with the axes exchanged,
    // the transpose of prolong() entry by entry (tests/vonkarman checks
    // <P x, y> = <x, P^T y>). Rows -3..n of the fine block lie in its zero
    // pads; non-free coarse dofs are masked.
    const InPlaneLevel& F = levels_[static_cast<std::size_t> (coarse - 1)];
    const InPlaneLevel& C = levels_[static_cast<std::size_t> (coarse)];
    const int nf = F.n, nc = C.n;
    constexpr T wn[3] = { T (0.5), T (1.0), T (0.5) };
    constexpr T wt[4] = { T (0.25), T (0.75), T (0.75), T (0.25) };
    for (int comp = 0; comp < 2; ++comp)
    {
        const std::uint8_t* freeC = comp == 0 ? C.freeU.data() : C.freeV.data();
        const T* src = rf + F.base (comp);
        T* dst = rc + C.base (comp);
        for (int J = j0; J < j1; ++J)
            for (int I = 0; I < nc; ++I)
            {
                const int k = J * nc + I;
                if (freeC[k] == 0)
                {
                    dst[k] = T (0);
                    continue;
                }
                T s = T (0);
                if (comp == 0)
                {
                    for (int b = 0; b < 4; ++b)
                    {
                        const T* row = src + (2 * J - 3 + b) * nf + 2 * I - 2;
                        s += wt[b] * (wn[0] * row[0] + wn[1] * row[1] + wn[2] * row[2]);
                    }
                }
                else
                {
                    for (int a = 0; a < 3; ++a)
                    {
                        const T* row = src + (2 * J - 2 + a) * nf + 2 * I - 3;
                        s += wn[a] * (wt[0] * row[0] + wt[1] * row[1] + wt[2] * row[2] + wt[3] * row[3]);
                    }
                }
                dst[k] = s;
            }
    }
}

template void InPlaneSolver::applyRows (int, const float*, float*, int, int, int) const;
template void InPlaneSolver::applyRows (int, const double*, double*, int, int, int) const;
template void InPlaneSolver::prolongRows (int, const float*, float*, int, int) const;
template void InPlaneSolver::restrictRows (int, const float*, float*, int, int) const;

void InPlaneSolver::prolong (int coarse, const double* xc, double* xf) const { prolongT (coarse, xc, xf); }
void InPlaneSolver::restrict_ (int coarse, const double* rf, double* rc) const { restrictT (coarse, rf, rc); }
void InPlaneSolver::prolong (int coarse, const float* xc, float* xf) const { prolongT (coarse, xc, xf); }
void InPlaneSolver::restrict_ (int coarse, const float* rf, float* rc) const { restrictT (coarse, rf, rc); }

// Which coarse dofs are free: those the prolongation maps onto at least one
// free fine dof, so the coarse space spans exactly the masked P.
void markCoarseFree (const InPlaneLevel& F, InPlaneLevel& C)
{
    C.freeU.assign (static_cast<std::size_t> (C.N), 0);
    C.freeV.assign (static_cast<std::size_t> (C.N), 0);
    for (int comp = 0; comp < 2; ++comp)
    {
        const std::uint8_t* freeF = comp == 0 ? F.freeU.data() : F.freeV.data();
        std::uint8_t* freeC = comp == 0 ? C.freeU.data() : C.freeV.data();
        for (int j = 0; j < F.n; ++j)
            for (int i = 0; i < F.n; ++i)
            {
                if (freeF[j * F.n + i] == 0)
                    continue;
                const AxisWeights ax = comp == 0 ? normalWeights (i) : tangentialWeights (i);
                const AxisWeights ay = comp == 0 ? tangentialWeights (j) : normalWeights (j);
                for (int b = 0; b < ay.count; ++b)
                    for (int a = 0; a < ax.count; ++a)
                        freeC[ay.index[b] * C.n + ax.index[a]] = 1;
            }
    }
}

} // namespace taikor::candidates::vk
