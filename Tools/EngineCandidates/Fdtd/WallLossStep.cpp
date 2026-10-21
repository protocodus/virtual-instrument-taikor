// The barrel wall's viscothermal boundary layer: one step and its energy
// (see WallLoss.h).

#include "WallLoss.h"

#include <cmath>

namespace taikor::candidates::fdtd
{

template <typename Real>
double WallLoss<Real>::faceDrag (int f, long long n) const noexcept
{
    const double* d = drag[static_cast<std::size_t> (n & 1)].data();
    double sum = 0.0;
    for (int q = faceStart[static_cast<std::size_t> (f)]; q < faceStart[static_cast<std::size_t> (f) + 1]; ++q)
        sum += faceWeights[static_cast<std::size_t> (q)] * d[faceElements[static_cast<std::size_t> (q)]];
    return sum;
}

template <typename Real>
double WallLoss<Real>::correction (const Cell& e, long long n) const noexcept
{
    // The cell's faces carry G = g - drag: the - side cell loses lambda^2 drag,
    // the + side cell gains it.
    double sum = 0.0;
    for (int q = e.faceBegin; q < e.faceEnd; ++q)
    {
        const int code = cellFaces[static_cast<std::size_t> (q)];
        const double d = faceDrag (code >> 1, n);
        sum += (code & 1) ? d : -d;
    }
    return lambda2 * sum;
}

template <typename Real>
double WallLoss<Real>::thermalUpdate (int t, long long n, double explicitValue, double previous) const noexcept
{
    // (1 + beta) Psi^{n+1} = explicit + beta Psi^{n-1}
    //                        + (c^2 k^2 / Delta^3) A C_th sum_i mu_i xi_i theta_i^{n-1} / (1 + k xi_i).
    const double* old = theta[static_cast<std::size_t> ((n + 1) & 1)].data() + static_cast<std::size_t> (t) * xi.size();
    double memory = 0.0;
    for (std::size_t i = 0; i < xi.size(); ++i)
        memory += muThermal[i] * xi[i] * relief[i] * old[i];
    const double b = c * c * k * k / (delta * delta * delta) * thermalScale * cTh * cellArea[static_cast<std::size_t> (t)] * memory;
    const double bt = beta[static_cast<std::size_t> (t)];
    return (explicitValue + bt * previous + b) / (1.0 + bt);
}

template <typename Real>
double WallLoss<Real>::thermalBranches (int t, long long n, double next, double previous) noexcept
{
    // theta_i^{n+1} = (2k p + (1 - k xi_i) theta_i^{n-1}) / (1 + k xi_i),
    // p = rho (Psi^{n+1} - Psi^{n-1}) / 2k, written over theta^{n-1}.
    double* slot = theta[static_cast<std::size_t> ((n + 1) & 1)].data() + static_cast<std::size_t> (t) * xi.size();
    const double p = rho * (next - previous) / (2.0 * k);
    double lostHere = 0.0;
    for (std::size_t i = 0; i < xi.size(); ++i)
    {
        const double old = slot[i];
        slot[i] = gain[i] * p + decay[i] * old;
        if (tracking)
        {
            const double rate = (slot[i] - old) / (2.0 * k);
            lostHere += muThermal[i] * rate * rate;
        }
    }
    return tracking ? k * thermalScale * cTh * cellArea[static_cast<std::size_t> (t)] * lostHere : 0.0;
}

template <typename Real>
double WallLoss<Real>::stepFriction (long long n, int begin, int end) noexcept
{
    const Real* cur = domain->level (n);
    const std::size_t poles = xi.size();
    double* dragNext = drag[static_cast<std::size_t> ((n + 1) & 1)].data();
    double* slots = phi[static_cast<std::size_t> ((n + 1) & 1)].data();
    double lostHere = 0.0;
    for (int e = begin; e < end; ++e)
    {
        const auto& element = elements[static_cast<std::size_t> (e)];
        // V_e^n = sum_f w_ef G_f^n, G_f = g_f - sum_e' w_e'f D_e'.
        double sensed = 0.0;
        for (int q = element.begin; q < element.end; ++q)
        {
            const int f = elementFaces[static_cast<std::size_t> (q)];
            const auto& face = faces[static_cast<std::size_t> (f)];
            const double flux = static_cast<double> (cur[face.high]) - cur[face.low] - faceDrag (f, n);
            sensed += elementWeights[static_cast<std::size_t> (q)] * flux;
        }
        double* slot = slots + static_cast<std::size_t> (e) * poles;
        double sum = 0.0, lostElement = 0.0;
        for (std::size_t i = 0; i < poles; ++i)
        {
            const double old = slot[i];
            slot[i] = gain[i] * sensed + decay[i] * old;
            sum += muViscous[i] * slot[i];
            if (tracking)
            {
                const double rate = (slot[i] - old) / (2.0 * k);
                lostElement += muViscous[i] * rate * rate;
            }
        }
        dragNext[e] = faceGain[static_cast<std::size_t> (e)] * sum;
        lostHere += faceGain[static_cast<std::size_t> (e)] * lostElement;
    }
    return tracking ? k * rho * delta * lostHere : 0.0;
}

template <typename Real>
double WallLoss<Real>::finishCells (long long n, int begin, int end) noexcept
{
    Real* next = const_cast<Real*> (domain->level (n + 1));
    const Real* prev = domain->level (n - 1);
    double lostHere = 0.0;
    for (int q = begin; q < end; ++q)
    {
        const auto& e = cells[static_cast<std::size_t> (q)];
        const double value = static_cast<double> (next[e.cell]) + correction (e, n);
        if (e.thermal >= 0)
        {
            next[e.cell] = static_cast<Real> (thermalUpdate (e.thermal, n, value, prev[e.cell]));
            lostHere += thermalBranches (e.thermal, n, next[e.cell], prev[e.cell]);
        }
        else
            next[e.cell] = static_cast<Real> (value);
    }
    return lostHere;
}

template <typename Real>
void WallLoss<Real>::slicePartial (long long n, int entry) noexcept
{
    Real* next = const_cast<Real*> (domain->level (n + 1));
    const Real* prev = domain->level (n - 1);
    const auto& e = cells[static_cast<std::size_t> (entry)];
    double value = static_cast<double> (next[e.cell]) + correction (e, n);
    if (e.thermal >= 0)
        value = thermalUpdate (e.thermal, n, value, prev[e.cell]);
    next[e.cell] = static_cast<Real> (value);
}

template <typename Real>
double WallLoss<Real>::sliceFactor (int entry) const noexcept
{
    const int t = cells[static_cast<std::size_t> (entry)].thermal;
    return t >= 0 ? 1.0 / (1.0 + beta[static_cast<std::size_t> (t)]) : 1.0;
}

template <typename Real>
double WallLoss<Real>::sliceFinish (long long n, int entry) noexcept
{
    const auto& e = cells[static_cast<std::size_t> (entry)];
    if (e.thermal < 0)
        return 0.0;
    return thermalBranches (e.thermal, n, domain->level (n + 1)[e.cell], domain->level (n - 1)[e.cell]);
}

template <typename Real>
double WallLoss<Real>::energy (long long n) const
{
    const Real* next = domain->level (n + 1);
    const Real* now = domain->level (n);
    const std::size_t poles = xi.size();
    const double* p1 = phi[static_cast<std::size_t> ((n + 1) & 1)].data();
    const double* p0 = phi[static_cast<std::size_t> (n & 1)].data();
    double faceTerm = 0.0, viscousStore = 0.0, thermalStore = 0.0;
    for (std::size_t f = 0; f < faces.size(); ++f)
    {
        const double g1 = static_cast<double> (next[faces[f].high]) - next[faces[f].low];
        const double g0 = static_cast<double> (now[faces[f].high]) - now[faces[f].low];
        faceTerm += (g1 - faceDrag (static_cast<int> (f), n + 1)) * (g0 - faceDrag (static_cast<int> (f), n)) - g1 * g0;
    }
    for (std::size_t e = 0; e < elements.size(); ++e)
    {
        double store = 0.0;
        for (std::size_t i = 0; i < poles; ++i)
            store += muViscous[i] * xi[i] * (p1[e * poles + i] * p1[e * poles + i] + p0[e * poles + i] * p0[e * poles + i]);
        viscousStore += faceGain[e] * store;
    }
    const double* t1 = theta[static_cast<std::size_t> ((n + 1) & 1)].data();
    const double* t0 = theta[static_cast<std::size_t> (n & 1)].data();
    for (std::size_t t = 0; t < cellArea.size(); ++t)
    {
        double store = 0.0;
        for (std::size_t i = 0; i < poles; ++i)
            store += muThermal[i] * xi[i] * (t1[t * poles + i] * t1[t * poles + i] + t0[t * poles + i] * t0[t * poles + i]);
        thermalStore += cellArea[t] * store;
    }
    return 0.5 * rho * delta * faceTerm + 0.25 * rho * delta * viscousStore + 0.25 * thermalScale * cTh * thermalStore;
}

#define TAIKOR_WALL_INSTANTIATE(Real)                                                       \
    template double WallLoss<Real>::faceDrag (int, long long) const noexcept;              \
    template double WallLoss<Real>::correction (const Cell&, long long) const noexcept;    \
    template double WallLoss<Real>::thermalUpdate (int, long long, double, double) const noexcept; \
    template double WallLoss<Real>::thermalBranches (int, long long, double, double) noexcept; \
    template double WallLoss<Real>::stepFriction (long long, int, int) noexcept;           \
    template double WallLoss<Real>::finishCells (long long, int, int) noexcept;            \
    template void WallLoss<Real>::slicePartial (long long, int) noexcept;                  \
    template double WallLoss<Real>::sliceFactor (int) const noexcept;                      \
    template double WallLoss<Real>::sliceFinish (long long, int) noexcept;                 \
    template double WallLoss<Real>::energy (long long) const;
TAIKOR_WALL_INSTANTIATE (float)
TAIKOR_WALL_INSTANTIATE (double)
#undef TAIKOR_WALL_INSTANTIATE

} // namespace taikor::candidates::fdtd
