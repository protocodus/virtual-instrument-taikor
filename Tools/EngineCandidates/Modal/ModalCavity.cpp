#include "ModalCavity.h"

#include <algorithm>
#include <cmath>

namespace taikor::candidates::modal
{
namespace
{
void jAndPrime (int m, double x, double& j, double& jp)
{
    thread_local std::vector<double> seq;
    special::besselJSequence (x, m + 1, seq);
    j = seq[static_cast<std::size_t> (m)];
    jp = 0.5 * ((m == 0 ? -seq[1] : seq[static_cast<std::size_t> (m - 1)]) - seq[static_cast<std::size_t> (m + 1)]);
}

// The shared decision's air constants for the Kirchhoff side-wall loss (the
// same in every candidate): air at 20 degC (the temperature of the spec's
// rho = 1.2041 kg/m^3), dynamic viscosity 1.81e-5 Pa s, Prandtl number 0.71,
// ratio of specific heats 1.40 (Kinsler, Frey, Coppens & Sanders, Fundamentals
// of Acoustics, 4th ed., Wiley 2000, Appendix A10; Pierce, Acoustics, ch. 10).
// Density and sound speed come from DrumSpec through AirSpec.
constexpr double airViscosity = 1.81e-5;
constexpr double airPrandtl = 0.71;
constexpr double airGamma = 1.40;

// The side-wall integrals of one axial index without J_m(gamma_q a)^2
// (ModalCavity.h): pi_p (pressure) and tau_p (tangential gradient).
void wallIntegrals (int m, int p, const AirSpec& spec, double& pressure, double& tangential)
{
    const double a = spec.radius, len = spec.length;
    const double eps = p == 0 ? 1.0 : 0.5;
    const double lam = angularNorm (m);
    const double lamSin = m == 0 ? 0.0 : pi;
    const double kz = p * pi / len;
    tangential = a * (double (m * m) / (a * a) * lamSin * len * eps + kz * kz * lam * len * (1.0 - eps));
    pressure = a * lam * len * eps;
}

// (c^2/2) sqrt(w nu/2) (tau/w^2 + (gamma-1) pi/(c^2 sqrt Pr)): the Kirchhoff
// power of a unit wall field (in units of rho c^2) at w.
double wallWeight (double pressure, double tangential, const AirSpec& spec, double w)
{
    const double c = spec.soundSpeed;
    const double nu = airViscosity / spec.density;
    return 0.5 * c * c * std::sqrt (w * nu / 2.0)
         * (tangential / (w * w) + (airGamma - 1.0) * pressure / (c * c * std::sqrt (airPrandtl)));
}

// coth and csch without overflow for large arguments.
double coth (double x) { return x > 20.0 ? 1.0 + 2.0 * std::exp (-2.0 * x) : 1.0 / std::tanh (x); }
double csch (double x) { return x > 20.0 ? 2.0 * std::exp (-x) : 1.0 / std::sinh (x); }
} // namespace

std::vector<double> besselPrimeZeros (int m, double xMax)
{
    std::vector<double> zeros;
    if (m == 0)
        zeros.push_back (0.0);
    const auto derivative = [m] (double x)
    {
        double j = 0.0, jp = 0.0;
        jAndPrime (m, x, j, jp);
        return jp;
    };
    constexpr double step = 0.25;
    // J_m' > 0 on (0, j'_m1) and j'_m1 > m for m >= 1, so the scan may start at m.
    double x0 = m == 0 ? step : static_cast<double> (m);
    double f0 = derivative (x0);
    for (double x1 = x0 + step; x1 <= xMax + step; x1 += step)
    {
        const double f1 = derivative (x1);
        if ((f0 < 0.0) != (f1 < 0.0) && f0 != 0.0)
        {
            double lo = x0, hi = x1, flo = f0;
            double x = 0.5 * (lo + hi);
            for (int iteration = 0; iteration < 60; ++iteration)
            {
                double j = 0.0, jp = 0.0;
                jAndPrime (m, x, j, jp);
                if ((jp < 0.0) == (flo < 0.0))
                {
                    lo = x;
                    flo = jp;
                }
                else
                    hi = x;
                // J'' = -J'/x - (1 - m^2/x^2) J (Bessel's equation).
                const double second = -jp / x - (1.0 - double (m * m) / (x * x)) * j;
                double next = second != 0.0 ? x - jp / second : 0.5 * (lo + hi);
                if (! (next > lo && next < hi))
                    next = 0.5 * (lo + hi);
                if (std::abs (next - x) < 1.0e-14 * x)
                {
                    x = next;
                    break;
                }
                x = next;
            }
            if (x < xMax)
                zeros.push_back (x);
        }
        x0 = x1;
        f0 = f1;
    }
    return zeros;
}

std::vector<CrossMode> crossModes (int m, double radius, double gammaMax)
{
    std::vector<CrossMode> out;
    for (double z : besselPrimeZeros (m, gammaMax * radius))
    {
        CrossMode c;
        c.gamma = z / radius;
        c.jmEdge = special::besselJ (m, z);
        const double shape = z == 0.0 ? 1.0 : 1.0 - double (m * m) / (z * z);
        c.normQ = angularNorm (m) * 0.5 * radius * radius * shape * c.jmEdge * c.jmEdge;
        out.push_back (c);
    }
    return out;
}

std::vector<AirMode> airModes (int m, const std::vector<CrossMode>& cross, const AirSpec& spec, double omegaMax)
{
    std::vector<AirMode> out;
    for (std::size_t q = 0; q < cross.size(); ++q)
        for (int p = 0;; ++p)
        {
            const double kz = p * pi / spec.length;
            const double omega = spec.soundSpeed * std::sqrt (cross[q].gamma * cross[q].gamma + kz * kz);
            if (omega > omegaMax)
                break;
            if (omega <= 0.0)
                continue;
            AirMode a;
            a.q = static_cast<int> (q);
            a.p = p;
            a.omega = omega;
            a.lambda = cross[q].normQ * spec.length * (p == 0 ? 1.0 : 0.5);
            out.push_back (a);
        }
    (void) m;
    return out;
}

InteriorBlock interiorBlock (int m, const std::vector<CrossMode>& cross, const std::vector<AirMode>& air,
                             const AirSpec& spec, const std::vector<std::vector<double>>& batterHankel,
                             const std::vector<std::vector<double>>& rearHankel)
{
    const int nb = static_cast<int> (batterHankel.size()), nr = static_cast<int> (rearHankel.size());
    const int nh = nb + nr;
    const double lam = angularNorm (m);
    const double rho = spec.density, c2 = spec.soundSpeed * spec.soundSpeed, len = spec.length;
    InteriorBlock out;
    out.heads = nh;
    out.air = static_cast<int> (air.size());
    out.staticStiffness.assign (static_cast<std::size_t> (nh * nh), 0.0);
    out.residualMass.assign (static_cast<std::size_t> (nh * nh), 0.0);
    out.coupling.assign (static_cast<std::size_t> (nh) * air.size(), 0.0);
    // Coupling vector of cavity mode (q, p): batter rows Lam W_i(gamma_q), rear rows (-1)^p Lam W_j(gamma_q).
    const auto couplingVector = [&] (std::size_t q, int p, std::vector<double>& c)
    {
        c.resize (static_cast<std::size_t> (nh));
        const double sign = (p & 1) ? -1.0 : 1.0;
        for (int i = 0; i < nb; ++i)
            c[static_cast<std::size_t> (i)] = lam * batterHankel[static_cast<std::size_t> (i)][q];
        for (int j = 0; j < nr; ++j)
            c[static_cast<std::size_t> (nb + j)] = sign * lam * rearHankel[static_cast<std::size_t> (j)][q];
    };
    const auto addOuter = [nh] (std::vector<double>& target, const std::vector<double>& c, double scale)
    {
        for (int i = 0; i < nh; ++i)
            for (int j = 0; j < nh; ++j)
                target[static_cast<std::size_t> (i * nh + j)] += scale * c[static_cast<std::size_t> (i)] * c[static_cast<std::size_t> (j)];
    };
    std::vector<double> c;
    for (std::size_t s = 0; s < air.size(); ++s)
    {
        couplingVector (static_cast<std::size_t> (air[s].q), air[s].p, c);
        addOuter (out.staticStiffness, c, rho * c2 / air[s].lambda);
        const double g = -air[s].omega * std::sqrt (rho * c2 / air[s].lambda);
        for (int i = 0; i < nh; ++i)
            out.coupling[static_cast<std::size_t> (i) * air.size() + s] = g * c[static_cast<std::size_t> (i)];
    }
    // The closed barrel's static compression (m = 0, q = 0, p = 0).
    if (m == 0 && ! cross.empty() && cross[0].gamma == 0.0)
    {
        couplingVector (0, 0, c);
        addOuter (out.staticStiffness, c, rho * c2 / (cross[0].normQ * len));
    }
    // Residual mass of the truncated modes, per cross-mode.
    for (std::size_t q = 0; q < cross.size(); ++q)
    {
        const double g = cross[q].gamma, nq = cross[q].normQ;
        double same = 0.0, across = 0.0;
        if (g == 0.0)
        {
            // Plane waves p >= 1: rho c^2 / (Lam_p Omega_p^2) = 2 rho L / (N_0 pi^2 p^2).
            same = pi * pi / 6.0;
            across = -pi * pi / 12.0;
            for (const auto& a : air)
                if (a.q == static_cast<int> (q))
                {
                    const double p2 = double (a.p) * a.p;
                    same -= 1.0 / p2;
                    across -= ((a.p & 1) ? -1.0 : 1.0) / p2;
                }
            const double scale = 2.0 * rho * len / (nq * pi * pi);
            same *= scale;
            across *= scale;
        }
        else
        {
            same = rho * coth (g * len) / (g * nq);
            across = rho * csch (g * len) / (g * nq);
            for (const auto& a : air)
                if (a.q == static_cast<int> (q))
                {
                    const double term = rho * c2 / (a.lambda * a.omega * a.omega);
                    same -= term;
                    across -= ((a.p & 1) ? -1.0 : 1.0) * term;
                }
        }
        for (int i = 0; i < nh; ++i)
        {
            const bool bi = i < nb;
            const double wi = bi ? batterHankel[static_cast<std::size_t> (i)][q] : rearHankel[static_cast<std::size_t> (i - nb)][q];
            for (int j = 0; j < nh; ++j)
            {
                const bool bj = j < nb;
                const double wj = bj ? batterHankel[static_cast<std::size_t> (j)][q] : rearHankel[static_cast<std::size_t> (j - nb)][q];
                out.residualMass[static_cast<std::size_t> (i * nh + j)] += (bi == bj ? same : across) * lam * lam * wi * wj;
            }
        }
    }
    return out;
}

double wallLossRate (int m, const CrossMode& cross, const AirMode& mode, const AirSpec& spec, double omega)
{
    double pressure = 0.0, tangential = 0.0;
    wallIntegrals (m, mode.p, spec, pressure, tangential);
    const double w = omega > 0.0 ? omega : mode.omega;
    // At the mode's own frequency the pressure amplitude is sqrt(rho c^2/Lam) Omega z;
    // at w it is sqrt(rho c^2/Lam) (w^2/Omega) z (ModalCavity.h), hence (w/Omega)^2.
    const double amplitude = (w / mode.omega) * (w / mode.omega);
    return wallWeight (pressure, tangential, spec, w) * cross.jmEdge * cross.jmEdge * amplitude / mode.lambda;
}

double kirchhoffPlaneWaveRate (const AirSpec& spec, double omega)
{
    const double nu = airViscosity / spec.density;
    return std::sqrt (omega * nu / 2.0) * (1.0 + (airGamma - 1.0) / std::sqrt (airPrandtl)) / spec.radius;
}

double wallLoss (int m, const std::vector<CrossMode>& cross, const std::vector<AirMode>& air, const AirSpec& spec,
                 const double* airCoordinates, double compression, double omega, std::vector<double>* half)
{
    int maxP = 0;
    for (const auto& s : air)
        maxP = std::max (maxP, s.p);
    std::vector<double> field (static_cast<std::size_t> (maxP + 1), 0.0);   // S_p
    const double w2 = omega * omega;
    for (std::size_t s = 0; s < air.size(); ++s)
    {
        const auto& mode = air[s];
        field[static_cast<std::size_t> (mode.p)] += cross[static_cast<std::size_t> (mode.q)].jmEdge * (w2 / mode.omega)
                                                    * airCoordinates[s] / std::sqrt (mode.lambda);
    }
    if (m == 0)
        field[0] += compression;   // J_0(0) = 1, cos(0) = 1
    double alpha = 0.0;
    for (int p = 0; p <= maxP; ++p)
    {
        double pressure = 0.0, tangential = 0.0;
        wallIntegrals (m, p, spec, pressure, tangential);
        const double weight = wallWeight (pressure, tangential, spec, omega) / w2;
        const double value = field[static_cast<std::size_t> (p)];
        alpha += weight * value * value;
        if (half != nullptr)
            half->push_back (std::sqrt (2.0 * weight) * value);
    }
    return alpha;
}

double evanescentInteriorMass (const HeadMode& mode, double radius, double k, double length, double density,
                               const std::vector<CrossMode>& cross)
{
    const double lam = angularNorm (mode.m);
    double sum = 0.0;
    for (const auto& q : cross)
    {
        if (q.gamma <= k)
            continue;
        const double kappa = std::sqrt (q.gamma * q.gamma - k * k);
        const double w = hankelTransform (mode, radius, q.gamma, q.jmEdge, 0.0);
        sum += density * coth (kappa * length) * lam * lam * w * w / (kappa * q.normQ);
    }
    return sum;
}

} // namespace taikor::candidates::modal
