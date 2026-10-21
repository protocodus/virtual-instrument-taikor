#pragma once

// A small in-vacuo stiff membrane stepped with the explicit SAV scheme, used
// only to exercise the potential in time (tests 6, 6b and 7). The linear part
// is the drum's lossless stiff membrane (5-point Laplacian, clamped 13-point
// biharmonic L^T L); the nonlinear part enters as in Bilbao, Ducceschi &
// Zama, J. Comput. Phys. 472 (2023), eq. 34 in second-order form:
//   sigma h^2 (w+ - 2 w + w-)/k^2 = T h^2 L w - D h^2 L^T L w + f
//                                    - g (psi+ + psi-)/2,
//   psi+ - psi- = g . (w+ - w-) / 2,   g = grad V(w) / sqrt(2 V(w) + eps_r),
// solved for w+ with Sherman-Morrison. The conserved energy is
//   sigma h^2 |w+ - w|^2 / (2 k^2) - T h^2 <w+, L w>/2 + D h^2 <L w+, L w>/2
//   + psi+^2 / 2.
// The optional stick is a lumped tip mass on a one-sided Hertz spring over
// the footprint disk, stepped explicitly (it only drives the head here).
// The optional hide loss is the shared two-parameter law of the letters,
// -2 sigma sigma0 w_t + 2 sigma sigma1 Lap w_t (HeadSpec::sigma0, sigma1),
// with the sigma0 term centred (implicit, diagonal) and the sigma1 term
// backward, (L w^n - L w^{n-1}) / k, so the step stays explicit.

#include "TestUtil.h"

namespace vktest
{

struct InVacuoHead
{
    MembraneGrid g;
    HeadSpec head;
    double k = 0.0;               // time step, s
    NonlinearPotential* potential = nullptr;
    double epsR = 1.0e-16;        // SAV regulariser, J (numerical)
    bool loss = false;            // hide loss sigma0, sigma1 on (energy then decays)
    std::vector<double> wm, w, wp, lap, lap2, lapm, grad, gvec, rhs;
    double psi = 0.0;             // psi^{n-1/2}
    double lastV = 0.0;

    // Stick (tip mass on a Hertz spring over the footprint disk).
    bool stick = false;
    double stickMass = 0.0, hertzK = 0.0, hertzAlpha = 1.5;
    double zm = 0.0, z = 0.0;     // tip height at n-1 and n (+z out of the drum)
    std::vector<int> footprint;
    double contactForce = 0.0;

    void init (const MembraneGrid& grid, const HeadSpec& spec, double sampleRate, NonlinearPotential* v)
    {
        g = grid;
        head = spec;
        k = 1.0 / sampleRate;
        potential = v;
        const auto s = static_cast<std::size_t> (g.size());
        for (auto* a : { &wm, &w, &wp, &lap, &lap2, &lapm, &grad, &gvec, &rhs })
            a->assign (s, 0.0);
    }

    // Initial displacement (at rest): w^0 = w^1 = w0, psi^{1/2} = sqrt(2V + eps).
    void setDisplacement (const std::vector<double>& w0)
    {
        w = w0;
        wm = w0;
        psi = potential != nullptr ? std::sqrt (2.0 * potential->evaluate (w.data(), grad.data()) + epsR) : 0.0;
    }

    void armStick (const StrokeSpec& s)
    {
        stick = true;
        stickMass = s.tipMass;
        hertzK = s.hertzStiffness;
        hertzAlpha = s.hertzExponent;
        footprint.clear();
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i)
                if (g.interior[static_cast<std::size_t> (g.index (i, j))] != 0
                    && std::hypot (g.x (i) - s.strikeX, g.y (j) - s.strikeY) < s.footprintRadius)
                    footprint.push_back (g.index (i, j));
        z = 0.0;
        zm = s.impactSpeed * k; // moving into the drum (-z)
    }

    // L w on every node (w zero off the interior), the 5-point Laplacian.
    void laplacian (const std::vector<double>& in, std::vector<double>& out) const
    {
        const double ih2 = 1.0 / (g.h * g.h);
        std::fill (out.begin(), out.end(), 0.0);
        for (int j = 1; j + 1 < g.ny; ++j)
            for (int i = 1; i + 1 < g.nx; ++i)
            {
                const int q = g.index (i, j);
                out[static_cast<std::size_t> (q)] = ih2 * (in[static_cast<std::size_t> (q - 1)] + in[static_cast<std::size_t> (q + 1)]
                    + in[static_cast<std::size_t> (q - g.nx)] + in[static_cast<std::size_t> (q + g.nx)] - 4.0 * in[static_cast<std::size_t> (q)]);
            }
    }

    double footprintMean (const std::vector<double>& field) const
    {
        double s = 0.0;
        for (int q : footprint)
            s += field[static_cast<std::size_t> (q)];
        return footprint.empty() ? 0.0 : s / static_cast<double> (footprint.size());
    }

    void step()
    {
        const auto s = static_cast<std::size_t> (g.size());
        const double h2 = g.h * g.h;
        const double mass = head.arealDensity * h2 / (k * k);
        const double s0 = loss ? head.sigma0 * k : 0.0; // sigma0 k
        const double alpha = mass * (1.0 + s0);
        laplacian (w, lap);
        laplacian (lap, lap2);
        if (loss)
            laplacian (wm, lapm);
        // SAV vector g at w^n.
        double gw = 0.0, gg = 0.0;
        if (potential != nullptr)
        {
            lastV = potential->evaluate (w.data(), grad.data());
            const double inv = 1.0 / std::sqrt (2.0 * lastV + epsR);
            for (std::size_t q = 0; q < s; ++q)
            {
                gvec[q] = grad[q] * inv;
                gw += gvec[q] * wm[q];
                gg += gvec[q] * gvec[q];
            }
        }
        contactForce = 0.0;
        if (stick)
        {
            const double eta = footprintMean (w) - z; // head surface above the tip = compression
            contactForce = eta > 0.0 ? hertzK * std::pow (eta, hertzAlpha) : 0.0;
            const double zNew = 2.0 * z - zm + k * k * contactForce / stickMass;
            zm = z;
            z = zNew;
        }
        const double perNode = stick && ! footprint.empty() ? contactForce / static_cast<double> (footprint.size()) : 0.0;
        for (std::size_t q = 0; q < s; ++q)
        {
            if (g.interior[q] == 0)
            {
                rhs[q] = 0.0;
                continue;
            }
            rhs[q] = mass * (2.0 * w[q] - (1.0 - s0) * wm[q]) + head.tension * h2 * lap[q] - head.bendingStiffness * h2 * lap2[q]
                   - gvec[q] * psi + 0.25 * gvec[q] * gw;
            if (loss)
                rhs[q] += 2.0 * head.arealDensity * head.sigma1 * h2 * (lap[q] - lapm[q]) / k;
        }
        for (int q : footprint)
            rhs[static_cast<std::size_t> (q)] -= perNode;
        // (alpha I + g g^T / 4)^-1 rhs by Sherman-Morrison.
        double gr = 0.0;
        for (std::size_t q = 0; q < s; ++q)
            gr += gvec[q] * rhs[q];
        const double factor = 0.25 * gr / (alpha * (alpha + 0.25 * gg));
        double dpsi = 0.0;
        for (std::size_t q = 0; q < s; ++q)
        {
            wp[q] = g.interior[q] != 0 ? rhs[q] / alpha - factor * gvec[q] : 0.0;
            dpsi += gvec[q] * (wp[q] - wm[q]);
        }
        psi += 0.5 * dpsi; // psi^{n+1/2}
        std::swap (wm, w);
        std::swap (w, wp); // now wm = w^n, w = w^{n+1}, wp = w^{n-1}
    }

    // Energy at n+1/2 after a step (wm = w^n, w = w^{n+1}), stick excluded.
    double energy()
    {
        const double h2 = g.h * g.h;
        laplacian (wm, lap);
        laplacian (w, lap2);
        double kin = 0.0, mem = 0.0, bend = 0.0;
        for (int q = 0; q < g.size(); ++q)
        {
            const auto u = static_cast<std::size_t> (q);
            const double dw = w[u] - wm[u];
            kin += dw * dw;
            mem -= w[u] * lap[u];
            bend += lap2[u] * lap[u];
        }
        return 0.5 * head.arealDensity * h2 * kin / (k * k) + 0.5 * head.tension * h2 * mem
             + 0.5 * head.bendingStiffness * h2 * bend + 0.5 * psi * psi;
    }

    // Smallest h for which the explicit stiff membrane is stable at this k
    // (von Neumann, largest 5-point eigenvalue 8 / h^2):
    //   k^2 (8 c^2 / h^2 + 64 kappa^2 / h^4) + 32 sigma1 k / h^2 <= 4,
    // c^2 = T / sigma, kappa^2 = D / sigma; the sigma1 term only with loss.
    static double minimumSpacing (const HeadSpec& spec, double sampleRate, bool withLoss = false)
    {
        const double k = 1.0 / sampleRate;
        const double a = 2.0 * spec.tension * k * k / spec.arealDensity + (withLoss ? 8.0 * spec.sigma1 * k : 0.0);
        const double b = 16.0 * spec.bendingStiffness * k * k / spec.arealDensity;
        // h^4 - a h^2 - b >= 0  ->  h^2 >= (a + sqrt(a^2 + 4 b)) / 2
        return std::sqrt (0.5 * (a + std::sqrt (a * a + 4.0 * b)));
    }
};

} // namespace vktest
