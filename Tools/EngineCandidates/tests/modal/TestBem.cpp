// The axisymmetric exterior solver (Modal/ModalBem.h) against exact fields:
// a sphere vibrating in a sectoral harmonic (every order m the drum uses, with
// the analytic spherical-Hankel impedance), and a ring source inside the drum
// body, whose field outside is known exactly, at the body's first interior
// Dirichlet frequencies, where the plain equation fails and CHIEF must hold.

#include "ModalBem.h"
#include "ModalMath.h"
#include "ModalTestUtil.h"

#include <cmath>
#include <complex>

using namespace taikor::candidates::modal;
using cplx = std::complex<double>;

namespace
{
// Spherical Hankel function of the second kind h_n(x) = j_n - j y_n and its derivative.
void sphericalHankel2 (int n, double x, cplx& h, cplx& hp)
{
    std::vector<cplx> seq (static_cast<std::size_t> (n + 2));
    seq[0] = cplx (0.0, 1.0) * std::exp (cplx (0.0, -x)) / x;
    // h_1 = (sin x / x^2 - cos x / x) - j (-cos x / x^2 - sin x / x).
    seq[1] = cplx (std::sin (x) / (x * x) - std::cos (x) / x, std::cos (x) / (x * x) + std::sin (x) / x);
    for (int l = 1; l <= n; ++l)
        seq[static_cast<std::size_t> (l + 1)] = (2.0 * l + 1.0) / x * seq[static_cast<std::size_t> (l)] - seq[static_cast<std::size_t> (l - 1)];
    h = seq[static_cast<std::size_t> (n)];
    hp = n == 0 ? -seq[1] : seq[static_cast<std::size_t> (n - 1)] - (n + 1.0) / x * h;
}

BemGeometry sphere (double radius, int segments, int chiefPoints)
{
    std::vector<std::array<double, 2>> nodes;
    for (int i = 0; i <= segments; ++i)
    {
        const double theta = pi * i / segments;
        nodes.push_back ({ radius * std::sin (theta), radius * std::cos (theta) });
    }
    nodes.front()[0] = nodes.back()[0] = 0.0;
    std::vector<std::array<double, 2>> chief;
    for (int i = 0; i < chiefPoints; ++i)
    {
        // Off-axis interior points away from symmetry planes.
        const double u = (i + 0.37) / chiefPoints;
        chief.push_back ({ radius * (0.15 + 0.6 * std::fmod (u * 7.31, 1.0)), radius * (0.8 * u - 0.43) });
    }
    return polylineGeometry (nodes, std::vector<int> (static_cast<std::size_t> (segments), 0),
                             std::vector<double> (static_cast<std::size_t> (segments), 1.0), chief);
}

// Relative rms error of the BEM surface pressure of a sphere vibrating as
// sin^m(theta) cos(m phi), against -j rho c h_m(ka)/h_m'(ka) v.
double sphereError (int m, double ka, bool chief, int segments = 160)
{
    const double radius = 0.4, rho = 1.2, c = 343.0, k = ka / radius;
    BemSolver solver (sphere (radius, segments, 24));
    BemProblem pr;
    pr.m = m;
    pr.activeParts = { 0 };
    pr.columns = 1;
    for (const auto& panel : solver.geometry().panels)
        for (double t : bemNodeOffsets (panel))
            pr.velocity.push_back (std::pow (panel.r (t) / radius, m));
    std::vector<BemProblem> problems { pr };
    solver.solve (k, rho, c, problems, chief, 8);
    cplx h, hp;
    sphericalHankel2 (m, ka, h, hp);
    const cplx ratio = cplx (0.0, -rho * c) * h / hp;
    double err = 0.0, ref = 0.0;
    for (std::size_t j = 0; j < solver.geometry().panels.size(); ++j)
    {
        const auto& panel = solver.geometry().panels[j];
        const double rm = panel.r (0.5 * panel.length), zm = panel.z (0.5 * panel.length);
        const cplx exact = ratio * std::pow (rm / std::hypot (rm, zm), m);
        err += std::norm (problems[0].pressure[j] - exact) * panel.length;
        ref += std::norm (exact) * panel.length;
    }
    return std::sqrt (err / ref);
}
} // namespace

MODAL_TEST (bem_sphere)
{
    // 160 panels (h = 7.9 mm on R = 0.4 m); the polygon and the constant
    // pressure per panel leave an O(h^2) error relative to the width of the
    // pattern, sin^m(theta) ~ exp(-m (theta - pi/2)^2 / 2): 400 panels at m = 60.
    for (const int m : { 0, 1, 3, 8, 30, 60 })
        for (const double ka : { 0.3, 2.0, 6.0 })
        {
            const double e = sphereError (m, ka, true, m > 30 ? 400 : 160);
            std::printf ("  m %d, ka %.1f: relative rms error %.2e\n", m, ka, e);
            const std::string what = "sphere m " + std::to_string (m) + " ka " + std::to_string (ka).substr (0, 3);
            modaltest::below (what.c_str(), e, 2e-3);
        }
    // The first interior Dirichlet frequency of the sphere for m = 0 (j_0(ka) = 0,
    // ka = pi) and m = 1 (j_1(ka) = 0, ka = 4.4934): the plain equation is singular.
    for (const auto& [m, ka] : { std::pair<int, double> { 0, pi }, std::pair<int, double> { 1, 4.493409458 } })
    {
        const double plain = sphereError (m, ka, false), chief = sphereError (m, ka, true);
        std::printf ("  m %d at the interior resonance ka %.4f: plain %.2e, with CHIEF %.2e\n", m, ka, plain, chief);
        modaltest::below ("sphere at an interior resonance, with CHIEF", chief, 2e-3);
    }
}

namespace
{
// The drum body with a ring source of order m inside it at (0.2, -0.4) m: its
// field g_m(source, y) is an exact exterior solution, so the solver given its
// normal derivative must return it on the surface (relative rms error).
double drumRingError (int m, double k, bool chief)
{
    const double a = 0.39, wall = 0.05226, len = 0.936008, rho = 1.2041, c = 343.0;
    BemSolver solver (drumGeometry (a, wall, len, 0.005, 0.01, 40));
    const double r0 = 0.2, z0 = -0.4;
    BemProblem pr;
    pr.m = m;
    pr.activeParts = { 0, 1, 2 };
    pr.columns = 2; // real and imaginary parts of v = dp/dn / (-j w rho)
    std::vector<cplx> g (static_cast<std::size_t> (m + 1)), h (g.size());
    for (const auto& panel : solver.geometry().panels)
        for (double t : bemNodeOffsets (panel))
        {
            ringKernels (r0, z0, panel.r (t), panel.z (t), panel.nr, panel.nz, k, m, 2, g.data(), h.data());
            const cplx v = h[static_cast<std::size_t> (m)] / cplx (0.0, -k * c * rho);
            pr.velocity.push_back (v.real());
            pr.velocity.push_back (v.imag());
        }
    std::vector<BemProblem> problems { pr };
    solver.solve (k, rho, c, problems, chief, 8);
    double err = 0.0, ref = 0.0;
    for (std::size_t j = 0; j < solver.geometry().panels.size(); ++j)
    {
        const auto& panel = solver.geometry().panels[j];
        ringKernels (r0, z0, panel.r (0.5 * panel.length), panel.z (0.5 * panel.length), panel.nr, panel.nz, k, m, 2,
                     g.data(), h.data());
        const cplx exact = g[static_cast<std::size_t> (m)];
        const cplx got = problems[0].pressure[2 * j] + cplx (0.0, 1.0) * problems[0].pressure[2 * j + 1];
        err += std::norm (got - exact) * panel.r (0.5 * panel.length) * panel.length;
        ref += std::norm (exact) * panel.r (0.5 * panel.length) * panel.length;
    }
    return std::sqrt (err / ref);
}
} // namespace

MODAL_TEST (bem_drum_body)
{
    // The body (outer radius b = a + wall, length L) is singular for the plain
    // equation at its interior Dirichlet eigenvalues k^2 = (j_m1 / b)^2 + (pi / L)^2.
    const double b = 0.39 + 0.05226, len = 0.936008;
    const double j01 = 2.404825557695773, j11 = 3.831705970207512;
    const double k0 = std::sqrt ((j01 / b) * (j01 / b) + (pi / len) * (pi / len));
    const double k1 = std::sqrt ((j11 / b) * (j11 / b) + (pi / len) * (pi / len));
    for (const auto& [m, k] : { std::pair<int, double> { 0, 1.2 }, std::pair<int, double> { 3, 20.0 },
                                std::pair<int, double> { 0, k0 }, std::pair<int, double> { 1, k1 } })
    {
        const double plain = drumRingError (m, k, false), chief = drumRingError (m, k, true);
        std::printf ("  m %d, k %.3f 1/m (%.1f Hz): relative rms error plain %.2e, with CHIEF %.2e\n", m, k,
                     k * 343.0 / (2.0 * pi), plain, chief);
        const std::string what = "drum body ring source m " + std::to_string (m) + ", with CHIEF";
        modaltest::below (what.c_str(), chief, 2e-3);
    }
}

MODAL_TEST (bem_compact_body)
{
    // At ka << 1 the closed drum is a compact body and each head's outward volume
    // velocity Q radiates as a free-field monopole, R = rho c k^2 Q^2 / (4 pi) (half
    // the baffled piston's rho c k^2 Q^2 / (2 pi)), and the two heads' fields add:
    // R_br -> the same value (A. D. Pierce, Acoustics, 3rd ed., 2019, sec. 4.3).
    // Unit outward piston velocity on each head in turn, at 5 Hz (k L = 0.086).
    const double a = 0.39, wall = 0.05226, len = 0.936008, rho = 1.2041, c = 343.0;
    const double k = 2.0 * pi * 5.0 / c;
    BemSolver solver (drumGeometry (a, wall, len, 0.01, 0.02, 40));
    BemProblem pr;
    pr.m = 0;
    pr.chief = false;
    pr.activeParts = { 0, 1 };
    pr.columns = 2;
    for (const auto& panel : solver.geometry().panels)
        if (panel.part <= 1)
            for (int q = 0; q < bemVelocityNodes; ++q)
            {
                pr.velocity.push_back (panel.part == 0 ? 1.0 : 0.0);
                pr.velocity.push_back (panel.part == 1 ? 1.0 : 0.0);
            }
    std::vector<BemProblem> problems { pr };
    solver.solve (k, rho, c, problems, false, 8);
    cplx z[2][2] = { { 0.0, 0.0 }, { 0.0, 0.0 } };
    for (std::size_t j = 0; j < solver.geometry().panels.size(); ++j)
    {
        const auto& p = solver.geometry().panels[j];
        if (p.part > 1)
            continue;
        const double area = 2.0 * pi * p.r (0.5 * p.length) * p.length; // piston: int p dS
        for (int col = 0; col < 2; ++col)
            z[p.part][col] += area * problems[0].pressure[j * 2 + static_cast<std::size_t> (col)];
    }
    const double q = pi * a * a;
    const double monopole = rho * c * k * k * q * q / (4.0 * pi);
    std::printf ("  R_bb / monopole %.4f, R_rr / monopole %.4f, R_br / monopole %.4f, R_rb / monopole %.4f\n",
                 z[0][0].real() / monopole, z[1][1].real() / monopole, z[0][1].real() / monopole, z[1][0].real() / monopole);
    std::printf ("  push-pull (dipole) resistance (R_bb + R_rr - R_br - R_rb) / monopole %.2e\n",
                 (z[0][0] + z[1][1] - z[0][1] - z[1][0]).real() / monopole);
    modaltest::nearRel ("R_bb at 5 Hz / free-field monopole", z[0][0].real(), monopole, 0.02);
    modaltest::nearRel ("R_rr at 5 Hz / free-field monopole", z[1][1].real(), monopole, 0.02);
    modaltest::nearRel ("R_br at 5 Hz / free-field monopole", z[0][1].real(), monopole, 0.02);
    modaltest::below ("reciprocity |Z_br - Z_rb| / |Z_bb|", std::abs (z[0][1] - z[1][0]) / std::abs (z[0][0]), 2e-3);
}
