#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace taikor::bachi
{
enum class Covering { bareWood, feltWrapped };

struct Profile
{
    float massKg = 0.165f;
    float tipRadiusMetres = 0.01125f;
    float lengthMetres = 0.42f;
    float woodModulusPa = 12.3e9f;
    Covering covering = Covering::bareWood;
};

// These are explicit reference sticks, not compulsory traditional pairings or
// measured contact calibrations. Mass is PER STICK (some makers quote pairs).
// Odaiko: hinoki 36 x 540 mm, 230-280 g; use the range's midpoint.
// https://taiko-shop.com/products/hinoki-odaiko-3654-bachi
// Nagado: Asano K-7514, 22.5 x 420 mm, 330 g per pair.
// https://asano.us/products/kashi-bachi
// Okedo: ho 18 x 380 mm, 40-55 g; use the range's midpoint.
// https://taiko-shop.com/products/katsugi-oke-1838-bachi-material-ho
// Shime: maker's short maple 3/4 x 14 1/4 inches, 4 oz per pair.
// https://users.lmi.net/taikousa/diyprice.html
// Longitudinal modulus priors: hinoki 6 GPa, Table 1 of
// https://www.jstage.jst.go.jp/article/aem/1/0/1_143/_pdf
// Ho 11.1-12.2 GPa static test range (midpoint below), Table 1 of
// https://www.sciencedirect.com/science/article/pii/S0167844225003672
// White oak 12.3 GPa / sugar maple 12.6 GPa, USDA Wood Handbook Table 5-3a:
// https://www.fpl.fs.usda.gov/documnts/fplgtr/fpl_gtr190.pdf
// The oak/maple species are reference material proxies; bulk axial moduli are
// NOT measured end-grain indentation moduli for these particular bachi. Their
// relative compliance is used below; the inherited absolute contact scale is
// still a model calibration, not a measured pressure/force match.
[[nodiscard]] constexpr Profile profileForFamily (int family) noexcept
{
    constexpr std::array<Profile, 4> profiles {{
        { 0.255f, 0.018f, 0.54f, 6.0e9f, Covering::bareWood },
        { 0.165f, 0.01125f, 0.42f, 12.3e9f, Covering::bareWood },
        { 0.0475f, 0.009f, 0.38f, 11.65e9f, Covering::bareWood },
        { 0.056699046f, 0.009525f, 0.36195f, 12.6e9f, Covering::bareWood }
    }};
    return profiles[static_cast<std::size_t> (std::clamp (family, 0, 3))];
}

[[nodiscard]] inline Profile blendProfiles (const Profile& from, const Profile& to,
                                             float amount) noexcept
{
    const float mix = std::isfinite (amount) ? std::clamp (amount, 0.0f, 1.0f) : 0.0f;
    const auto blend = [mix] (float a, float b) { return a + mix * (b - a); };
    return { blend (from.massKg, to.massKg),
             blend (from.tipRadiusMetres, to.tipRadiusMetres),
             blend (from.lengthMetres, to.lengthMetres),
             blend (from.woodModulusPa, to.woodModulusPa),
             mix < 0.5f ? from.covering : to.covering };
}

// Hertz K is proportional to reduced modulus times sqrt(tip curvature radius).
// Treat the rounded-end radius as the cylinder radius: an explicit geometric
// approximation, not a measured curvature. Equal Poisson ratios cancel in this
// relative ratio. The inherited hardness curve supplies the absolute K.
[[nodiscard]] inline float contactStiffnessScale (const Profile& stick,
                                                  float headModulusPa = 3.5e9f) noexcept
{
    const auto reference = profileForFamily (1);
    const double head = std::max (static_cast<double> (headModulusPa), 1.0);
    const double wood = std::max (static_cast<double> (stick.woodModulusPa), 1.0);
    const double reduced = head * wood / (head + wood);
    const double referenceReduced = head * reference.woodModulusPa
                                  / (head + reference.woodModulusPa);
    return static_cast<float> (reduced / referenceReduced * std::sqrt (
        std::max (static_cast<double> (stick.tipRadiusMetres), 1.0e-6)
            / reference.tipRadiusMetres));
}

[[nodiscard]] inline double contactExponent (Covering covering, float hardness) noexcept
{
    if (covering == Covering::bareWood)
        return 1.5;
    const double amount = std::isfinite (hardness)
        ? std::clamp (static_cast<double> (hardness), 0.0, 1.0) : 0.5;
    return 2.5 - amount;
}

// Normalized mean of a Helmholtz membrane mode over a circular contact patch:
// <phi>_disk = phi(center) * 2 J1(k b)/(k b). This follows by angular averaging
// the Bessel addition theorem and integrating r J0(k r). The SAME signed factor
// must multiply force distribution AND displacement sensing for reciprocity.
// Never take abs(): a footprint can cross a modal lobe and reverse its mean.
// This assumes the patch is wholly on the active membrane. Callers must bound
// b to the available rim clearance; this is not an outside-the-head extrapolation.
[[nodiscard]] inline double diskProjection (double waveNumberTimesRadius) noexcept
{
    const double x = std::abs (waveNumberTimesRadius);
    if (! std::isfinite (x))
        return 0.0;
    if (x < 8.0)
    {
        // 2 J1(x)/x = sum (-x^2/4)^n / (n! (n+1)!). Stable at x=0.
        const double squaredQuarter = 0.25 * x * x;
        double term = 1.0, sum = 1.0;
        for (int n = 1; n <= 32; ++n)
        {
            term *= -squaredQuarter / static_cast<double> (n * (n + 1));
            sum += term;
            if (std::abs (term) < 1.0e-17)
                break;
        }
        return sum;
    }

    // Gauss-Chebyshev quadrature of 2/pi int_-1^1 sqrt(1-u^2) cos(xu) du.
    // 64 nodes avoid the cancellation of the ascending Bessel series at the
    // largest resolved modal arguments (|k b| < 80). Geometry/setup only.
    constexpr int nodes = 64;
    constexpr double pi = 3.1415926535897932384626433832795;
    double sum = 0.0;
    for (int node = 1; node <= nodes; ++node)
    {
        const double angle = pi * static_cast<double> (node) / (nodes + 1);
        const double sine = std::sin (angle);
        sum += sine * sine * std::cos (x * std::cos (angle));
    }
    return 2.0 * sum / (nodes + 1);
}

[[nodiscard]] inline double contactRadius (const Profile& stick, double headRadius,
                                            double normalizedStrikeRadius) noexcept
{
    if (! std::isfinite (headRadius) || ! (headRadius > 0.0)
        || ! std::isfinite (normalizedStrikeRadius))
        return 0.0;
    const double clearance = headRadius
        * (1.0 - std::clamp (normalizedStrikeRadius, 0.0, 1.0));
    return std::clamp (static_cast<double> (stick.tipRadiusMetres), 0.0, clearance);
}

// General finite elliptical footprint for future oblique/tapered-tip contacts.
// The callback evaluates the full signed mode at Cartesian head coordinates;
// use one result on both sides of the mechanical transformer. Eight-point
// Gauss-Legendre in squared radius and 32 azimuths integrate the area mean.
// No allocation, and no callback invocation outside the requested ellipse.
template <typename ModalShape>
[[nodiscard]] double averageEllipse (ModalShape&& shape, double centerX,
                                      double centerY, double radiusX,
                                      double radiusY, double angle) noexcept
{
    constexpr std::array<double, 8> nodes {
        -0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
        -0.1834346424956498, 0.1834346424956498, 0.5255324099163290,
        0.7966664774136267, 0.9602898564975363
    };
    constexpr std::array<double, 8> weights {
        0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
        0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
        0.2223810344533745, 0.1012285362903763
    };
    constexpr double twoPi = 6.283185307179586476925286766559;
    const double cosine = std::cos (angle), sine = std::sin (angle);
    double sum = 0.0;
    for (std::size_t radial = 0; radial < nodes.size(); ++radial)
    {
        const double radius = std::sqrt (0.5 * (nodes[radial] + 1.0));
        for (int azimuth = 0; azimuth < 32; ++azimuth)
        {
            const double theta = twoPi * static_cast<double> (azimuth) / 32.0;
            const double x = radiusX * radius * std::cos (theta);
            const double y = radiusY * radius * std::sin (theta);
            sum += weights[radial] * shape (centerX + cosine * x - sine * y,
                                            centerY + sine * x + cosine * y);
        }
    }
    return sum / 64.0;
}
} // namespace taikor::bachi
