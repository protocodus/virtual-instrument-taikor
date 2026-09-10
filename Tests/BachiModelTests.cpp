#include "DSP/BachiModel.h"

#include <array>
#include <cmath>
#include <iostream>

namespace
{
int failures = 0;

void expect (bool condition, const char* message)
{
    if (! condition)
    {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

bool near (double a, double b, double tolerance = 1.0e-11)
{
    return std::abs (a - b) <= tolerance * std::max ({ 1.0, std::abs (a), std::abs (b) });
}

// Independent reference for spatial tests, with arguments bounded below ten.
double bessel (int order, double x)
{
    double term = std::pow (0.5 * x, order) / std::tgamma (order + 1.0);
    double sum = term;
    for (int n = 1; n < 80; ++n)
    {
        term *= -x * x / (4.0 * n * (n + order));
        sum += term;
    }
    return sum;
}

void testFiniteContactAgainstSpatialIntegration()
{
    using namespace taikor::bachi;
    expect (diskProjection (0.0) == 1.0, "point-contact limit must preserve the mode");
    expect (near (diskProjection (1.0), 0.880101171489867),
            "disk projection must agree with the tabulated Bessel value");
    expect (std::abs (diskProjection (3.8317059702075125)) < 1.0e-12,
            "a footprint spanning the first cancelling lobe must have zero mean");
    expect (diskProjection (5.0) < 0.0,
            "projection beyond the first zero must retain its physical sign");

    constexpr double waveNumber = 10.0;
    constexpr double centerX = 0.31, centerY = -0.23, footprint = 0.11;
    for (int order = 0; order <= 5; ++order)
        for (int sineBranch = 0; sineBranch < 2; ++sineBranch)
        {
            const auto shape = [order, sineBranch] (double x, double y)
            {
                const double angle = order * std::atan2 (y, x);
                return bessel (order, waveNumber * std::hypot (x, y))
                     * (sineBranch == 0 ? std::cos (angle) : std::sin (angle));
            };
            const double integrated = averageEllipse (
                shape, centerX, centerY, footprint, footprint, 0.63);
            const double projected = shape (centerX, centerY)
                                   * diskProjection (waveNumber * footprint);
            expect (near (integrated, projected, 2.0e-11),
                    "disk factor must agree with area integration of off-axis Bessel modes");
        }

    // Independent midpoint area integration for large arguments, where the
    // implementation switches from a series to a weighted line quadrature.
    constexpr int radialCount = 8192, angularCount = 512;
    constexpr double pi = 3.1415926535897932384626433832795;
    for (double x : { 8.0, 20.0, 50.0, 79.0 })
    {
        double mean = 0.0;
        for (int radial = 0; radial < radialCount; ++radial)
        {
            const double r = std::sqrt ((radial + 0.5) / radialCount);
            for (int angular = 0; angular < angularCount; ++angular)
                mean += std::cos (x * r * std::cos (
                    2.0 * pi * (angular + 0.5) / angularCount));
        }
        mean /= radialCount * angularCount;
        expect (near (diskProjection (x), mean, 4.0e-5),
                "large-footprint projection must agree with independent disc integration");
        expect (near (diskProjection (-x), diskProjection (x)),
                "a disk average must be even in wave number");
    }
}

void testEllipseGeometry()
{
    using namespace taikor::bachi;
    constexpr double cx = 0.27, cy = -0.31, a = 0.08, b = 0.025, angle = 0.72;
    const double cosine = std::cos (angle), sine = std::sin (angle);
    const double meanX2 = cx * cx + 0.25 * (a * a * cosine * cosine + b * b * sine * sine);
    const double meanXY = cx * cy + 0.25 * (a * a - b * b) * sine * cosine;
    expect (near (averageEllipse ([] (double, double) { return 1.0; }, cx, cy, a, b, angle), 1.0),
            "finite contact must preserve total force for a constant mode");
    expect (near (averageEllipse ([] (double x, double) { return x * x; }, cx, cy, a, b, angle), meanX2),
            "rotated ellipse must have the correct second moment");
    expect (near (averageEllipse ([] (double x, double y) { return x * y; }, cx, cy, a, b, angle), meanXY),
            "rotated ellipse must preserve the signed cross moment");
    const auto stick = profileForFamily (0);
    expect (near (contactRadius (stick, 0.75, 0.4), stick.tipRadiusMetres),
            "contact footprint must keep the physical stick radius away from the rim");
    expect (near (contactRadius (stick, 0.75, 0.99), 0.0075),
            "rim contact must fit entirely within the active head");
    expect (contactRadius (stick, 0.75, 1.0) == 0.0,
            "a contact centered on the boundary must not extend outside the head");
}

void testStickMechanics()
{
    using namespace taikor::bachi;
    for (int family = 0; family < 4; ++family)
    {
        const auto stick = profileForFamily (family);
        expect (stick.massKg > 0.02f && stick.massKg < 0.6f,
                "family mass must describe one plausible bachi, not a pair");
        expect (stick.tipRadiusMetres > 0.0f && stick.lengthMetres > 0.0f,
                "reference stick dimensions must be physical");
        for (float hardness : { 0.0f, 0.2f, 0.5f, 1.0f })
            expect (contactExponent (stick.covering, hardness) == 1.5,
                    "soft bare wood must not acquire a compacting-felt exponent");
    }
    expect (contactExponent (Covering::feltWrapped, 0.0f) == 2.5,
            "explicit felt must retain its steeper soft-contact law");
    expect (near (profileForFamily (1).massKg, 0.330 / 2.0, 1.0e-7),
            "the Asano mass listed per pair must be divided by two");

    auto reference = profileForFamily (1);
    expect (near (contactStiffnessScale (reference), 1.0),
            "the reference stick must preserve the inherited contact coefficient");
    auto wider = reference;
    wider.tipRadiusMetres *= 4.0f;
    expect (near (contactStiffnessScale (wider), 2.0),
            "Hertz coefficient must grow with square root of tip curvature radius");
    auto softer = reference;
    softer.woodModulusPa *= 0.5f;
    expect (contactStiffnessScale (softer) < contactStiffnessScale (reference),
            "increased wood compliance must soften the contact");

    const auto first = profileForFamily (0), last = profileForFamily (3);
    const auto midpoint = blendProfiles (first, last, 0.5f);
    expect (near (midpoint.massKg, (first.massKg + last.massKg) * 0.5f, 1.0e-7),
            "family morph must keep independent physical mass between its endpoints");
}
} // namespace

int main()
{
    testFiniteContactAgainstSpatialIntegration();
    testEllipseGeometry();
    testStickMechanics();
    if (failures != 0)
        return 1;
    std::cout << "Bachi profiles, finite contact geometry and power-law mechanics passed.\n";
    return 0;
}
