#include "DSP/PhysicalDrumProfile.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace
{
int failures = 0;
void expect (bool condition, const char* message)
{
    if (! condition)
    {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}
bool near (float a, float b, float tolerance = 2.0e-6f)
{
    return std::abs (a - b) <= tolerance * std::max ({ 1.0f, std::abs (a), std::abs (b) });
}

void testMaterialSeparation()
{
    using namespace taikor::physical;
    const auto oka = resolveHeads (2, 1.0f, 0.75f);
    const auto shime = resolveHeads (3, 1.0f, 0.75f);
    expect (oka.batter.thickness < shime.batter.thickness,
            "thin okedo and thick tsuke-shime must be independent natural-hide designs");
    expect (near (oka.batter.arealDensity / oka.batter.thickness, 1000.0f)
                && near (shime.batter.arealDensity / shime.batter.thickness, 1000.0f),
            "a thin natural family head must not inherit synthetic-film density");
    expect (oka.rear.thickness > oka.batter.thickness,
            "an independently specified rear head may be thicker than the batter");
    expect (! near (oka.rear.lossAngle, oka.batter.lossAngle, 1.0e-6f),
            "rear losses must be independently specified");

    const auto thin = resolveHeads (2, 1.0f, 0.30f).batter;
    const auto thick = resolveHeads (2, 1.0f, 0.90f).batter;
    expect (near (thin.youngsModulus, thick.youngsModulus),
            "thickness trim must not silently change natural-hide identity");
    const float ratio = thick.thickness / thin.thickness;
    expect (near (flexuralRigidity (thick) / flexuralRigidity (thin), ratio * ratio * ratio),
            "bending rigidity must follow the independent physical thickness cubed");
    expect (near (thick.arealDensity / thin.arealDensity, ratio),
            "mass must scale linearly with thickness, not with bending rigidity");
    const auto film = resolveHeads (2, 1.0f, 0.0f).batter;
    expect (near (film.arealDensity / film.thickness, 1390.0f),
            "the host minimum must still provide the explicit synthetic-film endpoint");
    for (int family = 0; family < 4; ++family)
        for (int trim = 0; trim <= 100; ++trim)
        {
            const auto heads = resolveHeads (family, 1.0f, static_cast<float> (trim) / 100.0f);
            for (const auto& head : { heads.batter, heads.rear })
                expect (std::isfinite (flexuralRigidity (head)) && flexuralRigidity (head) > 0.0f
                            && head.arealDensity > 0.0f && head.lossAngle > 0.0f
                            && head.viscousLossSeconds > 0.0f,
                        "all host material settings must describe passive finite heads");
        }
}

void testSpatialPerturbation()
{
    using namespace taikor::physical;
    for (int family = 0; family < 4; ++family)
        for (int order = 1; order <= 12; ++order)
            for (int radial = 1; radial <= 8; ++radial)
            {
                const auto basis = basisForMode (family, 1.0f, order, radial);
                expect (basis.frequencyScale[0] > 1.0f && basis.frequencyScale[1] < 1.0f,
                        "a principal pair must split around its unperturbed mode");
                expect (near (basis.frequencyScale[0] * basis.frequencyScale[0]
                                + basis.frequencyScale[1] * basis.frequencyScale[1], 2.0f),
                        "the perturbation must preserve the pair stiffness trace");
                for (int step = 0; step < 31; ++step)
                {
                    const float angle = static_cast<float> (step) * 0.213f;
                    const float a = basis.project (angle, 0), b = basis.project (angle, 1);
                    expect (near (a * a + b * b, 1.0f),
                            "rotating nodal diameters must preserve projection energy");
                    const auto again = basisForMode (family, 1.0f, order, radial);
                    expect (a == again.project (angle, 0) && b == again.project (angle, 1),
                            "a physical head basis must never drift between strikes");
                }
                const float theta = 0.27f, phi = -0.49f;
                const float covariance = basis.project (theta, 0) * basis.project (phi, 0)
                                       + basis.project (theta, 1) * basis.project (phi, 1);
                expect (near (covariance, std::cos (order * (theta - phi))),
                        "strike and microphone projections must share a reciprocal orthogonal basis");
            }
    const auto radial = basisForMode (2, 1.0f, 0, 4);
    expect (radial.project (0.7f, 0) == 1.0f && radial.project (0.7f, 1) == 0.0f,
            "an axisymmetric shape must remain angle-independent and unsplit");
    const auto reference = basisForMode (0, 1.0f, 1, 1);
    const auto tuned = basisForMode (3, 0.0f, 1, 1);
    expect (reference.frequencyScale == tuned.frequencyScale && reference.rotation == tuned.rotation,
            "1 Drum layout must retain the reference head's physical imperfections");
    const auto family = basisForMode (2, 1.0f, 1, 1);
    expect (reference.rotation != family.rotation && reference.frequencyScale != family.frequencyScale,
            "different family constructions must not reuse one global imperfection seed");
}
} // namespace

int main()
{
    testMaterialSeparation();
    testSpatialPerturbation();
    if (failures == 0)
        std::cout << "Physical drum profile tests passed\n";
    return failures == 0 ? 0 : 1;
}
