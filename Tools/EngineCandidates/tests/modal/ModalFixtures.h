#pragma once

// Drums and strokes for the tests. The nagado values are the ones the harness
// resolves from the shipping engine for the listening-test stroke (the engine
// itself always reads them from DrumSpec / StrokeSpec / MicSpec at run time).

#include "CandidateEngine.h"

#include <cmath>

namespace modaltest
{

inline taikor::candidates::DrumSpec nagadoDrum()
{
    taikor::candidates::DrumSpec d;
    d.octave = 1;
    d.radius = 0.39;
    d.depth = 0.936008;
    d.shellWallThickness = 0.05226;
    d.batter.tension = 6146.67;
    d.batter.arealDensity = 0.94;
    d.batter.thickness = 0.94e-3;
    d.batter.youngsModulus = 3.5e9;
    d.batter.poissonRatio = 0.3;
    d.batter.bendingStiffness = 0.266213;
    d.batter.sigma0 = 0.947723;
    d.batter.sigma1 = 0.00358265;
    d.rear.tension = 6269.61;
    d.rear.arealDensity = 0.89;
    d.rear.thickness = 0.89e-3;
    d.rear.youngsModulus = 3.4e9;
    d.rear.poissonRatio = 0.3;
    d.rear.bendingStiffness = 0.219496;
    d.rear.sigma0 = 0.97715;
    d.rear.sigma1 = 0.00400931;
    return d;
}

inline taikor::candidates::StrokeSpec nagadoStroke()
{
    taikor::candidates::StrokeSpec s;
    s.impactSpeed = 4.57943;
    s.strikeX = 0.0585;
    s.strikeY = 0.0;
    s.stickMass = 0.165;
    s.stickLength = 0.42;
    s.tipMass = 0.04125;
    s.footprintRadius = 0.01125;
    s.hertzStiffness = 1.08942e8;
    s.hertzExponent = 1.5;
    return s;
}

inline taikor::candidates::MicSpec nagadoMics()
{
    taikor::candidates::MicSpec m;
    m.capsule[0] = { 0.122351, 0.138576, 0.16 };
    m.capsule[1] = { 0.173493, 0.063822, 0.16 };
    return m;
}

// Suzuki & Hwang, "Coupling between two membranes of a Japanese drum",
// Acoust. Sci. & Tech. 29(3):215 (2008): R 0.20 m, L 0.50 m, cowhide
// 2000 kg/m^3 x 2 mm (4 kg/m^2), E 3.5 GPa, T 14 kN/m on both heads.
inline taikor::candidates::DrumSpec suzukiHwangDrum()
{
    taikor::candidates::DrumSpec d = nagadoDrum();
    d.radius = 0.20;
    d.depth = 0.50;
    for (auto* h : { &d.batter, &d.rear })
    {
        h->tension = 14000.0;
        h->arealDensity = 4.0;
        h->thickness = 0.002;
        h->youngsModulus = 3.5e9;
        h->poissonRatio = 0.3;
        h->bendingStiffness = 3.5e9 * 0.002 * 0.002 * 0.002 / (12.0 * (1.0 - 0.09));
    }
    return d;
}

} // namespace modaltest
