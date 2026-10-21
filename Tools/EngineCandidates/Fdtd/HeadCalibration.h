#pragma once

// Measuring a staircase head's discrete eigenfrequencies against the analytic
// clamped stiff membrane, and the numerical correction of the staircase that
// follows from it.

#include "Analytic.h"
#include "HeadGrid.h"

#include <vector>

namespace taikor::candidates::fdtd
{

// The analytic mode (m, n) sampled on the grid's interior nodes, cos(m theta)
// or sin(m theta), zero on clamped nodes.
[[nodiscard]] std::vector<double> modeShapeOnGrid (const HeadGrid& head, const StiffMode& mode, bool sine);

// The frequency of the strongest spectral peak of trace (sampled at 1/k) in
// [lowHz, highHz]: Hann-windowed DTFT, coarse scan then golden-section refine.
[[nodiscard]] double spectralPeakFrequency (const std::vector<double>& trace, double k,
                                            double lowHz, double highHz);
// The same over several traces, their windowed powers summed.
[[nodiscard]] double spectralPeakFrequency (const std::vector<std::vector<double>>& traces, double k,
                                            double lowHz, double highHz);

// The discrete in-vacuo frequency of the head's (0,1) mode: lossless run from
// the analytic (0,1) shape, spectral peak of the displacement projected on it.
[[nodiscard]] double discreteFundamental (const HeadSpec& spec, double radius, double k,
                                          double spacing, double inclusionRadius);

struct InclusionCalibration
{
    double inclusionRadius = 0.0;
    double analyticHz = 0.0;       // clamped stiff (0,1)
    double uncorrectedHz = 0.0;    // discrete (0,1) with inclusionRadius = radius
    double correctedHz = 0.0;      // discrete (0,1) with the chosen inclusionRadius
};

// The staircase disk of nodes with r < a behaves as a slightly larger disk
// (about 0.4 h larger at a/h ~ 80; research note "staircase-disk eigen test"
// of the engine-candidate study found every mode 0.39-0.81 % flat), so every
// mode is flat by a nearly uniform fraction. As a numerical correction (the
// kind FdtdContracts.h allows for makeDiskGrid) the inclusion radius is moved
// inward until the discrete (0,1) matches the analytic clamped value: the
// discrete frequency is measured at two radii, interpolated linearly in the
// radius, and verified once at the chosen radius. The residual and the error
// of every other mode are reported by the verification suite (test 1).
[[nodiscard]] InclusionCalibration calibrateInclusionRadius (const HeadSpec& spec, double radius,
                                                             double k, double spacing);

} // namespace taikor::candidates::fdtd
