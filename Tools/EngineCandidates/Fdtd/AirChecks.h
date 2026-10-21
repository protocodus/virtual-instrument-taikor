#pragma once

// Measurements of the air domain's absorbing layer, shared by the engine's
// self-check (diagnostics) and the verification suite (test 5).

namespace taikor::candidates::fdtd
{

struct ReflectionResult
{
    double peakReflectionDb = 0.0;   // max |p_layer - p_reference| / max |p_incident|, dB
    double energyReflectionDb = 0.0; // energy of the difference / energy of the incident, dB
    double seconds = 0.0;
};

// A Gaussian pressure blob (standard deviation `blobCells` cells) released at
// rest `sourceCells` cells from a PML face, recorded at receivers between the
// source and the layer (normal and oblique paths), against the same blob in a
// domain large enough that nothing returns within the recording.
[[nodiscard]] ReflectionResult measurePmlReflection (double delta, double k, double c, double rho,
                                                     int pmlCells, double pmlReflection,
                                                     double blobCells, int threads);

} // namespace taikor::candidates::fdtd
