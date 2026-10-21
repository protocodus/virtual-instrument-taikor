// The barrel wall's viscothermal boundary layer: finding the wall on the
// staircase grid (see WallLoss.h for the model, the scheme and its energy;
// WallLossStep.cpp steps it).

#include "WallLoss.h"

#include <algorithm>
#include <cmath>
#include <map>

namespace taikor::candidates::fdtd
{

template <typename Real>
void WallLoss<Real>::build (const AirDomain<Real>& air, double radius, double depth, const SquareRootFit& fit,
                            const Viscothermal& gas)
{
    domain = &air;
    xi = fit.xi;
    muThermal = fit.thermal;
    muViscous = fit.viscous;
    k = air.k;
    rho = air.density;
    c = air.c;
    delta = air.delta;
    lambda2 = air.lambda * air.lambda;
    const double nu = gas.viscosity / rho, chi = nu / gas.prandtl;
    rootNu = std::sqrt (nu);
    cTh = (gas.gamma - 1.0) * std::sqrt (chi) / (rho * c * c);
    cells.clear();
    faces.clear();
    elements.clear();
    elementFaces.clear();
    elementWeights.clear();
    cellFaces.clear();
    cellArea.clear();
    sliceLookup.clear();
    thermalArea = axialFriction = azimuthalFriction = 0.0;
    const int stride[3] = { 1, air.nx, air.nx * air.ny };
    const auto isAir = [&air] (int cell)
    {
        const CellType t = air.type (cell);
        return t == CellType::Regular || t == CellType::Irregular || t == CellType::Slice;
    };
    std::map<long long, int> faceOf; // (cell on the - side) * 3 + axis -> face
    const auto faceIndex = [&] (int low, int axis)
    {
        const long long key = static_cast<long long> (low) * 3 + axis;
        auto found = faceOf.find (key);
        if (found == faceOf.end())
        {
            found = faceOf.emplace (key, static_cast<int> (faces.size())).first;
            faces.push_back ({ low, low + stride[axis] });
        }
        return found->second;
    };
    // Wall area per wall cell, and the axial (z face) friction areas.
    std::map<int, double> wallOf;
    std::map<int, double> axialOf; // face -> area
    for (int kz = 1; kz < air.nz - 1; ++kz)
        for (int j = 1; j < air.ny - 1; ++j)
            for (int i = 1; i < air.nx - 1; ++i)
            {
                const int cell = air.index (i, j, kz);
                if (! isAir (cell) || std::hypot (air.x (i), air.y (j)) >= radius || air.z (kz) >= 0.0 || air.z (kz) <= -depth)
                    continue;
                const std::uint8_t mask = air.faces (cell);
                for (int d = 0; d < 2; ++d)
                    for (int side = -1; side <= 1; side += 2)
                    {
                        if (air.type (cell + side * stride[d]) != CellType::Solid)
                            continue;
                        // A face on the barrel: its true area from the wall's normal at the face centre.
                        const double xf = air.x (i) + (d == 0 ? 0.5 * side * delta : 0.0);
                        const double yf = air.y (j) + (d == 1 ? 0.5 * side * delta : 0.0);
                        const double area = delta * delta * std::hypot (xf, yf) / (std::abs (xf) + std::abs (yf));
                        wallOf[cell] += area;
                        thermalArea += area;
                        for (int s = 0; s < 2; ++s)
                            if (mask & (1u << (4 + s))) // the z faces (not the ends or a head)
                                axialOf[faceIndex (s == 0 ? cell - stride[2] : cell, 2)] += 0.5 * area;
                    }
            }
    // The staircase's extent is quantised (in each quadrant the wall spans
    // round(a / Delta) rows and columns), so its faces add up to ~1 % less
    // than the side wall (0.9902 of 2 pi a L at a / Delta = 31.25): a uniform
    // factor makes them carry exactly 2 pi a L, the area the continuum loss
    // integrates over (numerical, like the heads' inclusion radius).
    staircaseArea = thermalArea;
    const double normalise = thermalArea > 0.0 ? 2.0 * 3.14159265358979323846 * radius * depth / thermalArea : 1.0;
    thermalArea *= normalise;
    for (auto& [cell, area] : wallOf)
        area *= normalise;
    for (const auto& [face, area] : axialOf)
    {
        Element e;
        e.begin = static_cast<int> (elementFaces.size());
        elementFaces.push_back (face);
        elementWeights.push_back (1.0);
        e.end = static_cast<int> (elementFaces.size());
        e.area = area * normalise;
        e.axial = true;
        axialFriction += e.area;
        elements.push_back (e);
    }
    // Azimuthal elements: the tangential part of each wall cell's centre
    // velocity, the mean flux of its two x faces and of its two y faces with
    // a wall face counting as zero flux. That mean interpolates the normal
    // flow, which vanishes at the wall, to the centre half a cell away, so the
    // projection on the tangent leaves out the normal flow; the tangential
    // flow of an irrotational field is flat across the wall (d v_t / dn =
    // d v_n / dt = 0 on a flat wall), so the centre carries its wall value.
    for (const auto& [cell, area] : wallOf)
    {
        const std::uint8_t mask = air.faces (cell);
        const int i = cell % air.nx, j = (cell / air.nx) % air.ny;
        const double theta = std::atan2 (air.y (j), air.x (i));
        const double along[2] = { -std::sin (theta), std::cos (theta) };
        Element e;
        e.begin = static_cast<int> (elementFaces.size());
        for (int d = 0; d < 2; ++d)
            for (int s = 0; s < 2; ++s)
                if (mask & (1u << (2 * d + s)))
                {
                    elementFaces.push_back (faceIndex (s == 0 ? cell - stride[d] : cell, d));
                    elementWeights.push_back (0.5 * along[d]);
                }
        e.end = static_cast<int> (elementFaces.size());
        e.area = area;
        if (e.end > e.begin)
        {
            azimuthalFriction += e.area;
            elements.push_back (e);
        }
    }
    // Each face's elements (transpose), and every cell a friction face touches
    // or with a wall face, ordinary ones first.
    faceStart.assign (faces.size() + 1, 0);
    for (const int f : elementFaces)
        ++faceStart[static_cast<std::size_t> (f) + 1];
    for (std::size_t f = 0; f < faces.size(); ++f)
        faceStart[f + 1] += faceStart[f];
    faceElements.assign (elementFaces.size(), 0);
    faceWeights.assign (elementFaces.size(), 0.0);
    {
        std::vector<int> fill (faceStart.begin(), faceStart.end() - 1);
        for (std::size_t e = 0; e < elements.size(); ++e)
            for (int q = elements[e].begin; q < elements[e].end; ++q)
            {
                const int slot = fill[static_cast<std::size_t> (elementFaces[static_cast<std::size_t> (q)])]++;
                faceElements[static_cast<std::size_t> (slot)] = static_cast<int> (e);
                faceWeights[static_cast<std::size_t> (slot)] = elementWeights[static_cast<std::size_t> (q)];
            }
    }
    std::map<int, std::vector<int>> touching; // cell -> 2 * face + side
    for (const auto& [cell, area] : wallOf)
        touching[cell];
    for (std::size_t f = 0; f < faces.size(); ++f)
    {
        touching[faces[f].low].push_back (static_cast<int> (2 * f));
        touching[faces[f].high].push_back (static_cast<int> (2 * f + 1));
    }
    wallCells = static_cast<int> (wallOf.size());
    for (const bool slice : { false, true })
        for (const auto& [cell, list] : touching)
        {
            if ((air.type (cell) == CellType::Slice) != slice)
                continue;
            Cell entry;
            entry.cell = cell;
            entry.faceBegin = static_cast<int> (cellFaces.size());
            cellFaces.insert (cellFaces.end(), list.begin(), list.end());
            entry.faceEnd = static_cast<int> (cellFaces.size());
            if (const auto w = wallOf.find (cell); w != wallOf.end())
            {
                entry.thermal = static_cast<int> (cellArea.size());
                cellArea.push_back (w->second);
            }
            if (slice)
                sliceLookup.push_back ({ cell, static_cast<int> (cells.size()) });
            cells.push_back (entry);
        }
    ordinaryCells = static_cast<int> (cells.size() - sliceLookup.size());
    std::sort (sliceLookup.begin(), sliceLookup.end());
    const std::size_t poles = xi.size();
    for (int level = 0; level < 2; ++level)
    {
        theta[level].assign (cellArea.size() * poles, 0.0);
        phi[level].assign (elements.size() * poles, 0.0);
        drag[level].assign (elements.size(), 0.0);
    }
    decay.resize (poles);
    gain.resize (poles);
    relief.resize (poles);
    for (std::size_t i = 0; i < poles; ++i)
    {
        decay[i] = (1.0 - k * xi[i]) / (1.0 + k * xi[i]);
        gain[i] = 2.0 * k / (1.0 + k * xi[i]);
        relief[i] = 1.0 / (1.0 + k * xi[i]);
    }
    setConstants();
}

template <typename Real>
void WallLoss<Real>::setConstants()
{
    // beta = (c^2 k rho / (2 Delta^3)) A C_th sum_i mu_i / (1 + k xi_i): the
    // thermal branches' instantaneous conductance in the implicit cell update.
    double direct = 0.0;
    for (std::size_t i = 0; i < xi.size(); ++i)
        direct += muThermal[i] * relief[i];
    beta.resize (cellArea.size());
    for (std::size_t t = 0; t < cellArea.size(); ++t)
        beta[t] = c * c * k * rho / (2.0 * delta * delta * delta) * thermalScale * cTh * cellArea[t] * direct;
    faceGain.resize (elements.size());
    for (std::size_t e = 0; e < elements.size(); ++e)
        faceGain[e] = viscousScale * rootNu * elements[e].area / (delta * delta * delta);
}

template <typename Real>
void WallLoss<Real>::scale (double thermal, double viscous)
{
    thermalScale = thermal;
    viscousScale = viscous;
    setConstants();
}

template <typename Real>
int WallLoss<Real>::sliceEntry (int cell) const
{
    const auto it = std::lower_bound (sliceLookup.begin(), sliceLookup.end(), std::pair<int, int> (cell, -1));
    return it != sliceLookup.end() && it->first == cell ? it->second : -1;
}

template class WallLoss<float>;
template class WallLoss<double>;

} // namespace taikor::candidates::fdtd
