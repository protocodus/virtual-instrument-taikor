// The root scan of the clamped stiff membrane at high orders: no phantom roots
// where J_m underflows (x << m), every mode with b1 a > m and a positive norm,
// and the per-order count equal to an independent finer scan.

#include "ModalMath.h"
#include "ModalTestUtil.h"

#include <cmath>

using namespace taikor::candidates;
using namespace taikor::candidates::modal;

MODAL_TEST (heads_root_scan)
{
    HeadMaterial h;
    h.tension = 6146.67;
    h.arealDensity = 0.94;
    h.bending = 0.266213;
    h.radius = 0.39;
    const double omegaMax = 2.0 * pi * 21600.0;
    const auto modes = solveHeadModes (h, omegaMax);
    int total = 0, belowOrder = 0, badNorm = 0, badIndex = 0, countMismatch = 0;
    for (std::size_t m = 0; m < modes.size(); ++m)
    {
        int expectedN = 1;
        for (const auto& mode : modes[m])
        {
            ++total;
            if (mode.b1 * h.radius <= static_cast<double> (m))
                ++belowOrder;
            if (! (mode.normSq > 0.0))
                ++badNorm;
            if (mode.n != expectedN++)
                ++badIndex;
        }
        // Independent count: sign changes of b1 J_m'(x) - b2 J_m(x) I_m'(b2 a)/I_m(b2 a)
        // on a grid five times finer, from x = m to the band edge.
        const auto f = [&] (double x)
        {
            const double b1 = x / h.radius, b2 = std::sqrt (b1 * b1 + h.tension / h.bending);
            const double y = b2 * h.radius;
            const int order = static_cast<int> (m);
            const double ratio = 0.5 * (special::scaledBesselI (order == 0 ? 1 : order - 1, y)
                                        + special::scaledBesselI (order + 1, y)) / special::scaledBesselI (order, y);
            return b1 * special::besselJPrime (order, x) - b2 * special::besselJ (order, x) * ratio;
        };
        const double b1Max = std::sqrt ((std::sqrt (h.tension * h.tension + 4.0 * h.bending * h.arealDensity * omegaMax * omegaMax)
                                         - h.tension) / (2.0 * h.bending));
        int count = 0;
        double previous = f (std::max (static_cast<double> (m), 0.01));
        for (double x = std::max (static_cast<double> (m), 0.01) + 0.01; x <= b1Max * h.radius; x += 0.01)
        {
            const double value = f (x);
            if ((value < 0.0) != (previous < 0.0))
                ++count;
            previous = value;
        }
        if (count != static_cast<int> (modes[m].size()))
            ++countMismatch;
    }
    modaltest::note ("batter modes to 21.6 kHz: %.0f", total);
    modaltest::note ("orders: %.0f", static_cast<double> (modes.size()));
    modaltest::below ("modes with b1 a <= m", belowOrder, 0.0);
    modaltest::below ("modes with non-positive norm", badNorm, 0.0);
    modaltest::below ("modes out of radial sequence", badIndex, 0.0);
    modaltest::below ("orders whose count differs from the finer scan", countMismatch, 0.0);
}
