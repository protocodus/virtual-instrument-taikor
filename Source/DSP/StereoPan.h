#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace taikor
{
// Preserve the microphone pair at center and fold both microphones into the
// destination at either edge. The pan law preserves power for a mono source.
struct StereoPan
{
    // Double state lets slow smoothing reach its endpoint at high sample rates.
    double ll = 1.0, lr = 0.0, rl = 0.0, rr = 1.0;

    static StereoPan atPosition (float position) noexcept
    {
        if (! std::isfinite (position))
            return {};

        const float distance = std::clamp (std::abs (position), 0.0f, 1.0f);
        if (distance == 0.0f)
            return {};
        const double angle = distance * 0.7853981633974483;
        const double c = distance == 1.0f ? 0.7071067811865475 : std::cos (angle);
        const double s = distance == 1.0f ? c : std::sin (angle);
        return position < 0.0f ? StereoPan { c, s, 0.0, c - s }
                              : StereoPan { c - s, 0.0, s, c };
    }

    void apply (float& left, float& right) const noexcept
    {
        if (isCentered())
            return;

        if (! std::isfinite (left) || ! std::isfinite (right)
            || ! std::isfinite (ll) || ! std::isfinite (lr)
            || ! std::isfinite (rl) || ! std::isfinite (rr))
        {
            left = right = 0.0f;
            return;
        }

        const double originalLeft = left;
        const double originalRight = right;
        const double nextLeft = ll * originalLeft + lr * originalRight;
        const double nextRight = rl * originalLeft + rr * originalRight;
        constexpr double maximum = std::numeric_limits<float>::max();
        left = static_cast<float> (std::clamp (nextLeft, -maximum, maximum));
        right = static_cast<float> (std::clamp (nextRight, -maximum, maximum));
    }

    void approach (const StereoPan& target, double smoothing) noexcept
    {
        if (! std::isfinite (smoothing))
            return;
        smoothing = std::clamp (smoothing, 0.0, 1.0);
        const auto step = [smoothing] (double& current, double next)
        {
            if (! std::isfinite (next))
                return;
            if (! std::isfinite (current))
                current = next;
            current += smoothing * (next - current);
            if (std::abs (current - next) < 1.0e-10)
                current = next;
        };
        step (ll, target.ll);
        step (lr, target.lr);
        step (rl, target.rl);
        step (rr, target.rr);
    }

    bool isCentered() const noexcept
    { return ll == 1.0f && lr == 0.0f && rl == 0.0f && rr == 1.0f; }
};
} // namespace taikor
