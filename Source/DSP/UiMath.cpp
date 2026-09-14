#include "DSP/UiMath.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace taikor::ui
{
namespace
{
constexpr float minimumLinear = 1.0e-6f;
} // namespace

float clamp (float value, float low, float high) noexcept
{
    if (std::isnan (low) || std::isnan (high))
        return 0.0f;
    if (low > high)
        std::swap (low, high);
    if (! (value == value))
        return low;
    return value < low ? low : (value > high ? high : value);
}

float onePoleCoefficient (float seconds, float updateRateHz) noexcept
{
    if (! (updateRateHz > 0.0f))
        return 1.0f;
    if (! (seconds > 0.0f))
        return 1.0f;
    const float updates = seconds * updateRateHz;
    if (! std::isfinite (updates))
        return 0.0f;
    return clamp (-std::expm1 (-1.0f / updates), 0.0f, 1.0f);
}

float decayMultiplier (float decibels, float seconds, float updateRateHz) noexcept
{
    if (! (updateRateHz > 0.0f) || ! (seconds > 0.0f))
        return 0.0f;
    if (! std::isfinite (decibels))
        return 0.0f;
    const float updates = seconds * updateRateHz;
    if (updates == 0.0f)
        return decibels < 0.0f ? 0.0f : std::numeric_limits<float>::max();
    const float result = std::pow (10.0f, decibels / (20.0f * updates));
    if (std::isfinite (result))
        return result;
    return result > 0.0f ? std::numeric_limits<float>::max() : 0.0f;
}

float meterPositionForLinear (float linear, float floorDecibels) noexcept
{
    if (! std::isfinite (floorDecibels) || ! (floorDecibels < 0.0f))
        return 0.0f;
    const float bounded = linear > minimumLinear ? linear : minimumLinear;
    const float decibels = 20.0f * std::log10 (bounded);
    return clamp ((decibels - floorDecibels) / -floorDecibels, 0.0f, 1.0f);
}

float linearForMeterPosition (float position, float floorDecibels) noexcept
{
    if (! std::isfinite (floorDecibels) || ! (floorDecibels < 0.0f))
        return 0.0f;
    const float bounded = clamp (position, 0.0f, 1.0f);
    const float decibels = floorDecibels + bounded * -floorDecibels;
    return std::pow (10.0f, decibels / 20.0f);
}

void MeterBallistics::reset() noexcept
{
    level = 0.0f;
    peak = 0.0f;
    holdCountdown = 0.0f;
}

void MeterBallistics::update (float target, float attackCoefficient,
                              float releaseCoefficient, float peakFall,
                              float holdUpdates) noexcept
{
    if (! std::isfinite (level) || level < 0.0f)
        level = 0.0f;
    if (! std::isfinite (peak) || peak < 0.0f)
        peak = 0.0f;
    if (! std::isfinite (holdCountdown) || holdCountdown < 0.0f)
        holdCountdown = 0.0f;

    const float bounded = std::isfinite (target) && target > 0.0f ? target : 0.0f;
    const float coefficient = bounded > level ? attackCoefficient : releaseCoefficient;
    level += clamp (coefficient, 0.0f, 1.0f) * (bounded - level);

    if (level >= peak)
    {
        peak = level;
        holdCountdown = std::isfinite (holdUpdates) && holdUpdates > 0.0f
            ? holdUpdates : 0.0f;
    }
    else if (holdCountdown > 0.0f)
    {
        holdCountdown -= 1.0f;
    }
    else
    {
        peak *= clamp (peakFall, 0.0f, 1.0f);
        if (peak < level)
            peak = level;
    }
}

RowLayout rowLayout (int extent, int columns, int gap, int count) noexcept
{
    RowLayout layout;
    if (columns <= 0 || extent <= 0)
        return layout;

    const std::int64_t safeGap = std::max<std::int64_t> (gap, 0);
    const std::int64_t totalGap = safeGap * (static_cast<std::int64_t> (columns) - 1);
    const std::int64_t available = static_cast<std::int64_t> (extent) - totalGap;
    const std::int64_t cellSize = available > columns ? available / columns : 1;
    layout.cellSize = static_cast<int> (std::clamp<std::int64_t> (
        cellSize, 1, std::numeric_limits<int>::max()));

    const std::int64_t used = count > 0 && count <= columns
        ? static_cast<std::int64_t> (count) * layout.cellSize
              + safeGap * (static_cast<std::int64_t> (count) - 1)
        : static_cast<std::int64_t> (columns) * layout.cellSize + totalGap;
    const std::int64_t origin = (static_cast<std::int64_t> (extent) - used) / 2;
    layout.origin = static_cast<int> (std::clamp<std::int64_t> (
        origin, 0, std::numeric_limits<int>::max()));

    return layout;
}

int cellOffset (const RowLayout& layout, int gap, int index) noexcept
{
    if (index <= 0)
        return layout.origin;
    const std::int64_t safeGap = std::max<std::int64_t> (gap, 0);
    const std::int64_t offset = static_cast<std::int64_t> (layout.origin)
        + static_cast<std::int64_t> (index)
            * (static_cast<std::int64_t> (std::max (layout.cellSize, 1)) + safeGap);
    return static_cast<int> (std::clamp<std::int64_t> (
        offset, 0, std::numeric_limits<int>::max()));
}

HeadPoint headPointFor (float normalisedRadius, float angleRadians) noexcept
{
    if (! std::isfinite (angleRadians))
        return {};
    const float radius = clamp (normalisedRadius, 0.0f, 1.0f);
    HeadPoint point;
    point.x = radius * std::cos (angleRadians);
    // Screen coordinates run downwards, so the sine is negated once here
    // rather than at every call site.
    point.y = -radius * std::sin (angleRadians);
    return point;
}

float semitonesBetween (float frequencyHz, float referenceHz) noexcept
{
    if (! std::isfinite (frequencyHz) || ! std::isfinite (referenceHz)
        || ! (frequencyHz > 0.0f) || ! (referenceHz > 0.0f))
        return 0.0f;
    return 12.0f * (std::log2 (frequencyHz) - std::log2 (referenceHz));
}

float mix (float a, float b, float amount) noexcept
{
    if (! std::isfinite (a))
        a = 0.0f;
    if (! std::isfinite (b))
        b = 0.0f;
    const float boundedAmount = clamp (amount, 0.0f, 1.0f);
    const double result = static_cast<double> (a)
        + (static_cast<double> (b) - a) * boundedAmount;
    constexpr double maximum = std::numeric_limits<float>::max();
    return static_cast<float> (std::clamp (result, -maximum, maximum));
}

float smoothStep (float edge0, float edge1, float value) noexcept
{
    if (! std::isfinite (edge0) || ! std::isfinite (edge1))
        return 0.0f;
    if (edge0 == edge1)
        return value < edge0 ? 0.0f : 1.0f;
    if (! std::isfinite (value))
        return value > 0.0f ? 1.0f : 0.0f;
    const double t = std::clamp ((static_cast<double> (value) - edge0)
                                     / (static_cast<double> (edge1) - edge0),
                                 0.0, 1.0);
    return t * t * (3.0f - 2.0f * t);
}

} // namespace taikor::ui
