#pragma once

#include <array>
#include <cstddef>
#include <memory>

#if defined(__APPLE__) && defined(__aarch64__) && defined(__clang__) \
    && ! defined(TAIKOR_DISABLE_PACKED_MODAL_BANK)
 #define TAIKOR_USE_PACKED_MODAL_BANK 1
 #include <arm_neon.h>
#else
 #define TAIKOR_USE_PACKED_MODAL_BANK 0
#endif

namespace taikor::detail
{

struct ScalarModalResonator
{
    double a1 { 0.0 };
    double a2 { 0.0 };
    double b0 { 0.0 };
    double y1 { 0.0 };
    double y2 { 0.0 };

    void clear() noexcept { y1 = y2 = 0.0; }
    [[nodiscard]] float tick (float input) noexcept
    {
        const double output = b0 * static_cast<double> (input) - a1 * y1 - a2 * y2;
        y2 = y1;
        y1 = output;
        return static_cast<float> (output);
    }
};

// A normal value when copied out of a bank (for example by insertion sort),
// and a view into five packed arrays while it occupies a bank slot. Assignment
// copies VALUES, never bindings, so moving metadata cannot make two modes own
// the same state or leave a reference pointing into a temporary sort key.
class PackedModalResonator
{
private:
    ScalarModalResonator local {};

public:
    double& a1;
    double& a2;
    double& b0;
    double& y1;
    double& y2;

    PackedModalResonator() noexcept
        : a1 (local.a1), a2 (local.a2), b0 (local.b0),
          y1 (local.y1), y2 (local.y2) {}

    PackedModalResonator (double& first, double& second, double& input,
                          double& current, double& previous) noexcept
        : a1 (first), a2 (second), b0 (input), y1 (current), y2 (previous) {}

    PackedModalResonator (const PackedModalResonator& other) noexcept
        : PackedModalResonator() { *this = other; }

    PackedModalResonator& operator= (const PackedModalResonator& other) noexcept
    {
        if (this != &other)
        {
            a1 = other.a1; a2 = other.a2; b0 = other.b0;
            y1 = other.y1; y2 = other.y2;
        }
        return *this;
    }

    void bind (double& first, double& second, double& input,
               double& current, double& previous) noexcept
    {
        const ScalarModalResonator saved { a1, a2, b0, y1, y2 };
        // These objects are complete, non-const member objects of Mode, not
        // base subobjects. Reconstruction rebinds the references once at bank
        // construction; no render-time allocation or state synchronization.
        std::destroy_at (this);
        std::construct_at (this, first, second, input, current, previous);
        a1 = saved.a1; a2 = saved.a2; b0 = saved.b0;
        y1 = saved.y1; y2 = saved.y2;
    }

    void clear() noexcept { y1 = y2 = 0.0; }
    [[nodiscard]] float tick (float input) noexcept
    {
        const double output = b0 * static_cast<double> (input) - a1 * y1 - a2 * y2;
        y2 = y1;
        y1 = output;
        return static_cast<float> (output);
    }
};

template <typename Mode, std::size_t Count>
class PackedModalBank
{
private:
    // Independent contiguous streams, not a second copy of the live state.
    alignas (64) std::array<double, Count> a1 {};
    alignas (64) std::array<double, Count> a2 {};
    alignas (64) std::array<double, Count> b0 {};
    alignas (64) std::array<double, Count> y1 {};
    alignas (64) std::array<double, Count> y2 {};
    std::array<Mode, Count> modes {};

public:
    PackedModalBank() noexcept
    {
        for (std::size_t i = 0; i < Count; ++i)
            modes[i].resonator.bind (a1[i], a2[i], b0[i], y1[i], y2[i]);
    }

    PackedModalBank (const PackedModalBank& other) noexcept : PackedModalBank()
    {
        modes = other.modes;
    }

    PackedModalBank& operator= (const PackedModalBank& other) noexcept
    {
        if (this != &other)
            modes = other.modes;
        return *this;
    }

    [[nodiscard]] Mode& operator[] (std::size_t i) noexcept { return modes[i]; }
    [[nodiscard]] const Mode& operator[] (std::size_t i) const noexcept { return modes[i]; }
    [[nodiscard]] auto begin() noexcept { return modes.begin(); }
    [[nodiscard]] auto end() noexcept { return modes.end(); }
    [[nodiscard]] auto begin() const noexcept { return modes.begin(); }
    [[nodiscard]] auto end() const noexcept { return modes.end(); }
    [[nodiscard]] Mode* data() noexcept { return modes.data(); }
    [[nodiscard]] const Mode* data() const noexcept { return modes.data(); }
    [[nodiscard]] constexpr std::size_t size() const noexcept { return Count; }

    void tick (const std::array<float, Count>& inputs, int activeCount) noexcept
    {
        const auto count = activeCount > 0
            ? (static_cast<std::size_t> (activeCount) < Count
                ? static_cast<std::size_t> (activeCount) : Count) : 0;
        std::size_t i = 0;
#if TAIKOR_USE_PACKED_MODAL_BANK
        for (; i + 1 < count; i += 2)
        {
            const float64x2_t input {
                static_cast<double> (inputs[modes[i].physicalIndex]),
                static_cast<double> (inputs[modes[i + 1].physicalIndex]) };
            const auto current = vld1q_f64 (y1.data() + i);
            const auto previous = vld1q_f64 (y2.data() + i);
            // Match the established Apple Clang scalar contraction exactly:
            // fnmul(a1, y1), fmadd(b0, input, -a1*y1), fmsub(a2, y2).
            // Do not replace it with a differently grouped fused expression.
            const auto negativeFeedback = vnegq_f64 (
                vmulq_f64 (vld1q_f64 (a1.data() + i), current));
            auto output = vfmaq_f64 (negativeFeedback,
                                    vld1q_f64 (b0.data() + i), input);
            output = vfmsq_f64 (output, vld1q_f64 (a2.data() + i), previous);
            vst1q_f64 (y2.data() + i, current);
            vst1q_f64 (y1.data() + i, output);
        }
#endif
        // Odd active prefixes and non-Apple targets keep the scalar operation.
        for (; i < count; ++i)
            (void) modes[i].resonator.tick (inputs[modes[i].physicalIndex]);
    }
};

} // namespace taikor::detail
