#include "DSP/PackedModalBank.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <utility>

namespace
{
constexpr std::size_t count = 17;
struct Mode
{
    taikor::detail::PackedModalResonator resonator;
    std::uint16_t physicalIndex {};
};
using Bank = taikor::detail::PackedModalBank<Mode, count>;

bool same (double a, double b)
{
    return std::bit_cast<std::uint64_t> (a) == std::bit_cast<std::uint64_t> (b);
}

bool exercise (Bank& bank, int active, int steps)
{
    std::array<taikor::detail::ScalarModalResonator, count> reference {};
    std::array<float, count> inputs {};
    for (std::size_t i = 0; i < count; ++i)
    {
        const auto& r = bank[i].resonator;
        reference[i] = { r.a1, r.a2, r.b0, r.y1, r.y2 };
    }
    for (int step = 0; step < steps; ++step)
    {
        for (std::size_t i = 0; i < count; ++i)
            inputs[i] = static_cast<float> ((step * 13 + static_cast<int> (i) * 7) % 31 - 15) * 0.0001f;
        for (int i = 0; i < active; ++i)
            (void) reference[static_cast<std::size_t> (i)].tick (
                inputs[bank[static_cast<std::size_t> (i)].physicalIndex]);
        bank.tick (inputs, active);
        for (std::size_t i = 0; i < count; ++i)
        {
            const auto& r = bank[i].resonator;
            if (! same (r.y1, reference[i].y1) || ! same (r.y2, reference[i].y2))
                return false;
        }
    }
    return true;
}

bool check()
{
    Bank original;
    for (std::size_t i = 0; i < count; ++i)
    {
        auto& r = original[i].resonator;
        original[i].physicalIndex = static_cast<std::uint16_t> (count - 1 - i);
        r.a1 = -1.89 + static_cast<double> (i) * 0.013;
        r.a2 = 0.991 - static_cast<double> (i) * 0.001;
        r.b0 = 0.025 + static_cast<double> (i) * 0.001;
        r.y1 = 0.0007 * static_cast<double> (i + 1);
        r.y2 = -0.0003 * static_cast<double> (i + 1);
    }

    // Copies own independent state. The mutable metadata array can be sorted
    // using normal value semantics without rebinding its slot-owned states.
    for (int active = 0; active <= static_cast<int> (count); ++active)
    {
        Bank bank (original);
        if (&bank[0].resonator.y1 == &original[0].resonator.y1)
            return false;
        const auto untouched = original[0].resonator.y1;
        if (! exercise (bank, active, 4096) || ! same (original[0].resonator.y1, untouched))
            return false;

        bank = original;
        std::sort (bank.begin(), bank.end(), [] (const Mode& a, const Mode& b)
        {
            return a.physicalIndex < b.physicalIndex;
        });
        for (std::size_t i = 0; i < count; ++i)
            if (bank[i].physicalIndex != i
                || ! same (bank[i].resonator.y1, original[count - 1 - i].resonator.y1))
                return false;
        if (! exercise (bank, active, 128))
            return false;

        Bank moved (std::move (bank));
        if (&moved[0].resonator.y1 == &bank[0].resonator.y1
            || ! exercise (moved, active, 128))
            return false;

        Mode saved = original[0];
        const double savedState = saved.resonator.y1;
        moved[0] = saved;
        moved[0].resonator.y1 += 1.0;
        if (! same (saved.resonator.y1, savedState))
            return false;
        moved[0].resonator.clear();
        moved[0].resonator.a1 = -1.75;
        moved[0].resonator.b0 = 0.18;
        if (! exercise (moved, active, 128))
            return false;
    }
    return true;
}
} // namespace

int main()
{
    if (! check())
    {
        std::cerr << "FAIL: packed modal state ownership or recurrence changed\n";
        return 1;
    }
    std::cout << "Packed modal bank preserves scalar state across all 18 active prefixes, "
                 "long recurrence runs, sorting, copies, moves and direct state edits\n";
    return 0;
}
