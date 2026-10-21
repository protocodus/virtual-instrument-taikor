#pragma once

// A team of threads evaluating one von Karman potential together. The drum
// (Fdtd/FdtdEngineStep.cpp) hands the evaluate of its batter to `size`
// threads of its pool, member 0 to size - 1, each calling
// VonKarmanPotential::evaluateMember with its own Team; they split every pass
// over the grid by rows and meet at the team's barrier where a pass reads
// rows another member wrote. Reductions sum the members' partials in member
// order, so a take is reproducible for a given team size (a different size
// rounds differently). A team of one is the serial evaluate, bit for bit.

#include "Fdtd/ThreadPool.h"

#include <algorithm>
#include <array>

namespace taikor::candidates::vk
{

constexpr int maxTeam = 4;

// Reduction scratch shared by a team: kSlots rotating rows of maxTeam values
// (a barrier separates the reading of one slot from its next writing whenever
// kSlots >= 2; see Team::sum).
struct TeamShared
{
    static constexpr int kSlots = 4;
    std::array<std::array<double, maxTeam * 2>, kSlots> slots {};
};

struct Team
{
    int member = 0, size = 1;
    fdtd::Barrier* barrier = nullptr;
    TeamShared* shared = nullptr;
    int next = 0; // this member's reduction count (every member counts alike)

    [[nodiscard]] static Team serial() noexcept { return {}; }

    void sync() noexcept
    {
        if (size > 1)
            barrier->wait();
    }

    // This member's share [begin, end) of `total` rows.
    void rows (int total, int& begin, int& end) const noexcept
    {
        fdtd::splitRange (total, size, member, begin, end);
    }

    // Sums of a and b over the members (a barrier, and a view of every
    // member's writes before it).
    std::array<double, 2> sum2 (double a, double b) noexcept
    {
        if (size == 1)
            return { a, b };
        auto& slot = shared->slots[static_cast<std::size_t> (next++ % TeamShared::kSlots)];
        slot[static_cast<std::size_t> (2 * member)] = a;
        slot[static_cast<std::size_t> (2 * member + 1)] = b;
        barrier->wait();
        std::array<double, 2> total { 0.0, 0.0 };
        for (int m = 0; m < size; ++m)
        {
            total[0] += slot[static_cast<std::size_t> (2 * m)];
            total[1] += slot[static_cast<std::size_t> (2 * m + 1)];
        }
        return total;
    }

    double sum (double a) noexcept { return sum2 (a, 0.0)[0]; }

    double max (double a) noexcept
    {
        if (size == 1)
            return a;
        auto& slot = shared->slots[static_cast<std::size_t> (next++ % TeamShared::kSlots)];
        slot[static_cast<std::size_t> (2 * member)] = a;
        barrier->wait();
        double best = slot[0];
        for (int m = 1; m < size; ++m)
            best = std::max (best, slot[static_cast<std::size_t> (2 * m)]);
        return best;
    }
};

} // namespace taikor::candidates::vk
