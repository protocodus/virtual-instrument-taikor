#include "DSP/EnsembleEngine.h"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

namespace
{
using Clock = std::chrono::steady_clock;
struct Result
{
    double wallMs = 0.0, cpuMs = 0.0, p95Ms = 0.0, worstMs = 0.0;
    std::vector<float> left, right;
};

Result render (taikor::EnsembleEngine& engine, double rate, int block)
{
    const int samples = static_cast<int> (rate * 0.4);
    Result result;
    result.left.resize (static_cast<std::size_t> (samples));
    result.right.resize (static_cast<std::size_t> (samples));
    std::vector<double> durations;
    durations.reserve (static_cast<std::size_t> ((samples + block - 1) / block));
    engine.reset();
    for (int drum = 0; drum < 4; ++drum)
        engine.trigger (taikor::Articulation::Don, drum, 0.85f);
    const auto cpuStart = std::clock();
    const auto start = Clock::now();
    for (int offset = 0; offset < samples; offset += block)
    {
        const auto blockStart = Clock::now();
        engine.process (result.left.data() + offset, result.right.data() + offset,
                        std::min (block, samples - offset));
        durations.push_back (std::chrono::duration<double, std::milli> (Clock::now() - blockStart).count());
    }
    result.wallMs = std::chrono::duration<double, std::milli> (Clock::now() - start).count();
    result.cpuMs = 1000.0 * static_cast<double> (std::clock() - cpuStart) / CLOCKS_PER_SEC;
    std::sort (durations.begin(), durations.end());
    result.p95Ms = durations[static_cast<std::size_t> (0.95 * static_cast<double> (durations.size() - 1))];
    result.worstMs = durations.back();
    return result;
}

bool same (const Result& a, const Result& b)
{
    for (std::size_t sample = 0; sample < a.left.size(); ++sample)
        if (std::bit_cast<std::uint32_t> (a.left[sample]) != std::bit_cast<std::uint32_t> (b.left[sample])
            || std::bit_cast<std::uint32_t> (a.right[sample]) != std::bit_cast<std::uint32_t> (b.right[sample])
            || ! std::isfinite (b.left[sample]) || ! std::isfinite (b.right[sample]))
            return false;
    return true;
}

double median (std::array<double, 3> values)
{
    std::sort (values.begin(), values.end());
    return values[1];
}
} // namespace

int main()
{
    std::cout << "OFFLINE benchmark: wall time is bounce speed; process CPU time includes every worker.\n";
    std::cout << "players,rate,block,workers,serial_wall_ms,parallel_wall_ms,wall_reduction_percent,"
                 "serial_cpu_ms,parallel_cpu_ms,serial_p95_ms,parallel_p95_ms,parallel_worst_ms,bit_exact\n";
    for (const int size : { 4, 8 })
        for (const auto configuration : { std::pair { 48000.0, 256 }, std::pair { 96000.0, 64 } })
        {
            const auto [rate, block] = configuration;
            auto serial = std::make_unique<taikor::EnsembleEngine>();
            auto parallel = std::make_unique<taikor::EnsembleEngine>();
            taikor::EngineParameters parameters;
            parameters.ensembleSize = size;
            parameters.ensembleVariation = 0.4f;
            for (auto* engine : { serial.get(), parallel.get() })
            {
                engine->setParameters (parameters);
                engine->prepare (rate, block);
            }
            parallel->prepareOfflineRendering (3);
            parallel->setOfflineRendering (true);
            (void) render (*serial, rate, block);
            (void) render (*parallel, rate, block);
            std::array<double, 3> serialWall {}, parallelWall {}, serialCpu {}, parallelCpu {}, serialP95 {}, parallelP95 {};
            double worst = 0.0;
            for (int round = 0; round < 3; ++round)
            {
                Result before, after;
                if (round % 2 == 0)
                {
                    before = render (*serial, rate, block);
                    after = render (*parallel, rate, block);
                }
                else
                {
                    after = render (*parallel, rate, block);
                    before = render (*serial, rate, block);
                }
                if (! same (before, after))
                {
                    std::cerr << "FAIL: offline parallel audio differs from serial\n";
                    return 1;
                }
                serialWall[round] = before.wallMs;
                parallelWall[round] = after.wallMs;
                serialCpu[round] = before.cpuMs;
                parallelCpu[round] = after.cpuMs;
                serialP95[round] = before.p95Ms;
                parallelP95[round] = after.p95Ms;
                worst = std::max (worst, after.worstMs);
            }
            std::cout << std::fixed << std::setprecision (4) << size << ',' << rate << ',' << block << ','
                      << parallel->getOfflineWorkerCount() << ',' << median (serialWall) << ',' << median (parallelWall) << ','
                      << 100.0 * (1.0 - median (parallelWall) / median (serialWall)) << ','
                      << median (serialCpu) << ',' << median (parallelCpu) << ','
                      << median (serialP95) << ',' << median (parallelP95) << ',' << worst << ",true\n";
        }
}
