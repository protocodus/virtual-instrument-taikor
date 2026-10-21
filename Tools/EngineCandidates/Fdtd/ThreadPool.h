#pragma once

// A persistent worker pool with a reusable barrier, for the finite-difference
// drum's time step: every thread runs the same step loop and meets the others
// at a few barriers per step (FdtdDrum.cpp). The machine this renders on is
// shared, so waiting threads spin briefly and then sleep on the atomic
// (C++20 atomic wait), rather than spinning without bound while a descheduled
// peer holds everyone up.

#include <atomic>
#include <cstdlib>
#include <functional>
#include <string>
#include <thread>
#include <vector>

namespace taikor::candidates::fdtd
{

class Barrier
{
public:
    explicit Barrier (int count) : parties (count) {}

    void wait() noexcept
    {
        const unsigned gen = generation.load (std::memory_order_acquire);
        if (arrived.fetch_add (1, std::memory_order_acq_rel) == parties - 1)
        {
            arrived.store (0, std::memory_order_relaxed);
            generation.fetch_add (1, std::memory_order_acq_rel);
            generation.notify_all();
            return;
        }
        for (int spin = 0; spin < spinLimit; ++spin)
            if (generation.load (std::memory_order_acquire) != gen)
                return;
        while (generation.load (std::memory_order_acquire) == gen)
            generation.wait (gen, std::memory_order_acquire);
    }

    [[nodiscard]] int size() const noexcept { return parties; }

private:
    // Numerical/performance only: about 10-20 us of spinning on an M1 core
    // before the thread sleeps.
    static constexpr int spinLimit = 20000;
    int parties;
    std::atomic<int> arrived { 0 };
    std::atomic<unsigned> generation { 0 };
};

// Runs job(thread) on `count` threads (the caller is thread 0) and returns
// when all have finished.
class ThreadPool
{
public:
    explicit ThreadPool (int count) : threads (count < 1 ? 1 : count), start (threads), finish (threads)
    {
        for (int t = 1; t < threads; ++t)
            workers.emplace_back ([this, t] { loop (t); });
    }

    ~ThreadPool()
    {
        stopping = true;
        start.wait();
        for (auto& worker : workers)
            worker.join();
    }

    ThreadPool (const ThreadPool&) = delete;
    ThreadPool& operator= (const ThreadPool&) = delete;

    [[nodiscard]] int size() const noexcept { return threads; }

    void run (const std::function<void (int)>& job)
    {
        current = &job;
        start.wait();
        job (0);
        finish.wait();
        current = nullptr;
    }

private:
    void loop (int t)
    {
        for (;;)
        {
            start.wait();
            if (stopping)
                return;
            (*current) (t);
            finish.wait();
        }
    }

    int threads;
    Barrier start, finish;
    std::vector<std::thread> workers;
    const std::function<void (int)>* current = nullptr;
    std::atomic<bool> stopping { false };
};

// Thread count for the candidate renders: TAIKOR_CANDIDATE_THREADS, default 8.
[[nodiscard]] inline int candidateThreadCount()
{
    if (const char* text = std::getenv ("TAIKOR_CANDIDATE_THREADS"))
    {
        const int value = std::atoi (text);
        if (value >= 1 && value <= 64)
            return value;
    }
    return 8;
}

// Splits [0, total) into `parts` contiguous ranges; returns range `index`.
inline void splitRange (int total, int parts, int index, int& begin, int& end) noexcept
{
    begin = static_cast<int> ((static_cast<long long> (total) * index) / parts);
    end = static_cast<int> ((static_cast<long long> (total) * (index + 1)) / parts);
}

} // namespace taikor::candidates::fdtd
