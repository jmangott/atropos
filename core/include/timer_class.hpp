#ifndef TIMER_CLASS_HPP
#define TIMER_CLASS_HPP

#include <cmath>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

// Adapted from `sldg` (Einkemmer 2020)

/// This timer class measures the elapsed time between two events. Timers can be
/// started and stopped repeatedly. The total time as well as the average time
/// between two events can be queried using the total() and average() methods,
/// respectively.
struct ntime
{
    std::chrono::time_point<std::chrono::high_resolution_clock> t_start;
    bool running;
    double elapsed;
    unsigned counter;
    double elapsed_sq;

    ntime()
    {
        counter = 0;
        elapsed = 0.0;
        running = false;
        elapsed_sq = 0.0;
    }

    void reset()
    {
        counter = 0;
        elapsed = 0.0;
        running = false;
        elapsed_sq = 0.0;
    }

    void start()
    {
        t_start = std::chrono::high_resolution_clock::now();
        // clock_gettime(CLOCK_MONOTONIC, &t_start);
        running = true;
    }

    /// The stop method returns the elapsed time since the last call of start().
    double stop()
    {
        if (running == false)
        {
            std::cout << "WARNING: ntime::stop() has been called without calling "
                 << "ntime::start() first." << std::endl;
            return 0.0;
        }
        else
        {
            std::chrono::time_point<std::chrono::high_resolution_clock> t_end;
            t_end = std::chrono::high_resolution_clock::now();
            auto duration = t_end - t_start;
            auto sec = std::chrono::duration_cast<std::chrono::seconds>(duration);
            duration -= sec;
            auto nsec = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
            double t = (double)sec.count() + nsec.count() / 1e9;
            counter++;
            elapsed += t;
            elapsed_sq += t * t;
            return t;
        }
    }

    double total()
    {
        return elapsed;
    }

    double average()
    {
        return elapsed / double(counter);
    }

    double deviation()
    {
        if (counter == 1)
        {
            return 0.0;
        }
        else
        {
            return sqrt(elapsed_sq / double(counter) - average() * average());
        }
    }

    unsigned count()
    {
        return counter;
    }
};

namespace get_time
{
    extern std::map<std::string, ntime> timers;

    bool is_master();
    
    void reset();

    void print();

    std::string sorted_output();

    void start(std::string name);

    void stop(std::string name);

    double total(std::string name);

    double average(std::string name);

    double deviation(std::string name);
}

#endif