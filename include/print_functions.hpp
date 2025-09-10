#ifndef PRINT_FUNCTIONS_HPP
#define PRINT_FUNCTIONS_HPP

#include <chrono>
#include <iostream>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "integrators.hpp"

#ifdef __OPENMP__
#include <omp.h>
#endif

template <class Rep, std::intmax_t num, std::intmax_t denom>
auto ChronoBurst(std::chrono::duration<Rep, std::ratio<num, denom>> d)
{
    const auto hrs = duration_cast<std::chrono::hours>(d);
    const auto mins = duration_cast<std::chrono::minutes>(d - hrs);
    const auto secs = duration_cast<std::chrono::seconds>(d - hrs - mins);
    const auto ms = duration_cast<std::chrono::milliseconds>(d - hrs - mins - secs);

    return std::make_tuple(hrs, mins, secs, ms);
}

// TODO: output width is not always the same
// TODO: calculate `time_left` with a moving mean (of the last n steps)
template <typename T>
void PrintProgressBar(const Index ts, const Index kNsteps, const T t_start, const double norm)
{
    int bar_width = 30;
    double progress = (ts + 1.0) / kNsteps;
    int pos = bar_width * progress;
    string time_unit;
    string progress_bar(bar_width, ' ');

    auto t_stop(std::chrono::high_resolution_clock::now());
    auto duration(std::chrono::duration_cast<std::chrono::seconds>(t_stop - t_start));
    // time_per_step = duration.count() / (ts + 1.0);
    auto time_per_step = duration / (ts + 1.0);
    auto time_left = time_per_step * (kNsteps - 1.0 - ts);
    double time_per_step_count;

    if (time_per_step.count() < 0.01)
    {
        time_unit = "ms";
        time_per_step_count = time_per_step.count() * 1000.0;
    }
    else if (time_per_step.count() < 60.0)
    {
        time_unit = "s";
        time_per_step_count = time_per_step.count();
    }
    else
    {
        time_unit = "min";
        time_per_step_count = time_per_step.count() / 60.0;
    }

    const auto [hrs, mins, secs, ms] = ChronoBurst(time_left);

    std::fill(std::begin(progress_bar), std::begin(progress_bar) + pos, '#');
    std::fill(std::begin(progress_bar) + pos, std::end(progress_bar), '-');

    printf("[%*s], step: %ti/%ti, time per step: %.2f%*s, time left: %2.2lli:%2.2lli:%2.2lli, progress: %4.2f%%, |norm(P)-1|: %3.2e\r", bar_width, progress_bar.c_str(), ts + 1, kNsteps, time_per_step_count, (int)time_unit.size(), time_unit.c_str(), (long long int)hrs.count(), (long long int)mins.count(), (long long int)secs.count(), progress * 100, std::abs(norm - 1.0));
    fflush(stdout);
}

struct diagnostics
{
    diagnostics(const integrator_base &_integrator, const std::chrono::nanoseconds _t_elapsed, const double _tau, const double _dm_max)
    : integrator(_integrator)
    , t_elapsed(_t_elapsed)
    , tau(_tau)
    , dm_max(_dm_max)
    {}

    const integrator_base integrator;
    const std::chrono::nanoseconds t_elapsed;
    const double tau;
    const double dm_max;
};

// Print diagnostic information
std::ostream& operator<<(std::ostream &os, const diagnostics &dgn);

#endif