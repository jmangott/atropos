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

// Print progress bar
void PrintProgressBar(const Index ts, const Index kNsteps,
                      const std::chrono::system_clock::time_point t_start,
                      const double norm);

struct diagnostics {
    diagnostics(const integrator_base& _integrator,
                const std::chrono::nanoseconds _t_elapsed, const double _tau,
                const double _dm_max)
        : integrator(_integrator), t_elapsed(_t_elapsed), tau(_tau), dm_max(_dm_max)
    {
    }

    const integrator_base integrator;
    const std::chrono::nanoseconds t_elapsed;
    const double tau;
    const double dm_max;
};

// Print diagnostic information
std::ostream& operator<<(std::ostream& os, const diagnostics& dgn);

#endif
