#ifndef PRINT_FUNCTIONS_HPP
#define PRINT_FUNCTIONS_HPP

#include <chrono>
#include <iostream>
#include <sstream>
#include <string>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "integration_methods.hpp"

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
void PrintProgressBar(const Index ts, const Index kNsteps, const std::chrono::system_clock::time_point t_start, const double norm);

// Print diagnostic information
void PrintDiagnostics(const integration_method &method, const std::chrono::nanoseconds t_elapsed, const double tau, const double dm_max);

#endif