#ifndef PRINT_FUNCTIONS_HPP
#define PRINT_FUNCTIONS_HPP

#include <chrono>
#include <iostream>
#include <sstream>
#include <string>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

// Print progress bar
void PrintProgressBar(Index ts, Index kNsteps, std::chrono::_V2::system_clock::time_point start_time);

#endif