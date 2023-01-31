#ifndef WEIGHT_FUNCTIONS_HPP
#define WEIGHT_FUNCTIONS_HPP

#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"
#include "kl_step_functions.hpp"
#include "reaction_class.hpp"

// Calculate the integration weight for coefficients B1, B2, C1, C2, D1 and D2
// NOTE: the result has to be multiplied with h1_mult or h2_mult!
void CalculateWeight(multi_array<double, 3> &w_x, mysys reaction_system, grid_info grid);

void CalculateWeightDep(vector<multi_array<double, 2>> &w_x_dep, mysys reaction_system, grid_info grid, partition_info partition);

#endif