#ifndef KL_STEP_FUNCTIONS_HPP
#define KL_STEP_FUNCTIONS_HPP

#include <array>
#include <memory>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"
#include "reaction_class.hpp"

// Calculate the integration weight for coefficients C1 and D1 (`id` = 1) or C2 and D2 (`id` = 2), depending on `alpha1` and `mu`
multi_array<double, 1> CalculateWeightX(int id, multi_array<Index, 1> vec_index_c, mysys reaction_system, grid_info grid, Index mu);

// Calculate coefficients C1 and D1 (`id` = 1) or C2 and D2 (`id` = 2) for all values of `dep_vec` for a given reaction `mu`
void CalculateCoefficientsX(int id, multi_array<double, 2> &c_coeff, multi_array<double, 2> &d_coeff, const lr2<double> &lr_sol, blas_ops blas, const multi_array<double, 2> &xx2_shift, multi_array<Index, 1> vec_index1, mysys reaction_system, grid_info grid, Index mu);

// Perform L-Step (`id` = 1) or K-Step (`id` = 2) with time step size `tau`
void PerformKLStep(int id, std::vector<Index> sigma1, vector<Index> sigma2, lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, double tau);

#endif