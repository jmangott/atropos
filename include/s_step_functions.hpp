#ifndef S_STEP_FUNCTIONS_HPP
#define S_STEP_FUNCTIONS_HPP

#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"
#include "kl_step_functions.hpp"
#include "reaction_class.hpp"
#include "timer_class.hpp"

// Calculate the coefficients E and F (depending on `mu`)
void CalculateCoefficientsS(multi_array<double, 4> &e_coeff_tot, multi_array<double, 4> &f_coeff_tot, const std::vector<multi_array<double, 3>> &c_coeff_dep, const std::vector<multi_array<double, 3>> &d_coeff_dep, std::vector<Index> sigma1, std::vector<Index> sigma2, const lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, partition_info<1> partition1, partition_info<2> partition2, const vector<multi_array<double, 2>> &w_x_dep);

// Perform S-Step with time step size `tau`
void PerformSStep(multi_array<double, 2> &s_dot, const multi_array<double, 2> &s, const multi_array<double, 4> &e_coeff, const multi_array<double, 4> &f_coeff, std::vector<Index> sigma1, std::vector<Index> sigma2, blas_ops blas, mysys reaction_system, grid_info grid, partition_info<1> partition1, partition_info<2> partition2, const vector<multi_array<double, 2>> &w_x_dep, double tau);

#endif