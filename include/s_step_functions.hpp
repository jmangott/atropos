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


// Calculate the auxiliary coefficients B (`b_coeff_vec`) and B_tilde (`b_coeff_vec_shift`) 
// for coefficients E and F (depending on `mu`)
void CalculateCoefficientsB(multi_array<double, 3> &b_coeff_vec_shift, multi_array<double, 3> &b_coeff_vec, const lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, partition_info1<1> partition1, partition_info1<2> partition2, Index mu, std::vector<Index> sigma2, const vector<multi_array<double, 2>> &w_x);

// Calculate the coefficients E and F (depending on `mu`)
void CalculateCoefficientsS(multi_array<double, 4> &e_coeff_tot, multi_array<double, 4> &f_coeff_tot, const multi_array<double, 3> &b_coeff_vec_shift, const multi_array<double, 3> &b_coeff_vec, const lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, partition_info1<1> partition1, Index mu, std::vector<Index> sigma1);

// Perform S-Step with time step size `tau`
void PerformSStep(std::vector<Index> sigma1, std::vector<Index> sigma2, lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, partition_info1<1> partition1, partition_info1<2> partition2, const vector<multi_array<double, 2>> &w_x, double tau);

#endif