#ifndef K_STEP_FUNCTIONS_HPP
#define K_STEP_FUNCTIONS_HPP

#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"
#include "k_step_functions.hpp"
#include "reaction_class.hpp"


// Calculate the integration weight for coefficients C2 and D2 (depending on `alpha1` and `mu`)
multi_array<double, 1> CalculateWeightX2(multi_array<Index, 1> n_xx1, multi_array<Index, 1> n_xx2, multi_array<double, 1> lim_xx1, multi_array<double, 1> lim_xx2, Index dxx2_mult, multi_array<double, 1> h_xx2, multi_array<Index, 1> vec_index1, mysys reaction_system, Index mu);


// Calculate coefficients C2 and D2 for all values of `dep_vec` for a given reaction `mu`
void CalculateCoefficientsX2(multi_array<double, 2> &c2, multi_array<double, 2> &d2, multi_array<Index, 1> n_xx1, multi_array<Index, 1> n_xx2, multi_array<double, 1> lim_xx1, multi_array<double, 1> lim_xx2, multi_array<double, 1> h_xx2, lr2<double> &lr_sol, blas_ops blas, Index shift, multi_array<Index, 1> vec_index1, mysys reaction_system, Index mu);


// Perform K-Step with time step size `tau`
void PerformKStep(multi_array<Index, 1> n_xx1, multi_array<Index, 1> n_xx2, multi_array<double, 1> h_xx2, multi_array<double, 1> lim_xx1, multi_array<double, 1> lim_xx2, std::vector<Index> sigma1, std::vector<Index> sigma2, lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, double tau);

#endif