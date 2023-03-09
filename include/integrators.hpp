#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP

#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"
#include "io_functions.hpp"
#include "kl_step_functions.hpp"
#include "partition_class.hpp"
#include "print_functions.hpp"
#include "s_step_functions.hpp"
#include "timer_class.hpp"
#include "weight_functions.hpp"

double CalculateNorm(lr2<double> &lr_sol, grid_info &grid);


// Perform integration for one timestep `tau` with Lie-Trotter splitting and explicit Euler
void IntegrateFirstOrder(lr2<double> &lr_sol, const std::vector<multi_array<double, 2>> &w_x_dep, std::vector<multi_array<double, 3>> &c_coeff1, std::vector<multi_array<double, 3>> &d_coeff1, std::vector<multi_array<double, 3>> &c_coeff2, std::vector<multi_array<double, 3>> &d_coeff2, multi_array<double, 5> &e_coeff, multi_array<double, 5> &f_coeff, const std::vector<Index> sigma1, const std::vector<Index> sigma2, mysys &mysystem, grid_info &grid, partition_info<1> &partition1, partition_info<2> &partition2, std::function<double(double *, double *)> ip_xx1, std::function<double(double *, double *)> ip_xx2, blas_ops &blas, double tau);


// Perform integration for one timestep `tau` with Strang splitting and explicit Euler with substeps
void IntegrateSecondOrder(lr2<double> &lr_sol, const vector<multi_array<double, 2>> &w_x_dep, std::vector<multi_array<double, 3>> &c_coeff1, std::vector<multi_array<double, 3>> &d_coeff1, std::vector<multi_array<double, 3>> &c_coeff2, std::vector<multi_array<double, 3>> &d_coeff2, multi_array<double, 5> &e_coeff, multi_array<double, 5> &f_coeff, const vector<Index> sigma1, const vector<Index> sigma2, mysys &mysystem, grid_info &grid, partition_info<1> &partition1, partition_info<2> &partition2, std::function<double(double *, double *)> ip_xx1, std::function<double(double *, double *)> ip_xx2, blas_ops &blas, double tau, Index n_substeps);

#endif