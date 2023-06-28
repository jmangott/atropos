#ifndef TEST_PARAMETERS_HPP
#define TEST_PARAMETERS_HPP

#include <cmath>
#include <iostream>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"
#include "kl_step_functions.hpp"
#include "partition_class.hpp"
#include "reaction_class.hpp"
#include "weight_functions.hpp"

vector<string> test_names = {"S1, S2"};
mysys test_system(test_names);
myreact test_react0(
    {-1, 0}, {0}, [](std::vector<double> y)
    { return y[0]; },
    test_system);
myreact test_react1(
    {0, -1}, {1}, [](std::vector<double> y)
    { return y[1]; },
    test_system);
myreact test_react2(
    {1, 0}, {1}, [](std::vector<double> y)
    { return 1.0 / (1.0 + y[1]); },
    test_system);
myreact test_react3(
    {0, 1}, {0}, [](std::vector<double> y)
    { return 1.0 / (1.0 + y[0]); },
    test_system);
vector<Index> sigma1(test_system.mu()), sigma2(test_system.mu());

const Index kR = 2;
const Index kN = 2;
const Index kBinsize = 1;
std::vector<double> kLiml1 {0.0};
std::vector<double> kLiml2 {0.0};

const Index kM1 = 1;
const Index kM2 = 1;
Index nsteps = 1;
double tstar = 1.0;
double tau = tstar / nsteps;

double min_prop, max_prop;

grid_info grid(kM1, kM2, kR, kN, kBinsize, kLiml1, kLiml2);

// Low rank structure (for storing X1, X2 and S)
lr2<double> lr_sol(kR, {grid.dx1, grid.dx2});

// Inner products
std::function<double(double *, double *)> ip_xx1 = inner_product_from_const_weight((double) grid.h1_mult, grid.dx1);
std::function<double(double *, double *)> ip_xx2 = inner_product_from_const_weight((double) grid.h2_mult, grid.dx2);
blas_ops blas;

std::vector<multi_array<double, 2>> w_x_dep(test_system.mu());

// Coefficients
vector<multi_array<double, 3>> c_coeff1(test_system.mu()), d_coeff1(test_system.mu());
vector<multi_array<double, 3>> c_coeff2(test_system.mu()), d_coeff2(test_system.mu());
multi_array<double, 4> e_coeff({grid.r, grid.r, grid.r, grid.r}), f_coeff({grid.r, grid.r, grid.r, grid.r});

// Perform partition of the network
partition_info<1> partition1(grid, test_system);
partition_info<2> partition2(grid, test_system);

inline void InitializeTest(lr2<double> &lr_sol, grid_info grid, std::function<double(double *, double *)> ip_xx1, std::function<double(double *, double *)> ip_xx2, blas_ops blas, std::vector<multi_array<double, 2>> &w_x_dep)
{
    multi_array<double, 1> xx1({grid.dx1});
    multi_array<double, 1> xx2({grid.dx2});

    std::vector<const double *> x1, x2;
    x1.reserve(1);
    x2.reserve(1);

    CalculateShiftAmount(sigma1, sigma2, test_system, grid);

    // Allocate memory
    for (Index mu = 0; mu < test_system.mu(); mu++)
    {
        w_x_dep[mu].resize({partition1.dx_dep(mu), partition2.dx_dep(mu)});
        c_coeff1[mu].resize({partition1.dx_dep(mu), grid.r, grid.r});
        d_coeff1[mu].resize({partition1.dx_dep(mu), grid.r, grid.r});
        c_coeff2[mu].resize({partition2.dx_dep(mu), grid.r, grid.r});
        d_coeff2[mu].resize({partition2.dx_dep(mu), grid.r, grid.r});
    }

    xx1(0) = std::exp(-std::pow(-0.5, 2));
    xx1(1) = std::exp(-std::pow( 0.5, 2));
    xx2(0) = std::exp(-std::pow(-0.5, 2));
    xx2(1) = std::exp(-std::pow( 0.5, 2));

    // Point to the beginning of every column of X1 and X2
    double *it1 = xx1.begin();
    double *it2 = xx2.begin();
    x1.push_back(it1);
    x2.push_back(it2);

    initialize(lr_sol, x1, x2, ip_xx1, ip_xx2, blas);
    double norm = 1.0 / std::sqrt(2);
    double two_sqrt_e = 2.0 / std::exp(0.5);

    // Reset X1, S and X2 to guarantee that they always stay the same
    lr_sol.X(0, 0) = norm;
    lr_sol.X(0, 1) = norm;
    lr_sol.X(1, 0) = norm;
    lr_sol.X(1, 1) = -norm;

    lr_sol.V(0, 0) = norm;
    lr_sol.V(0, 1) = norm;
    lr_sol.V(1, 0) = norm;
    lr_sol.V(1, 1) = -norm;

    lr_sol.S(0, 0) = two_sqrt_e;
    lr_sol.S(1, 1) = 0.0;

    // Calculate the integration weights
    CalculateWeightDep(w_x_dep, min_prop, max_prop, test_system, grid, partition1, partition2);
}

#endif