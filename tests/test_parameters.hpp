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
#include "reaction_class.hpp"

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
vector<Index> sigma1, sigma2;

const Index kR = 2;
const Index kN = 2;
const Index kK = 1;
const Index kM1 = 1;
const Index kM2 = 1;
Index nsteps = 1;
double tstar = 1.0;
double tau = tstar / nsteps;

grid_info grid(kM1, kM2, kR, kN, kK);

// Low rank structure (for storing X1, X2 and S)
lr2<double> lr_sol(kR, {grid.dx1, grid.dx2});

// Inner products
std::function<double(double *, double *)> ip_xx1 = inner_product_from_const_weight(grid.h1(0), grid.dx1);
std::function<double(double *, double *)> ip_xx2 = inner_product_from_const_weight(grid.h2(0), grid.dx2);
blas_ops blas;


inline void InitializeTest(lr2<double> &lr_sol, grid_info grid, std::function<double(double *, double *)> ip_xx1, std::function<double(double *, double *)> ip_xx2, blas_ops blas)
{
    multi_array<double, 1> xx1({grid.dx1});
    multi_array<double, 1> xx2({grid.dx2});

    std::vector<const double *> x1, x2;
    x1.reserve(1);
    x2.reserve(1);

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
}

#endif