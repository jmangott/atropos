#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <cmath>
#include <iostream>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"
#include "k_step_functions.hpp"
#include "reaction_class.hpp"

using std::cout;
using std::endl;
using std::vector;

TEST_CASE("k_step", "[k_step]")
{
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

    Index r = 2;
    Index n = 2;
    Index k = 1;
    Index m1 = 1;
    Index m2 = 1;
    Index nsteps = 1;
    double tstar = 1.0;
    double tau = tstar / nsteps;

    multi_array<Index, 1> n_xx1({1}), n_xx2({1});
    multi_array<Index, 1> k_xx1({1}), k_xx2({1});
    multi_array<double, 1> lim_xx1({1}), lim_xx2({1});
    multi_array<double, 1> h_xx1({1}), h_xx2({1});
    Index dxx1_mult, dxx2_mult;

    InitializeAuxiliaryObjects(n_xx1, n_xx2, k_xx1, k_xx2, lim_xx1, lim_xx2, h_xx1, h_xx2, dxx1_mult, dxx2_mult, m1, m2, n, k);

    // SECTION("InitializeAuxiliaryObjects")
    // {
        multi_array<double, 1> lim_xx1_comparison({1}), lim_xx2_comparison({1});
        multi_array<double, 1> h_xx1_comparison({1}), h_xx2_comparison({1});
        Index dxx1_mult_comparison = 2;
        Index dxx2_mult_comparison = 2;
        lim_xx1_comparison(0) = 1.0;
        lim_xx2_comparison(0) = 1.0;
        h_xx1_comparison(0) = 1.0;
        h_xx2_comparison(0) = 1.0;

        REQUIRE(bool(lim_xx1 == lim_xx1_comparison));
        REQUIRE(bool(lim_xx2 == lim_xx2_comparison));
        REQUIRE(bool(h_xx1 == h_xx1_comparison));
        REQUIRE(bool(h_xx2 == h_xx2_comparison));
        REQUIRE(bool(dxx1_mult == dxx1_mult_comparison));
        REQUIRE(bool(dxx2_mult == dxx2_mult_comparison));
    // }

    // Coefficients
    multi_array<double, 2> c2({r, r});
    multi_array<double, 2> d2({r, r});

    // Integration weight
    multi_array<double, 1> w_x2({dxx2_mult});

    // Temporary objects for multiplication
    multi_array<double, 2> tmp_x({dxx1_mult, r});

    // Objects for setting up X1 and X2 for t = 0
    multi_array<double, 1> xx1({dxx1_mult});
    multi_array<double, 1> xx2({dxx1_mult});
    vector<const double *> x1, x2;

    // Low rank structure (for storing X1, X2 and S)
    lr2<double> lr_sol(r, {dxx1_mult, dxx2_mult});

    // Inner products
    std::function<double(double *, double *)> ip_xx1;
    std::function<double(double *, double *)> ip_xx2;
    blas_ops blas;

    multi_array<double, 1> state_vec1({1}), state_vec2({1});
    multi_array<Index, 1> vec_index1({1}), vec_index2({1});

    for (Index i = 0; i < n_xx1(0); i++)
    {
        vec_index1 = CombIndexToVecIndex(i, n_xx1);
        state_vec1 = VecIndexToState(vec_index1, n_xx1, lim_xx1);
        xx1(i) = std::exp(-std::pow(state_vec1(0) - 0.5, 2));
    }
    for (Index i = 0; i < n_xx2(0); i++)
    {
        vec_index2 = CombIndexToVecIndex(i, n_xx2);
        state_vec2 = VecIndexToState(vec_index2, n_xx2, lim_xx2);
        xx2(i) = std::exp(-std::pow(state_vec2(0) - 0.5, 2));
    }

    // Point to the beginning of every column of X1 and X2
    double *it1 = xx1.begin();
    double *it2 = xx2.begin();
    x1.push_back(it1);
    x2.push_back(it2);

    // Set up the low-rank structure and the inner products
    ip_xx1 = inner_product_from_const_weight(h_xx1(0), dxx1_mult);
    ip_xx2 = inner_product_from_const_weight(h_xx2(0), dxx2_mult);
    initialize(lr_sol, x1, x2, ip_xx1, ip_xx2, blas);

    multi_array<double, 2> x1x2_comparison({2, 2}), k_comparison({2, 2});
    double norm = 1.0 / std::sqrt(2);
    double norm_2e = std::sqrt(2.0 / std::exp(1.0));

    x1x2_comparison(0, 0) = norm;
    x1x2_comparison(0, 1) = norm;
    x1x2_comparison(1, 0) = norm;
    x1x2_comparison(1, 1) = -norm;

    k_comparison(0, 0) = norm_2e;
    k_comparison(0, 1) = 0.0;
    k_comparison(1, 0) = norm_2e;
    k_comparison(1, 1) = 0.0;

    REQUIRE(bool(lr_sol.X == x1x2_comparison));
    REQUIRE(bool(lr_sol.V == x1x2_comparison));

    tmp_x = lr_sol.X;
    blas.matmul(tmp_x, lr_sol.S, lr_sol.X); // lr_sol.X contains now K

    REQUIRE(bool(lr_sol.X == k_comparison));

    CalculateShiftAmount(sigma1, sigma2, test_system, n_xx1, n_xx2, k_xx1, k_xx2);

    // SECTION("CalculateShiftAmount")
    // {
        vector<Index> sigma1_comparison, sigma2_comparison;
        sigma1_comparison = {-1, 0, 1, 0};
        sigma2_comparison = {0, -1, 0, 1};

        REQUIRE(bool(sigma1 == sigma1_comparison));
        REQUIRE(bool(sigma2 == sigma2_comparison));
    // }

    // SECTION("CalculateWeightX2")
    // {
        multi_array<double, 1> weight0({2}), weight1({2}), weight2({2}), weight3({2});
        multi_array<double, 1> weight0_comparison({2}), weight1_comparison({2}), weight2_comparison({2}), weight3_comparison({2});

        // CASE 0: x_1 = 0.0, CASE 1: x_1 = 1.0
        vec_index1(0) = 0;
        state_vec1 = VecIndexToState(vec_index1, n_xx1, lim_xx1);

        weight0 = CalculateWeightX2(n_xx1, n_xx2, lim_xx1, lim_xx2, dxx2_mult, h_xx2, vec_index1, test_system, 0);
        weight1 = CalculateWeightX2(n_xx1, n_xx2, lim_xx1, lim_xx2, dxx2_mult, h_xx2, vec_index1, test_system, 1);
        weight2 = CalculateWeightX2(n_xx1, n_xx2, lim_xx1, lim_xx2, dxx2_mult, h_xx2, vec_index1, test_system, 2);
        weight3 = CalculateWeightX2(n_xx1, n_xx2, lim_xx1, lim_xx2, dxx2_mult, h_xx2, vec_index1, test_system, 3);

        weight0_comparison(0) = state_vec1(0);
        weight0_comparison(1) = state_vec1(0);
        weight1_comparison(0) = 0.0;
        weight1_comparison(1) = 1.0;
        weight2_comparison(0) = 1.0;
        weight2_comparison(1) = 0.5;
        weight3_comparison(0) = 1.0 / (1.0 + state_vec1(0));
        weight3_comparison(1) = 1.0 / (1.0 + state_vec1(0));

        REQUIRE(bool(weight0 == weight0_comparison));
        REQUIRE(bool(weight1 == weight1_comparison));
        REQUIRE(bool(weight2 == weight2_comparison));
        REQUIRE(bool(weight3 == weight3_comparison));
    // }

    // SECTION("CalculateCoefficientsX2")
    // {
        // CASE 0: x_1 = 0.0, CASE 1: x_1 = 1.0
        vec_index1(0) = 0;
        state_vec1 = VecIndexToState(vec_index1, n_xx1, lim_xx1);

        multi_array<double, 2> c2_0({2, 2}), c2_1({2, 2}), c2_2({2, 2}), c2_3({2, 2});
        multi_array<double, 2> d2_0({2, 2}), d2_1({2, 2}), d2_2({2, 2}), d2_3({2, 2});

        multi_array<double, 2> c2_0_comparison({2, 2}), c2_1_comparison({2, 2}), c2_2_comparison({2, 2}), c2_3_comparison({2, 2});
        multi_array<double, 2> d2_0_comparison({2, 2}), d2_1_comparison({2, 2}), d2_2_comparison({2, 2}), d2_3_comparison({2, 2});

        CalculateCoefficientsX2(c2_0, d2_0, n_xx1, n_xx2, lim_xx1, lim_xx2, h_xx2, lr_sol, blas, sigma2[0], vec_index1, test_system, 0);
        CalculateCoefficientsX2(c2_1, d2_1, n_xx1, n_xx2, lim_xx1, lim_xx2, h_xx2, lr_sol, blas, sigma2[1], vec_index1, test_system, 1);
        CalculateCoefficientsX2(c2_2, d2_2, n_xx1, n_xx2, lim_xx1, lim_xx2, h_xx2, lr_sol, blas, sigma2[2], vec_index1, test_system, 2);
        CalculateCoefficientsX2(c2_3, d2_3, n_xx1, n_xx2, lim_xx1, lim_xx2, h_xx2, lr_sol, blas, sigma2[3], vec_index1, test_system, 3);

        c2_0_comparison(0, 0) = state_vec1(0);
        c2_0_comparison(0, 1) = 0.0;
        c2_0_comparison(1, 0) = state_vec1(0);
        c2_0_comparison(1, 1) = 0.0;

        c2_1_comparison(0, 0) = 0.5;
        c2_1_comparison(0, 1) =-0.5;
        c2_1_comparison(1, 0) =-0.5;
        c2_1_comparison(1, 1) = 0.5;

        c2_2_comparison(0, 0) = 0.75;
        c2_2_comparison(0, 1) = 0.25;
        c2_2_comparison(1, 0) = 0.25;
        c2_2_comparison(1, 1) = 0.75;

        c2_3_comparison(0, 0) = 1.0 / (1.0 + state_vec1(0));
        c2_3_comparison(0, 1) = 0.0;
        c2_3_comparison(1, 0) = 1.0 / (1.0 + state_vec1(0));
        c2_3_comparison(1, 1) = 0.0;

        d2_0_comparison = c2_0_comparison;
        d2_1_comparison = c2_1_comparison;
        d2_2_comparison = c2_2_comparison;

        d2_3_comparison(0, 0) = 1.0 / (1.0 + state_vec1(0));
        d2_3_comparison(0, 1) = 0.0;
        d2_3_comparison(1, 0) = 0.0;
        d2_3_comparison(1, 1) = 1.0 / (1.0 + state_vec1(0));

        REQUIRE(bool(c2_0 == c2_0_comparison));
        REQUIRE(bool(c2_1 == c2_1_comparison));
        REQUIRE(bool(c2_2 == c2_2_comparison));
        REQUIRE(bool(c2_3 == c2_3_comparison));

        REQUIRE(bool(d2_0 == d2_0_comparison));
        REQUIRE(bool(d2_1 == d2_1_comparison));
        REQUIRE(bool(d2_2 == d2_2_comparison));
        REQUIRE(bool(d2_3 == d2_3_comparison));
    // }

    // SECTION("PerformKStep")
    // {
        multi_array<double, 2> ktau_comparison({2, 2});
        ktau_comparison = lr_sol.X;
        ktau_comparison(1, 0) += -norm_2e;
        cout << "K_comp" << endl;
        cout << ktau_comparison(0, 0) << " " << ktau_comparison(0, 1) << endl;
        cout << ktau_comparison(1, 0) << " " << ktau_comparison(1, 1) << endl;
        cout << endl;

        PerformKStep(n_xx1, n_xx2, h_xx2, lim_xx1, lim_xx2, sigma1, sigma2, lr_sol, blas, test_system, tau);
        cout << "K(tau):" << endl;
        cout << lr_sol.X(0, 0) << " " << lr_sol.X(0, 1) << endl;
        cout << lr_sol.X(1, 0) << " " << lr_sol.X(1, 1) << endl;
        cout << endl;
        REQUIRE(bool(lr_sol.X == ktau_comparison));
    // }
}