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
#include "reaction_class.hpp"
#include "test_parameters.hpp"

using std::cout;
using std::endl;
using std::vector;

TEST_CASE("k_step", "[k_step]")
{
    // Temporary objects for multiplication
    multi_array<double, 2> tmp_x1({grid.dx1, grid.r});

    InitializeTest(lr_sol, grid, ip_xx1, ip_xx2, blas, w_x_dep);

    // SECTION("InitializeAuxiliaryObjects")
    // {
        multi_array<double, 1> lim1_comparison({1}), lim2_comparison({1});
        multi_array<double, 1> h1_comparison({1}), h2_comparison({1});
        Index dx1_comparison = 2;
        Index dx2_comparison = 2;
        lim1_comparison(0) = 1.0;
        lim2_comparison(0) = 1.0;
        h1_comparison(0) = 1.0;
        h2_comparison(0) = 1.0;

        REQUIRE(bool(grid.limr1 == lim1_comparison));
        REQUIRE(bool(grid.limr2 == lim2_comparison));
        REQUIRE(bool(grid.h1 == h1_comparison));
        REQUIRE(bool(grid.h2 == h2_comparison));
        REQUIRE(bool(grid.dx1 == dx1_comparison));
        REQUIRE(bool(grid.dx2 == dx2_comparison));
    // }

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

    tmp_x1 = lr_sol.X;
    blas.matmul(tmp_x1, lr_sol.S, lr_sol.X); // lr_sol.X contains now K

    REQUIRE(bool(lr_sol.X == k_comparison));

    CalculateShiftAmount(sigma1, sigma2, test_system, grid);

    multi_array<double, 1> state_vec1({1});
    multi_array<Index, 1> vec_index1({1});

    // SECTION("CalculateShiftAmount")
    // {
        vector<Index> sigma1_comparison(test_system.mu()), sigma2_comparison(test_system.mu());
        sigma1_comparison = {-1, 0, 1, 0};
        sigma2_comparison = {0, -1, 0, 1};

        REQUIRE(bool(sigma1 == sigma1_comparison));
        REQUIRE(bool(sigma2 == sigma2_comparison));
    // }

    // SECTION("CalculateWeightX")
    // {
        multi_array<double, 2> weight0_comparison({2, 1}), weight1_comparison({1, 2}), weight2_comparison({1, 2}), weight3_comparison({2, 1});

        weight0_comparison(0, 0) = 0.0;
        weight0_comparison(1, 0) = 1.0;

        weight1_comparison(0, 0) = 0.0;
        weight1_comparison(0, 1) = 1.0;

        weight2_comparison(0, 0) = 1.0;
        weight2_comparison(0, 1) = 0.5;

        weight3_comparison(0, 0) = 1.0;
        weight3_comparison(1, 0) = 0.5;

        REQUIRE(bool(w_x_dep[0] == weight0_comparison));
        REQUIRE(bool(w_x_dep[1] == weight1_comparison));
        REQUIRE(bool(w_x_dep[2] == weight2_comparison));
        REQUIRE(bool(w_x_dep[3] == weight3_comparison));
    // }

    // SECTION("CalculateCoefficientsX")
    // {
        multi_array<double, 2> c2_0_0({2, 2}), c2_0_1({2, 2}), c2_1({2, 2}), c2_2({2, 2}), c2_3_0({2, 2}), c2_3_1({2, 2});
        multi_array<double, 2> d2_0_0({2, 2}), d2_0_1({2, 2}), d2_1({2, 2}), d2_2({2, 2}), d2_3_0({2, 2}), d2_3_1({2, 2});

        multi_array<double, 2> c2_0_0_comparison({2, 2}), c2_0_1_comparison({2, 2}), c2_1_comparison({2, 2}), c2_2_comparison({2, 2}), c2_3_0_comparison({2, 2}), c2_3_1_comparison({2, 2});
        multi_array<double, 2> d2_0_0_comparison({2, 2}), d2_0_1_comparison({2, 2}), d2_1_comparison({2, 2}), d2_2_comparison({2, 2}), d2_3_0_comparison({2, 2}), d2_3_1_comparison({2, 2});

        multi_array<double, 2> xx2_shift0(lr_sol.V.shape());
        multi_array<double, 2> xx2_shift1(lr_sol.V.shape());
        multi_array<double, 2> xx2_shift2(lr_sol.V.shape());
        multi_array<double, 2> xx2_shift3(lr_sol.V.shape());

        ShiftMultiArrayRows(2, xx2_shift0, lr_sol.V, -sigma2[0], test_system.reactions[0]->minus_nu, grid, test_system);
        ShiftMultiArrayRows(2, xx2_shift1, lr_sol.V, -sigma2[1], test_system.reactions[1]->minus_nu, grid, test_system);
        ShiftMultiArrayRows(2, xx2_shift2, lr_sol.V, -sigma2[2], test_system.reactions[2]->minus_nu, grid, test_system);
        ShiftMultiArrayRows(2, xx2_shift3, lr_sol.V, -sigma2[3], test_system.reactions[3]->minus_nu, grid, test_system);

        CalculateCoefficientsX<1>(c2_0_0, d2_0_0, lr_sol, blas, xx2_shift0, 0, test_system, grid, partition2, 0, w_x_dep);
        CalculateCoefficientsX<1>(c2_0_1, d2_0_1, lr_sol, blas, xx2_shift0, 1, test_system, grid, partition2, 0, w_x_dep);
        CalculateCoefficientsX<1>(c2_1, d2_1, lr_sol, blas, xx2_shift1, 0, test_system, grid, partition2, 1, w_x_dep);
        CalculateCoefficientsX<1>(c2_2, d2_2, lr_sol, blas, xx2_shift2, 0, test_system, grid, partition2, 2, w_x_dep);
        CalculateCoefficientsX<1>(c2_3_0, d2_3_0, lr_sol, blas, xx2_shift3, 0, test_system, grid, partition2, 3, w_x_dep);
        CalculateCoefficientsX<1>(c2_3_1, d2_3_1, lr_sol, blas, xx2_shift3, 1, test_system, grid, partition2, 3, w_x_dep);

        c2_0_0_comparison(0, 0) = 0.0;
        c2_0_0_comparison(0, 1) = 0.0;
        c2_0_0_comparison(1, 0) = 0.0;
        c2_0_0_comparison(1, 1) = 0.0;

        c2_0_1_comparison(0, 0) = 1.0;
        c2_0_1_comparison(0, 1) = 0.0;
        c2_0_1_comparison(1, 0) = 0.0;
        c2_0_1_comparison(1, 1) = 1.0;

        c2_1_comparison(0, 0) = 0.5;
        c2_1_comparison(0, 1) =-0.5;
        c2_1_comparison(1, 0) = 0.5;
        c2_1_comparison(1, 1) =-0.5;

        c2_2_comparison(0, 0) = 0.75;
        c2_2_comparison(0, 1) = 0.25;
        c2_2_comparison(1, 0) = 0.25;
        c2_2_comparison(1, 1) = 0.75;

        c2_3_0_comparison(0, 0) = 0.5;
        c2_3_0_comparison(0, 1) = 0.5;
        c2_3_0_comparison(1, 0) =-0.5;
        c2_3_0_comparison(1, 1) =-0.5;

        c2_3_1_comparison(0, 0) = 0.25;
        c2_3_1_comparison(0, 1) = 0.25;
        c2_3_1_comparison(1, 0) =-0.25;
        c2_3_1_comparison(1, 1) =-0.25;

        d2_0_0_comparison = c2_0_0_comparison;
        d2_0_1_comparison = c2_0_1_comparison;

        d2_1_comparison(0, 0) = 0.5;
        d2_1_comparison(0, 1) =-0.5;
        d2_1_comparison(1, 0) =-0.5;
        d2_1_comparison(1, 1) = 0.5;

        d2_2_comparison = c2_2_comparison;

        d2_3_0_comparison(0, 0) = 1.0;
        d2_3_0_comparison(0, 1) = 0.0;
        d2_3_0_comparison(1, 0) = 0.0;
        d2_3_0_comparison(1, 1) = 1.0;

        d2_3_1_comparison(0, 0) = 0.5;
        d2_3_1_comparison(0, 1) = 0.0;
        d2_3_1_comparison(1, 0) = 0.0;
        d2_3_1_comparison(1, 1) = 0.5;

        REQUIRE(bool(c2_0_0 == c2_0_0_comparison));
        REQUIRE(bool(c2_0_1 == c2_0_1_comparison));
        REQUIRE(bool(c2_1 == c2_1_comparison));
        REQUIRE(bool(c2_2 == c2_2_comparison));
        REQUIRE(bool(c2_3_0 == c2_3_0_comparison));
        REQUIRE(bool(c2_3_1 == c2_3_1_comparison));

        REQUIRE(bool(d2_0_0 == d2_0_0_comparison));
        REQUIRE(bool(d2_0_1 == d2_0_1_comparison));
        REQUIRE(bool(d2_1 == d2_1_comparison));
        REQUIRE(bool(d2_2 == d2_2_comparison));
        REQUIRE(bool(d2_3_0 == d2_3_0_comparison));
        REQUIRE(bool(d2_3_1 == d2_3_1_comparison));
    // }

    // SECTION("PerformKStep")
    // {
        multi_array<double, 2> ktau_comparison({2, 2});
        ktau_comparison = lr_sol.X;
        ktau_comparison(0, 0) -= tau * 0.25 * norm_2e;
        ktau_comparison(0, 1) += tau * 0.25 * norm_2e;
        ktau_comparison(1, 0) -= tau * 1.25 * norm_2e;
        ktau_comparison(1, 1) += tau * 0.75 * norm_2e;

        PerformKLStep<1>(sigma1, sigma2, lr_sol, blas, test_system, grid, partition1, partition2, w_x_dep, tau);
        REQUIRE(bool(lr_sol.X == ktau_comparison));
    // }
}