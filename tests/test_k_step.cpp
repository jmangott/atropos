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
#include "test_parameters.hpp"

using std::cout;
using std::endl;
using std::vector;

TEST_CASE("k_step", "[k_step]")
{
    // Temporary objects for multiplication
    multi_array<double, 2> tmp_x({grid.dx1, r});

    InitializeTest<kM1, kM2>(lr_sol, grid, ip_xx1, ip_xx2, blas);

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

        REQUIRE(bool(grid.lim1 == lim1_comparison));
        REQUIRE(bool(grid.lim2 == lim2_comparison));
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

    tmp_x = lr_sol.X;
    blas.matmul(tmp_x, lr_sol.S, lr_sol.X); // lr_sol.X contains now K

    REQUIRE(bool(lr_sol.X == k_comparison));

    CalculateShiftAmount(sigma1, sigma2, test_system, grid);

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
        multi_array<double, 1> state_vec1({1}), state_vec2({1});
        multi_array<Index, 1> vec_index1({1}), vec_index2({1});

        // CASE 0: x_1 = 0.0, CASE 1: x_1 = 1.0
        for (Index i = 0; i < 2; i++)
        {
            vec_index1(0) = i;
            state_vec1 = VecIndexToState(vec_index1, grid.n1, grid.lim1);

            weight0 = CalculateWeightX2(vec_index1, test_system, grid, 0);
            weight1 = CalculateWeightX2(vec_index1, test_system, grid, 1);
            weight2 = CalculateWeightX2(vec_index1, test_system, grid, 2);
            weight3 = CalculateWeightX2(vec_index1, test_system, grid, 3);

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
        }
    // }

        // SECTION("CalculateCoefficientsX2")
        // {
        // Coefficients
        // CASE 0: x_1 = 0.0, CASE 1: x_1 = 1.0
        for (Index i = 0; i < 2; i++)
        {
            vec_index1(0) = i;
            state_vec1 = VecIndexToState(vec_index1, grid.n1, grid.lim1);

            multi_array<double, 2> c2_0({2, 2}), c2_1({2, 2}), c2_2({2, 2}), c2_3({2, 2});
            multi_array<double, 2> d2_0({2, 2}), d2_1({2, 2}), d2_2({2, 2}), d2_3({2, 2});

            multi_array<double, 2> c2_0_comparison({2, 2}), c2_1_comparison({2, 2}), c2_2_comparison({2, 2}), c2_3_comparison({2, 2});
            multi_array<double, 2> d2_0_comparison({2, 2}), d2_1_comparison({2, 2}), d2_2_comparison({2, 2}), d2_3_comparison({2, 2});

            multi_array<double, 2> xx2_shift0(lr_sol.V.shape());
            multi_array<double, 2> xx2_shift1(lr_sol.V.shape());
            multi_array<double, 2> xx2_shift2(lr_sol.V.shape());
            multi_array<double, 2> xx2_shift3(lr_sol.V.shape());

            ShiftMultiArrayRows(xx2_shift0, lr_sol.V, -sigma2[0]);
            ShiftMultiArrayRows(xx2_shift1, lr_sol.V, -sigma2[1]);
            ShiftMultiArrayRows(xx2_shift2, lr_sol.V, -sigma2[2]);
            ShiftMultiArrayRows(xx2_shift3, lr_sol.V, -sigma2[3]);

            CalculateCoefficientsX2(c2_0, d2_0, lr_sol, blas, xx2_shift0, vec_index1, test_system, grid, 0);
            CalculateCoefficientsX2(c2_1, d2_1, lr_sol, blas, xx2_shift1, vec_index1, test_system, grid, 1);
            CalculateCoefficientsX2(c2_2, d2_2, lr_sol, blas, xx2_shift2, vec_index1, test_system, grid, 2);
            CalculateCoefficientsX2(c2_3, d2_3, lr_sol, blas, xx2_shift3, vec_index1, test_system, grid, 3);

            c2_0_comparison(0, 0) = state_vec1(0);
            c2_0_comparison(0, 1) = 0.0;
            c2_0_comparison(1, 0) = 0.0;
            c2_0_comparison(1, 1) = state_vec1(0);

            c2_1_comparison(0, 0) = 0.5;
            c2_1_comparison(0, 1) =-0.5;
            c2_1_comparison(1, 0) = 0.5;
            c2_1_comparison(1, 1) =-0.5;

            c2_2_comparison(0, 0) = 0.75;
            c2_2_comparison(0, 1) = 0.25;
            c2_2_comparison(1, 0) = 0.25;
            c2_2_comparison(1, 1) = 0.75;

            c2_3_comparison(0, 0) = 1.0 / (1.0 + state_vec1(0));
            c2_3_comparison(0, 1) = 0.0;
            c2_3_comparison(1, 0) =-1.0 / (1.0 + state_vec1(0));
            c2_3_comparison(1, 1) = 0.0;

            d2_0_comparison = c2_0_comparison;

            d2_1_comparison(0, 0) = 0.5;
            d2_1_comparison(0, 1) =-0.5;
            d2_1_comparison(1, 0) =-0.5;
            d2_1_comparison(1, 1) = 0.5;

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
        }
    // }

    // SECTION("PerformKStep")
    // {
        multi_array<double, 2> ktau_comparison({2, 2});
        ktau_comparison = lr_sol.X;
        ktau_comparison(0, 0) += tau * norm_2e;
        ktau_comparison(1, 1) += tau * 0.5 * norm_2e;

        PerformKStep(sigma1, sigma2, lr_sol, blas, test_system, grid, tau);
        REQUIRE(bool(lr_sol.X == ktau_comparison));
    // }
}