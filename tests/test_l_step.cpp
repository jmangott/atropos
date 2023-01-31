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

TEST_CASE("l_step", "[l_step]")
{
    // Temporary objects for multiplication
    multi_array<double, 2> tmp_x2({grid.dx1, grid.r});

    InitializeTest(lr_sol, grid, ip_xx1, ip_xx2, blas, w_x);

    multi_array<double, 2> x1x2_comparison({2, 2}), l_comparison({2, 2});
    double norm = 1.0 / std::sqrt(2);
    double norm_2e = std::sqrt(2.0 / std::exp(1.0));

    x1x2_comparison(0, 0) = norm;
    x1x2_comparison(0, 1) = norm;
    x1x2_comparison(1, 0) = norm;
    x1x2_comparison(1, 1) = -norm;

    l_comparison(0, 0) = norm_2e;
    l_comparison(0, 1) = 0.0;
    l_comparison(1, 0) = norm_2e;
    l_comparison(1, 1) = 0.0;

    tmp_x2 = lr_sol.V;
    blas.matmul_transb(tmp_x2, lr_sol.S, lr_sol.V); // lr_sol.V contains now L

    REQUIRE(bool(lr_sol.V == l_comparison));

    CalculateShiftAmount(sigma1, sigma2, test_system, grid);

    // SECTION("CalculateWeightX")
    // {
        multi_array<double, 1> weight0_comparison({2}), weight1_comparison({2}), weight2_comparison({2}), weight3_comparison({2});
        multi_array<double, 1> state_vec2({1});
        multi_array<Index, 1> vec_index2({1});

        // CASE 0: x_2 = 0.0, CASE 1: x_2 = 1.0
        for (Index i = 0; i < 2; i++)
        {
            vec_index2(0) = i;
            state_vec2 = VecIndexToState(vec_index2, grid.n2, grid.lim2);

            weight0_comparison(0) = 0.0;
            weight0_comparison(1) = 1.0;
            weight1_comparison(0) = state_vec2(0);
            weight1_comparison(1) = state_vec2(0);
            weight2_comparison(0) = 1.0 / (1.0 + state_vec2(0));
            weight2_comparison(1) = 1.0 / (1.0 + state_vec2(0));
            weight3_comparison(0) = 1.0;
            weight3_comparison(1) = 0.5;

            REQUIRE(bool(w_x(0, i, 0) == weight0_comparison(0)));
            REQUIRE(bool(w_x(1, i, 0) == weight0_comparison(1)));
            REQUIRE(bool(w_x(0, i, 1) == weight1_comparison(0)));
            REQUIRE(bool(w_x(1, i, 1) == weight1_comparison(1)));
            REQUIRE(bool(w_x(0, i, 2) == weight2_comparison(0)));
            REQUIRE(bool(w_x(1, i, 2) == weight2_comparison(1)));
            REQUIRE(bool(w_x(0, i, 3) == weight3_comparison(0)));
            REQUIRE(bool(w_x(1, i, 3) == weight3_comparison(1)));
        }
    // }

        // SECTION("CalculateCoefficientsX")
        // {
        // Coefficients
        // CASE 0: x_2 = 0.0, CASE 1: x_2 = 1.0
        for (Index i = 0; i < 2; i++)
        {
            vec_index2(0) = i;
            state_vec2 = VecIndexToState(vec_index2, grid.n2, grid.lim2);

            multi_array<double, 2> c1_0({2, 2}), c1_1({2, 2}), c1_2({2, 2}), c1_3({2, 2});
            multi_array<double, 2> d1_0({2, 2}), d1_1({2, 2}), d1_2({2, 2}), d1_3({2, 2});

            multi_array<double, 2> c1_0_comparison({2, 2}), c1_1_comparison({2, 2}), c1_2_comparison({2, 2}), c1_3_comparison({2, 2});
            multi_array<double, 2> d1_0_comparison({2, 2}), d1_1_comparison({2, 2}), d1_2_comparison({2, 2}), d1_3_comparison({2, 2});

            multi_array<double, 2> xx1_shift0(lr_sol.X.shape());
            multi_array<double, 2> xx1_shift1(lr_sol.X.shape());
            multi_array<double, 2> xx1_shift2(lr_sol.X.shape());
            multi_array<double, 2> xx1_shift3(lr_sol.X.shape());

            ShiftMultiArrayRows(1, xx1_shift0, lr_sol.X, -sigma1[0], test_system.reactions[0]->minus_nu, grid, test_system);
            ShiftMultiArrayRows(1, xx1_shift1, lr_sol.X, -sigma1[1], test_system.reactions[1]->minus_nu, grid, test_system);
            ShiftMultiArrayRows(1, xx1_shift2, lr_sol.X, -sigma1[2], test_system.reactions[2]->minus_nu, grid, test_system);
            ShiftMultiArrayRows(1, xx1_shift3, lr_sol.X, -sigma1[3], test_system.reactions[3]->minus_nu, grid, test_system);

            CalculateCoefficientsX(2, c1_0, d1_0, lr_sol, blas, xx1_shift0, i, test_system, grid, 0, w_x);
            CalculateCoefficientsX(2, c1_1, d1_1, lr_sol, blas, xx1_shift1, i, test_system, grid, 1, w_x);
            CalculateCoefficientsX(2, c1_2, d1_2, lr_sol, blas, xx1_shift2, i, test_system, grid, 2, w_x);
            CalculateCoefficientsX(2, c1_3, d1_3, lr_sol, blas, xx1_shift3, i, test_system, grid, 3, w_x);

            c1_0_comparison(0, 0) = 0.5;
            c1_0_comparison(0, 1) =-0.5;
            c1_0_comparison(1, 0) = 0.5;
            c1_0_comparison(1, 1) =-0.5;

            c1_1_comparison(0, 0) = state_vec2(0);
            c1_1_comparison(0, 1) = 0.0;
            c1_1_comparison(1, 0) = 0.0;
            c1_1_comparison(1, 1) = state_vec2(0);

            c1_2_comparison(0, 0) = 0.5 / (1.0 + state_vec2(0));
            c1_2_comparison(0, 1) = 0.5 / (1.0 + state_vec2(0));
            c1_2_comparison(1, 0) =-0.5 / (1.0 + state_vec2(0));
            c1_2_comparison(1, 1) =-0.5 / (1.0 + state_vec2(0));

            c1_3_comparison(0, 0) = 0.75;
            c1_3_comparison(0, 1) = 0.25;
            c1_3_comparison(1, 0) = 0.25;
            c1_3_comparison(1, 1) = 0.75;

            d1_0_comparison(0, 0) = 0.5;
            d1_0_comparison(0, 1) =-0.5;
            d1_0_comparison(1, 0) =-0.5;
            d1_0_comparison(1, 1) = 0.5;

            d1_1_comparison = c1_1_comparison;

            d1_2_comparison(0, 0) = 1.0 / (1.0 + state_vec2(0));
            d1_2_comparison(0, 1) = 0.0;
            d1_2_comparison(1, 0) = 0.0;
            d1_2_comparison(1, 1) = 1.0 / (1.0 + state_vec2(0));

            d1_3_comparison = c1_3_comparison;

            REQUIRE(bool(c1_0 == c1_0_comparison));
            REQUIRE(bool(c1_1 == c1_1_comparison));
            REQUIRE(bool(c1_2 == c1_2_comparison));
            REQUIRE(bool(c1_3 == c1_3_comparison));

            REQUIRE(bool(d1_0 == d1_0_comparison));
            REQUIRE(bool(d1_1 == d1_1_comparison));
            REQUIRE(bool(d1_2 == d1_2_comparison));
            REQUIRE(bool(d1_3 == d1_3_comparison));
        }
    // }

    // SECTION("PerformLStep")
    // {
        multi_array<double, 2> ltau_comparison({2, 2});
        ltau_comparison = lr_sol.V;
        ltau_comparison(0, 0) -= tau * 0.25 * norm_2e;
        ltau_comparison(0, 1) += tau * 0.25 * norm_2e;
        ltau_comparison(1, 0) -= tau * 1.25 * norm_2e;
        ltau_comparison(1, 1) += tau * 0.75 * norm_2e;

        PerformKLStep<2>(sigma1, sigma2, lr_sol, blas, test_system, grid, partition2, w_x, tau);
        cout << ltau_comparison(0, 0) << " " << ltau_comparison(0, 1) << endl;
        cout << ltau_comparison(1, 0) << " " << ltau_comparison(1, 1) << endl;
        cout << endl;

        cout << lr_sol.V(0, 0) << " " << lr_sol.V(0, 1) << endl;
        cout << lr_sol.V(1, 0) << " " << lr_sol.V(1, 1) << endl;

        REQUIRE(bool(lr_sol.V == ltau_comparison));
    // }
}