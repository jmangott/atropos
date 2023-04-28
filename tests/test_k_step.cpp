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
        Index dx1_comparison = 2;
        Index dx2_comparison = 2;
        
        REQUIRE(bool(grid.dx1 == dx1_comparison));
        REQUIRE(bool(grid.dx2 == dx2_comparison));
    // }

    multi_array<double, 2> k_comparison({2, 2});
    double norm = 1.0 / std::sqrt(2);
    double norm_2e = std::sqrt(2.0 / std::exp(1.0));

    k_comparison(0, 0) = norm_2e;
    k_comparison(0, 1) = 0.0;
    k_comparison(1, 0) = norm_2e;
    k_comparison(1, 1) = 0.0;

    tmp_x1 = lr_sol.X;
    blas.matmul(tmp_x1, lr_sol.S, lr_sol.X); // lr_sol.X contains now K

    REQUIRE(bool(lr_sol.X == k_comparison));

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
        vector<multi_array<double, 3>> c_coeff1_comparison(test_system.mu()), d_coeff1_comparison(test_system.mu());
        for (Index mu = 0; mu < test_system.mu(); mu++)
        {
            c_coeff1_comparison[mu].resize({partition1.dx_dep(mu), grid.r, grid.r});
            d_coeff1_comparison[mu].resize({partition1.dx_dep(mu), grid.r, grid.r});
        }

        CalculateCoefficientsKL<1>(c_coeff1, d_coeff1, sigma2, lr_sol, blas, test_system, grid, partition1, partition2, w_x_dep);

        c_coeff1_comparison[0](0, 0, 0) = 0.0;
        c_coeff1_comparison[0](0, 0, 1) = 0.0;
        c_coeff1_comparison[0](0, 1, 0) = 0.0;
        c_coeff1_comparison[0](0, 1, 1) = 0.0;

        c_coeff1_comparison[0](1, 0, 0) = 1.0;
        c_coeff1_comparison[0](1, 0, 1) = 0.0;
        c_coeff1_comparison[0](1, 1, 0) = 0.0;
        c_coeff1_comparison[0](1, 1, 1) = 1.0;

        c_coeff1_comparison[1](0, 0, 0) = 0.5;
        c_coeff1_comparison[1](0, 0, 1) = -0.5;
        c_coeff1_comparison[1](0, 1, 0) = 0.5;
        c_coeff1_comparison[1](0, 1, 1) = -0.5;

        c_coeff1_comparison[2](0, 0, 0) = 0.75;
        c_coeff1_comparison[2](0, 0, 1) = 0.25;
        c_coeff1_comparison[2](0, 1, 0) = 0.25;
        c_coeff1_comparison[2](0, 1, 1) = 0.75;

        c_coeff1_comparison[3](0, 0, 0) = 0.5;
        c_coeff1_comparison[3](0, 0, 1) = 0.5;
        c_coeff1_comparison[3](0, 1, 0) = -0.5;
        c_coeff1_comparison[3](0, 1, 1) = -0.5;

        c_coeff1_comparison[3](1, 0, 0) = 0.25;
        c_coeff1_comparison[3](1, 0, 1) = 0.25;
        c_coeff1_comparison[3](1, 1, 0) = -0.25;
        c_coeff1_comparison[3](1, 1, 1) = -0.25;

        d_coeff1_comparison[0](0, 0, 0) = 0.0;
        d_coeff1_comparison[0](0, 0, 1) = 0.0;
        d_coeff1_comparison[0](0, 1, 0) = 0.0;
        d_coeff1_comparison[0](0, 1, 1) = 0.0;

        d_coeff1_comparison[0](1, 0, 0) = 1.0;
        d_coeff1_comparison[0](1, 0, 1) = 0.0;
        d_coeff1_comparison[0](1, 1, 0) = 0.0;
        d_coeff1_comparison[0](1, 1, 1) = 1.0;

        d_coeff1_comparison[1](0, 0, 0) = 0.5;
        d_coeff1_comparison[1](0, 0, 1) = -0.5;
        d_coeff1_comparison[1](0, 1, 0) = -0.5;
        d_coeff1_comparison[1](0, 1, 1) = 0.5;

        d_coeff1_comparison[2](0, 0, 0) = 0.75;
        d_coeff1_comparison[2](0, 0, 1) = 0.25;
        d_coeff1_comparison[2](0, 1, 0) = 0.25;
        d_coeff1_comparison[2](0, 1, 1) = 0.75;

        d_coeff1_comparison[3](0, 0, 0) = 1.0;
        d_coeff1_comparison[3](0, 0, 1) = 0.0;
        d_coeff1_comparison[3](0, 1, 0) = 0.0;
        d_coeff1_comparison[3](0, 1, 1) = 1.0;

        d_coeff1_comparison[3](1, 0, 0) = 0.5;
        d_coeff1_comparison[3](1, 0, 1) = 0.0;
        d_coeff1_comparison[3](1, 1, 0) = 0.0;
        d_coeff1_comparison[3](1, 1, 1) = 0.5;

        for (Index mu = 0; mu < test_system.mu(); mu++)
        {
            REQUIRE(bool(c_coeff1[mu] == c_coeff1_comparison[mu]));
            REQUIRE(bool(d_coeff1[mu] == d_coeff1_comparison[mu]));
        }
    // }

    // SECTION("PerformKStep")
    // {
        multi_array<double, 2> ktau_comparison({2, 2});
        ktau_comparison(0, 0) =-tau * 0.25 * norm_2e;
        ktau_comparison(0, 1) = tau * 0.25 * norm_2e;
        ktau_comparison(1, 0) =-tau * 1.25 * norm_2e;
        ktau_comparison(1, 1) = tau * 0.75 * norm_2e;

        PerformKLStep<1>(tmp_x1, lr_sol.X, c_coeff1, d_coeff1, sigma1, blas, test_system, grid, partition1, w_x_dep, tau);
        REQUIRE(bool(tmp_x1 == ktau_comparison));
    // }
}