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

    InitializeTest(lr_sol, grid, ip_xx1, ip_xx2, blas, w_x_dep);

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

    multi_array<double, 1> state_vec2({1});
    multi_array<Index, 1> vec_index2({1});

    // SECTION("CalculateCoefficientsX")
    // {
        vector<multi_array<double, 3>> c_coeff2_comparison(test_system.mu()), d_coeff2_comparison(test_system.mu());
        for (Index mu = 0; mu < test_system.mu(); mu++)
        {
            c_coeff2_comparison[mu].resize({partition2.dx_dep(mu), grid.r, grid.r});
            d_coeff2_comparison[mu].resize({partition2.dx_dep(mu), grid.r, grid.r});
        }

        CalculateCoefficientsKL<2>(c_coeff2, d_coeff2, sigma1, lr_sol, blas, test_system, grid, partition2, partition1, w_x_dep);

        c_coeff2_comparison[0](0, 0, 0) = 0.5;
        c_coeff2_comparison[0](0, 0, 1) = -0.5;
        c_coeff2_comparison[0](0, 1, 0) = 0.5;
        c_coeff2_comparison[0](0, 1, 1) = -0.5;

        c_coeff2_comparison[1](0, 0, 0) = 0.0;
        c_coeff2_comparison[1](0, 0, 1) = 0.0;
        c_coeff2_comparison[1](0, 1, 0) = 0.0;
        c_coeff2_comparison[1](0, 1, 1) = 0.0;

        c_coeff2_comparison[1](1, 0, 0) = 1.0;
        c_coeff2_comparison[1](1, 0, 1) = 0.0;
        c_coeff2_comparison[1](1, 1, 0) = 0.0;
        c_coeff2_comparison[1](1, 1, 1) = 1.0;

        c_coeff2_comparison[2](0, 0, 0) = 0.5;
        c_coeff2_comparison[2](0, 0, 1) = 0.5;
        c_coeff2_comparison[2](0, 1, 0) = -0.5;
        c_coeff2_comparison[2](0, 1, 1) = -0.5;

        c_coeff2_comparison[2](1, 0, 0) = 0.25;
        c_coeff2_comparison[2](1, 0, 1) = 0.25;
        c_coeff2_comparison[2](1, 1, 0) = -0.25;
        c_coeff2_comparison[2](1, 1, 1) = -0.25;

        c_coeff2_comparison[3](0, 0, 0) = 0.75;
        c_coeff2_comparison[3](0, 0, 1) = 0.25;
        c_coeff2_comparison[3](0, 1, 0) = 0.25;
        c_coeff2_comparison[3](0, 1, 1) = 0.75;

        d_coeff2_comparison[0](0, 0, 0) = 0.5;
        d_coeff2_comparison[0](0, 0, 1) = -0.5;
        d_coeff2_comparison[0](0, 1, 0) = -0.5;
        d_coeff2_comparison[0](0, 1, 1) = 0.5;

        d_coeff2_comparison[1](0, 0, 0) = 0.0;
        d_coeff2_comparison[1](0, 0, 1) = 0.0;
        d_coeff2_comparison[1](0, 1, 0) = 0.0;
        d_coeff2_comparison[1](0, 1, 1) = 0.0;

        d_coeff2_comparison[1](1, 0, 0) = 1.0;
        d_coeff2_comparison[1](1, 0, 1) = 0.0;
        d_coeff2_comparison[1](1, 1, 0) = 0.0;
        d_coeff2_comparison[1](1, 1, 1) = 1.0;

        d_coeff2_comparison[2](0, 0, 0) = 1.0;
        d_coeff2_comparison[2](0, 0, 1) = 0.0;
        d_coeff2_comparison[2](0, 1, 0) = 0.0;
        d_coeff2_comparison[2](0, 1, 1) = 1.0;

        d_coeff2_comparison[2](1, 0, 0) = 0.5;
        d_coeff2_comparison[2](1, 0, 1) = 0.0;
        d_coeff2_comparison[2](1, 1, 0) = 0.0;
        d_coeff2_comparison[2](1, 1, 1) = 0.5;

        d_coeff2_comparison[3](0, 0, 0) = 0.75;
        d_coeff2_comparison[3](0, 0, 1) = 0.25;
        d_coeff2_comparison[3](0, 1, 0) = 0.25;
        d_coeff2_comparison[3](0, 1, 1) = 0.75;

        for (Index mu = 0; mu < test_system.mu(); mu++)
        {
            REQUIRE(bool(c_coeff2[mu] == c_coeff2_comparison[mu]));
            REQUIRE(bool(d_coeff2[mu] == d_coeff2_comparison[mu]));
        }
    // }

    // SECTION("PerformLStep")
    // {
        multi_array<double, 2> ltau_comparison({2, 2});
        ltau_comparison(0, 0) =-tau * 0.25 * norm_2e;
        ltau_comparison(0, 1) = tau * 0.25 * norm_2e;
        ltau_comparison(1, 0) =-tau * 1.25 * norm_2e;
        ltau_comparison(1, 1) = tau * 0.75 * norm_2e;

        PerformKLStep<2>(tmp_x2, lr_sol.V, c_coeff2, d_coeff2, sigma2, blas, test_system, grid, partition2, w_x_dep, tau);

        REQUIRE(bool(tmp_x2 == ltau_comparison));
    // }
}