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
#include "s_step_functions.hpp"
#include "test_parameters.hpp"

using std::cout;
using std::endl;
using std::vector;

TEST_CASE("s_step", "[s_step]")
{
    // Temporary objects for multiplication
    multi_array<double, 2> tmp_s({grid.r, grid.r});
    double inv_sqrt_e = 1 / std::sqrt(std::exp(1.0));
    
    InitializeTest(lr_sol, grid, ip_xx1, ip_xx2, blas, w_x_dep);

    // SECTION("CalculateCoefficientsB")
    // {
        multi_array<double, 3> b_coeff_vec_shift0({partition1.dx_dep(0), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec_shift1({partition1.dx_dep(1), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec_shift2({partition1.dx_dep(2), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec_shift3({partition1.dx_dep(3), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec0({partition1.dx_dep(0), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec1({partition1.dx_dep(1), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec2({partition1.dx_dep(2), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec3({partition1.dx_dep(3), grid.r, grid.r});

        multi_array<double, 3> b_coeff_vec_shift0_comparison({partition1.dx_dep(0), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec_shift1_comparison({partition1.dx_dep(1), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec_shift2_comparison({partition1.dx_dep(2), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec_shift3_comparison({partition1.dx_dep(3), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec0_comparison({partition1.dx_dep(0), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec1_comparison({partition1.dx_dep(1), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec2_comparison({partition1.dx_dep(2), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec3_comparison({partition1.dx_dep(3), grid.r, grid.r});

        CalculateCoefficientsB(b_coeff_vec_shift0, b_coeff_vec0, lr_sol, blas, test_system, grid, partition1, partition2, 0, sigma2, w_x_dep);
        CalculateCoefficientsB(b_coeff_vec_shift1, b_coeff_vec1, lr_sol, blas, test_system, grid, partition1, partition2, 1, sigma2, w_x_dep);
        CalculateCoefficientsB(b_coeff_vec_shift2, b_coeff_vec2, lr_sol, blas, test_system, grid, partition1, partition2, 2, sigma2, w_x_dep);
        CalculateCoefficientsB(b_coeff_vec_shift3, b_coeff_vec3, lr_sol, blas, test_system, grid, partition1, partition2, 3, sigma2, w_x_dep);

        b_coeff_vec_shift0_comparison(0, 0, 0) = 0.0;
        b_coeff_vec_shift0_comparison(0, 0, 1) = 0.0;
        b_coeff_vec_shift0_comparison(0, 1, 0) = 0.0;
        b_coeff_vec_shift0_comparison(0, 1, 1) = 0.0;
        b_coeff_vec_shift0_comparison(1, 0, 0) = 1.0;
        b_coeff_vec_shift0_comparison(1, 0, 1) = 0.0;
        b_coeff_vec_shift0_comparison(1, 1, 0) = 0.0;
        b_coeff_vec_shift0_comparison(1, 1, 1) = 1.0;

        b_coeff_vec_shift1_comparison(0, 0, 0) = 0.5;
        b_coeff_vec_shift1_comparison(0, 0, 1) = -0.5;
        b_coeff_vec_shift1_comparison(0, 1, 0) = 0.5;
        b_coeff_vec_shift1_comparison(0, 1, 1) = -0.5;

        b_coeff_vec_shift2_comparison(0, 0, 0) = 0.75;
        b_coeff_vec_shift2_comparison(0, 0, 1) = 0.25;
        b_coeff_vec_shift2_comparison(0, 1, 0) = 0.25;
        b_coeff_vec_shift2_comparison(0, 1, 1) = 0.75;

        b_coeff_vec_shift3_comparison(0, 0, 0) = 0.5;
        b_coeff_vec_shift3_comparison(0, 0, 1) = 0.5;
        b_coeff_vec_shift3_comparison(0, 1, 0) =-0.5;
        b_coeff_vec_shift3_comparison(0, 1, 1) =-0.5;
        b_coeff_vec_shift3_comparison(1, 0, 0) = 0.25;
        b_coeff_vec_shift3_comparison(1, 0, 1) = 0.25;
        b_coeff_vec_shift3_comparison(1, 1, 0) = -0.25;
        b_coeff_vec_shift3_comparison(1, 1, 1) = -0.25;
        // TODO: What about b_coeff_vec?
        REQUIRE(bool(b_coeff_vec_shift0 == b_coeff_vec_shift0_comparison));
        REQUIRE(bool(b_coeff_vec_shift1 == b_coeff_vec_shift1_comparison));
        REQUIRE(bool(b_coeff_vec_shift2 == b_coeff_vec_shift2_comparison));
        REQUIRE(bool(b_coeff_vec_shift3 == b_coeff_vec_shift3_comparison));
    // }

    // SECTION(CalculateCoefficientsS)
    // {
        multi_array<double, 5> e_coeff({test_system.mu(), grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 5> f_coeff({test_system.mu(), grid.r, grid.r, grid.r, grid.r});

        CalculateCoefficientsS(e_coeff, f_coeff, sigma1, sigma2, lr_sol, blas, test_system, grid, partition1, partition2, w_x_dep);

        multi_array<double, 5> e_coeff_comparison({test_system.mu(), grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 5> f_coeff_comparison({test_system.mu(), grid.r, grid.r, grid.r, grid.r});

        // e_coeff_comparison, mu = 0
        e_coeff_comparison(0, 0, 0, 0, 0) = 0.5;
        e_coeff_comparison(0, 0, 0, 1, 0) =-0.5;
        e_coeff_comparison(0, 1, 0, 0, 0) = 0.5;
        e_coeff_comparison(0, 1, 0, 1, 0) =-0.5;

        e_coeff_comparison(0, 0, 0, 0, 1) = 0.0;
        e_coeff_comparison(0, 0, 0, 1, 1) = 0.0;
        e_coeff_comparison(0, 1, 0, 0, 1) = 0.0;
        e_coeff_comparison(0, 1, 0, 1, 1) = 0.0;

        e_coeff_comparison(0, 0, 1, 0, 0) = e_coeff_comparison(0, 0, 0, 0, 1);
        e_coeff_comparison(0, 0, 1, 1, 0) = e_coeff_comparison(0, 0, 0, 1, 1);
        e_coeff_comparison(0, 1, 1, 0, 0) = e_coeff_comparison(0, 1, 0, 0, 1);
        e_coeff_comparison(0, 1, 1, 1, 0) = e_coeff_comparison(0, 1, 0, 1, 1);

        e_coeff_comparison(0, 0, 1, 0, 1) = e_coeff_comparison(0, 0, 0, 0, 0);
        e_coeff_comparison(0, 0, 1, 1, 1) = e_coeff_comparison(0, 0, 0, 1, 0);
        e_coeff_comparison(0, 1, 1, 0, 1) = e_coeff_comparison(0, 1, 0, 0, 0);
        e_coeff_comparison(0, 1, 1, 1, 1) = e_coeff_comparison(0, 1, 0, 1, 0);

        // f_coeff_comparison, mu = 0
        f_coeff_comparison(0, 0, 0, 0, 0) = 0.5;
        f_coeff_comparison(0, 0, 0, 1, 0) =-0.5;
        f_coeff_comparison(0, 1, 0, 0, 0) =-0.5;
        f_coeff_comparison(0, 1, 0, 1, 0) = 0.5;

        f_coeff_comparison(0, 0, 0, 0, 1) = 0.0;
        f_coeff_comparison(0, 0, 0, 1, 1) = 0.0;
        f_coeff_comparison(0, 1, 0, 0, 1) = 0.0;
        f_coeff_comparison(0, 1, 0, 1, 1) = 0.0;

        f_coeff_comparison(0, 0, 1, 0, 0) = f_coeff_comparison(0, 0, 0, 0, 1);
        f_coeff_comparison(0, 0, 1, 1, 0) = f_coeff_comparison(0, 0, 0, 1, 1);
        f_coeff_comparison(0, 1, 1, 0, 0) = f_coeff_comparison(0, 1, 0, 0, 1);
        f_coeff_comparison(0, 1, 1, 1, 0) = f_coeff_comparison(0, 1, 0, 1, 1);

        f_coeff_comparison(0, 0, 1, 0, 1) = f_coeff_comparison(0, 0, 0, 0, 0);
        f_coeff_comparison(0, 0, 1, 1, 1) = f_coeff_comparison(0, 0, 0, 1, 0);
        f_coeff_comparison(0, 1, 1, 0, 1) = f_coeff_comparison(0, 1, 0, 0, 0);
        f_coeff_comparison(0, 1, 1, 1, 1) = f_coeff_comparison(0, 1, 0, 1, 0);

        // e_coeff_comparison, mu = 1
        e_coeff_comparison(1, 0, 0, 0, 0) = 0.5;
        e_coeff_comparison(1, 0, 0, 1, 0) = 0.0;
        e_coeff_comparison(1, 1, 0, 0, 0) = 0.0;
        e_coeff_comparison(1, 1, 0, 1, 0) = 0.5;

        e_coeff_comparison(1, 0, 0, 0, 1) = -e_coeff_comparison(1, 0, 0, 0, 0);
        e_coeff_comparison(1, 0, 0, 1, 1) = -e_coeff_comparison(1, 0, 0, 1, 0);
        e_coeff_comparison(1, 1, 0, 0, 1) = -e_coeff_comparison(1, 1, 0, 0, 0);
        e_coeff_comparison(1, 1, 0, 1, 1) = -e_coeff_comparison(1, 1, 0, 1, 0);

        e_coeff_comparison(1, 0, 1, 0, 0) = e_coeff_comparison(1, 0, 0, 0, 0);
        e_coeff_comparison(1, 0, 1, 1, 0) = e_coeff_comparison(1, 0, 0, 1, 0);
        e_coeff_comparison(1, 1, 1, 0, 0) = e_coeff_comparison(1, 1, 0, 0, 0);
        e_coeff_comparison(1, 1, 1, 1, 0) = e_coeff_comparison(1, 1, 0, 1, 0);

        e_coeff_comparison(1, 0, 1, 0, 1) = -e_coeff_comparison(1, 0, 0, 0, 0);
        e_coeff_comparison(1, 0, 1, 1, 1) = -e_coeff_comparison(1, 0, 0, 1, 0);
        e_coeff_comparison(1, 1, 1, 0, 1) = -e_coeff_comparison(1, 1, 0, 0, 0);
        e_coeff_comparison(1, 1, 1, 1, 1) = -e_coeff_comparison(1, 1, 0, 1, 0);

        // f_coeff_comparison, mu = 1
        f_coeff_comparison(1, 0, 0, 0, 0) = 0.5;
        f_coeff_comparison(1, 0, 0, 1, 0) = 0.0;
        f_coeff_comparison(1, 1, 0, 0, 0) = 0.0;
        f_coeff_comparison(1, 1, 0, 1, 0) = 0.5;

        f_coeff_comparison(1, 0, 0, 0, 1) = -f_coeff_comparison(1, 0, 0, 0, 0);
        f_coeff_comparison(1, 0, 0, 1, 1) = -f_coeff_comparison(1, 0, 0, 1, 0);
        f_coeff_comparison(1, 1, 0, 0, 1) = -f_coeff_comparison(1, 1, 0, 0, 0);
        f_coeff_comparison(1, 1, 0, 1, 1) = -f_coeff_comparison(1, 1, 0, 1, 0);

        f_coeff_comparison(1, 0, 1, 0, 0) = -f_coeff_comparison(1, 0, 0, 0, 0);
        f_coeff_comparison(1, 0, 1, 1, 0) = -f_coeff_comparison(1, 0, 0, 1, 0);
        f_coeff_comparison(1, 1, 1, 0, 0) = -f_coeff_comparison(1, 1, 0, 0, 0);
        f_coeff_comparison(1, 1, 1, 1, 0) = -f_coeff_comparison(1, 1, 0, 1, 0);

        f_coeff_comparison(1, 0, 1, 0, 1) = f_coeff_comparison(1, 0, 0, 0, 0);
        f_coeff_comparison(1, 0, 1, 1, 1) = f_coeff_comparison(1, 0, 0, 1, 0);
        f_coeff_comparison(1, 1, 1, 0, 1) = f_coeff_comparison(1, 1, 0, 0, 0);
        f_coeff_comparison(1, 1, 1, 1, 1) = f_coeff_comparison(1, 1, 0, 1, 0);

        // e_coeff_comparison, mu = 2
        e_coeff_comparison(2, 0, 0, 0, 0) = 0.375;
        e_coeff_comparison(2, 0, 0, 1, 0) = 0.375;
        e_coeff_comparison(2, 1, 0, 0, 0) =-0.375;
        e_coeff_comparison(2, 1, 0, 1, 0) =-0.375;

        e_coeff_comparison(2, 0, 0, 0, 1) = 0.125;
        e_coeff_comparison(2, 0, 0, 1, 1) = 0.125;
        e_coeff_comparison(2, 1, 0, 0, 1) =-0.125;
        e_coeff_comparison(2, 1, 0, 1, 1) =-0.125;

        e_coeff_comparison(2, 0, 1, 0, 0) = e_coeff_comparison(2, 0, 0, 0, 1);
        e_coeff_comparison(2, 0, 1, 1, 0) = e_coeff_comparison(2, 0, 0, 1, 1);
        e_coeff_comparison(2, 1, 1, 0, 0) = e_coeff_comparison(2, 1, 0, 0, 1);
        e_coeff_comparison(2, 1, 1, 1, 0) = e_coeff_comparison(2, 1, 0, 1, 1);

        e_coeff_comparison(2, 0, 1, 0, 1) = e_coeff_comparison(2, 0, 0, 0, 0);
        e_coeff_comparison(2, 0, 1, 1, 1) = e_coeff_comparison(2, 0, 0, 1, 0);
        e_coeff_comparison(2, 1, 1, 0, 1) = e_coeff_comparison(2, 1, 0, 0, 0);
        e_coeff_comparison(2, 1, 1, 1, 1) = e_coeff_comparison(2, 1, 0, 1, 0);

        // f_coeff_comparison, mu = 2
        f_coeff_comparison(2, 0, 0, 0, 0) = 0.75;
        f_coeff_comparison(2, 0, 0, 1, 0) = 0.0;
        f_coeff_comparison(2, 1, 0, 0, 0) = 0.0;
        f_coeff_comparison(2, 1, 0, 1, 0) = 0.75;

        f_coeff_comparison(2, 0, 0, 0, 1) = 0.25;
        f_coeff_comparison(2, 0, 0, 1, 1) = 0.0;
        f_coeff_comparison(2, 1, 0, 0, 1) = 0.0;
        f_coeff_comparison(2, 1, 0, 1, 1) = 0.25;

        f_coeff_comparison(2, 0, 1, 0, 0) = f_coeff_comparison(2, 0, 0, 0, 1);
        f_coeff_comparison(2, 0, 1, 1, 0) = f_coeff_comparison(2, 0, 0, 1, 1);
        f_coeff_comparison(2, 1, 1, 0, 0) = f_coeff_comparison(2, 1, 0, 0, 1);
        f_coeff_comparison(2, 1, 1, 1, 0) = f_coeff_comparison(2, 1, 0, 1, 1);

        f_coeff_comparison(2, 0, 1, 0, 1) = f_coeff_comparison(2, 0, 0, 0, 0);
        f_coeff_comparison(2, 0, 1, 1, 1) = f_coeff_comparison(2, 0, 0, 1, 0);
        f_coeff_comparison(2, 1, 1, 0, 1) = f_coeff_comparison(2, 1, 0, 0, 0);
        f_coeff_comparison(2, 1, 1, 1, 1) = f_coeff_comparison(2, 1, 0, 1, 0);

        // e_coeff_comparison, mu = 3
        e_coeff_comparison(3, 0, 0, 0, 0) = 0.375;
        e_coeff_comparison(3, 0, 0, 1, 0) = 0.125;
        e_coeff_comparison(3, 1, 0, 0, 0) = 0.125;
        e_coeff_comparison(3, 1, 0, 1, 0) = 0.375;

        e_coeff_comparison(3, 0, 0, 0, 1) = e_coeff_comparison(3, 0, 0, 0, 0);
        e_coeff_comparison(3, 0, 0, 1, 1) = e_coeff_comparison(3, 0, 0, 1, 0);
        e_coeff_comparison(3, 1, 0, 0, 1) = e_coeff_comparison(3, 1, 0, 0, 0);
        e_coeff_comparison(3, 1, 0, 1, 1) = e_coeff_comparison(3, 1, 0, 1, 0);

        e_coeff_comparison(3, 0, 1, 0, 1) = -e_coeff_comparison(3, 0, 0, 0, 0);
        e_coeff_comparison(3, 0, 1, 1, 1) = -e_coeff_comparison(3, 0, 0, 1, 0);
        e_coeff_comparison(3, 1, 1, 0, 1) = -e_coeff_comparison(3, 1, 0, 0, 0);
        e_coeff_comparison(3, 1, 1, 1, 1) = -e_coeff_comparison(3, 1, 0, 1, 0);

        e_coeff_comparison(3, 0, 1, 0, 0) = -e_coeff_comparison(3, 0, 0, 0, 0);
        e_coeff_comparison(3, 0, 1, 1, 0) = -e_coeff_comparison(3, 0, 0, 1, 0);
        e_coeff_comparison(3, 1, 1, 0, 0) = -e_coeff_comparison(3, 1, 0, 0, 0);
        e_coeff_comparison(3, 1, 1, 1, 0) = -e_coeff_comparison(3, 1, 0, 1, 0);

        // f_coeff_comparison, mu = 3
        f_coeff_comparison(3, 0, 0, 0, 0) = 0.75;
        f_coeff_comparison(3, 0, 0, 1, 0) = 0.25;
        f_coeff_comparison(3, 1, 0, 0, 0) = 0.25;
        f_coeff_comparison(3, 1, 0, 1, 0) = 0.75;

        f_coeff_comparison(3, 0, 0, 0, 1) = 0.0;
        f_coeff_comparison(3, 0, 0, 1, 1) = 0.0;
        f_coeff_comparison(3, 1, 0, 0, 1) = 0.0;
        f_coeff_comparison(3, 1, 0, 1, 1) = 0.0;

        f_coeff_comparison(3, 0, 1, 0, 0) = f_coeff_comparison(3, 0, 0, 0, 1);
        f_coeff_comparison(3, 0, 1, 1, 0) = f_coeff_comparison(3, 0, 0, 1, 1);
        f_coeff_comparison(3, 1, 1, 0, 0) = f_coeff_comparison(3, 1, 0, 0, 1);
        f_coeff_comparison(3, 1, 1, 1, 0) = f_coeff_comparison(3, 1, 0, 1, 1);

        f_coeff_comparison(3, 0, 1, 0, 1) = f_coeff_comparison(3, 0, 0, 0, 0);
        f_coeff_comparison(3, 0, 1, 1, 1) = f_coeff_comparison(3, 0, 0, 1, 0);
        f_coeff_comparison(3, 1, 1, 0, 1) = f_coeff_comparison(3, 1, 0, 0, 0);
        f_coeff_comparison(3, 1, 1, 1, 1) = f_coeff_comparison(3, 1, 0, 1, 0);

        REQUIRE(bool(e_coeff == e_coeff_comparison));
        REQUIRE(bool(f_coeff == f_coeff_comparison));
    // }

    // SECTION("PerformSStep")
    // {
        multi_array<double, 2> stau_comparison({2, 2});
        stau_comparison(0, 0) = tau * 1.5 * inv_sqrt_e;
        stau_comparison(0, 1) =-tau * inv_sqrt_e;
        stau_comparison(1, 0) =-tau * inv_sqrt_e;
        stau_comparison(1, 1) = tau * 0.5 * inv_sqrt_e;

        PerformSStep(tmp_s, lr_sol.S, e_coeff, f_coeff, sigma1, sigma2, blas, test_system, grid, partition1, partition2, w_x_dep, -tau);

        REQUIRE(bool(tmp_s == stau_comparison));
    // }
}