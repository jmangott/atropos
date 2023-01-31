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
    double inv_sqrt_e = 1 / std::sqrt(std::exp(1.0));
    
    InitializeTest(lr_sol, grid, ip_xx1, ip_xx2, blas, w_x);
    CalculateShiftAmount(sigma1, sigma2, test_system, grid);

    // SECTION("CalculateCoefficientsB")
    // {
        multi_array<double, 3> b_coeff_vec_shift0({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec_shift1({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec_shift2({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec_shift3({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec0({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec1({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec2({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec3({grid.dx1, grid.r, grid.r});

        multi_array<double, 3> b_coeff_vec_shift0_comparison({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec_shift1_comparison({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec_shift2_comparison({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec_shift3_comparison({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec0_comparison({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec1_comparison({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec2_comparison({grid.dx1, grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec3_comparison({grid.dx1, grid.r, grid.r});

        CalculateCoefficientsB(b_coeff_vec_shift0, b_coeff_vec0, lr_sol, blas, test_system, grid, partition, 0, sigma2, w_x_dep);
        CalculateCoefficientsB(b_coeff_vec_shift1, b_coeff_vec1, lr_sol, blas, test_system, grid, partition, 1, sigma2, w_x_dep);
        CalculateCoefficientsB(b_coeff_vec_shift2, b_coeff_vec2, lr_sol, blas, test_system, grid, partition, 2, sigma2, w_x_dep);
        CalculateCoefficientsB(b_coeff_vec_shift3, b_coeff_vec3, lr_sol, blas, test_system, grid, partition, 3, sigma2, w_x_dep);

        b_coeff_vec_shift0_comparison(0, 0, 0) = 0.0;
        b_coeff_vec_shift0_comparison(0, 0, 1) = 0.0;
        b_coeff_vec_shift0_comparison(0, 1, 0) = 0.0;
        b_coeff_vec_shift0_comparison(0, 1, 1) = 0.0;
        b_coeff_vec_shift0_comparison(1, 0, 0) = 1.0;
        b_coeff_vec_shift0_comparison(1, 0, 1) = 0.0;
        b_coeff_vec_shift0_comparison(1, 1, 0) = 0.0;
        b_coeff_vec_shift0_comparison(1, 1, 1) = 1.0;

        b_coeff_vec_shift1_comparison(0, 0, 0) = 0.5;
        b_coeff_vec_shift1_comparison(0, 0, 1) =-0.5;
        b_coeff_vec_shift1_comparison(0, 1, 0) = 0.5;
        b_coeff_vec_shift1_comparison(0, 1, 1) =-0.5;
        b_coeff_vec_shift1_comparison(1, 0, 0) = 0.5;
        b_coeff_vec_shift1_comparison(1, 0, 1) =-0.5;
        b_coeff_vec_shift1_comparison(1, 1, 0) = 0.5;
        b_coeff_vec_shift1_comparison(1, 1, 1) =-0.5;

        b_coeff_vec_shift2_comparison(0, 0, 0) = 0.75;
        b_coeff_vec_shift2_comparison(0, 0, 1) = 0.25;
        b_coeff_vec_shift2_comparison(0, 1, 0) = 0.25;
        b_coeff_vec_shift2_comparison(0, 1, 1) = 0.75;
        b_coeff_vec_shift2_comparison(1, 0, 0) = 0.75;
        b_coeff_vec_shift2_comparison(1, 0, 1) = 0.25;
        b_coeff_vec_shift2_comparison(1, 1, 0) = 0.25;
        b_coeff_vec_shift2_comparison(1, 1, 1) = 0.75;

        b_coeff_vec_shift3_comparison(0, 0, 0) = 0.5;
        b_coeff_vec_shift3_comparison(0, 0, 1) = 0.5;
        b_coeff_vec_shift3_comparison(0, 1, 0) =-0.5;
        b_coeff_vec_shift3_comparison(0, 1, 1) =-0.5;
        b_coeff_vec_shift3_comparison(1, 0, 0) = 0.25;
        b_coeff_vec_shift3_comparison(1, 0, 1) = 0.25;
        b_coeff_vec_shift3_comparison(1, 1, 0) =-0.25;
        b_coeff_vec_shift3_comparison(1, 1, 1) =-0.25;

        REQUIRE(bool(b_coeff_vec_shift0 == b_coeff_vec_shift0_comparison));
        REQUIRE(bool(b_coeff_vec_shift1 == b_coeff_vec_shift1_comparison));
        REQUIRE(bool(b_coeff_vec_shift2 == b_coeff_vec_shift2_comparison));
        REQUIRE(bool(b_coeff_vec_shift3 == b_coeff_vec_shift3_comparison));
    // }

    // SECTION(CalculateCoefficientsS)
    // {
        multi_array<double, 4> e_coeff0({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> e_coeff1({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> e_coeff2({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> e_coeff3({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> f_coeff0({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> f_coeff1({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> f_coeff2({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> f_coeff3({grid.r, grid.r, grid.r, grid.r});

        CalculateCoefficientsS(e_coeff0, f_coeff0, b_coeff_vec_shift0, b_coeff_vec0, lr_sol, blas, test_system, grid, partition, 0, sigma1);
        CalculateCoefficientsS(e_coeff1, f_coeff1, b_coeff_vec_shift1, b_coeff_vec1, lr_sol, blas, test_system, grid, partition, 1, sigma1);
        CalculateCoefficientsS(e_coeff2, f_coeff2, b_coeff_vec_shift2, b_coeff_vec2, lr_sol, blas, test_system, grid, partition, 2, sigma1);
        CalculateCoefficientsS(e_coeff3, f_coeff3, b_coeff_vec_shift3, b_coeff_vec3, lr_sol, blas, test_system, grid, partition, 3, sigma1);

        multi_array<double, 4> e_coeff0_comparison({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> e_coeff1_comparison({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> e_coeff2_comparison({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> e_coeff3_comparison({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> f_coeff0_comparison({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> f_coeff1_comparison({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> f_coeff2_comparison({grid.r, grid.r, grid.r, grid.r});
        multi_array<double, 4> f_coeff3_comparison({grid.r, grid.r, grid.r, grid.r});

        // e_coeff0_comparison
        e_coeff0_comparison(0, 0, 0, 0) = 0.5;
        e_coeff0_comparison(0, 0, 1, 0) =-0.5;
        e_coeff0_comparison(1, 0, 0, 0) = 0.5;
        e_coeff0_comparison(1, 0, 1, 0) =-0.5;

        e_coeff0_comparison(0, 0, 0, 1) = 0.0;
        e_coeff0_comparison(0, 0, 1, 1) = 0.0;
        e_coeff0_comparison(1, 0, 0, 1) = 0.0;
        e_coeff0_comparison(1, 0, 1, 1) = 0.0;

        e_coeff0_comparison(0, 1, 0, 0) = e_coeff0_comparison(0, 0, 0, 1);
        e_coeff0_comparison(0, 1, 1, 0) = e_coeff0_comparison(0, 0, 1, 1);
        e_coeff0_comparison(1, 1, 0, 0) = e_coeff0_comparison(1, 0, 0, 1);
        e_coeff0_comparison(1, 1, 1, 0) = e_coeff0_comparison(1, 0, 1, 1);

        e_coeff0_comparison(0, 1, 0, 1) = e_coeff0_comparison(0, 0, 0, 0);
        e_coeff0_comparison(0, 1, 1, 1) = e_coeff0_comparison(0, 0, 1, 0);
        e_coeff0_comparison(1, 1, 0, 1) = e_coeff0_comparison(1, 0, 0, 0);
        e_coeff0_comparison(1, 1, 1, 1) = e_coeff0_comparison(1, 0, 1, 0);

        // f_coeff0_comparison
        f_coeff0_comparison = e_coeff0_comparison;
        f_coeff0_comparison(0, 0, 0, 0) = 0.5;
        f_coeff0_comparison(0, 0, 1, 0) = -0.5;
        f_coeff0_comparison(1, 0, 0, 0) = -0.5;
        f_coeff0_comparison(1, 0, 1, 0) = 0.5;

        f_coeff0_comparison(0, 1, 0, 1) = f_coeff0_comparison(0, 0, 0, 0);
        f_coeff0_comparison(0, 1, 1, 1) = f_coeff0_comparison(0, 0, 1, 0);
        f_coeff0_comparison(1, 1, 0, 1) = f_coeff0_comparison(1, 0, 0, 0);
        f_coeff0_comparison(1, 1, 1, 1) = f_coeff0_comparison(1, 0, 1, 0);

        // e_coeff1_comparison
        e_coeff1_comparison(0, 0, 0, 0) = 0.5;
        e_coeff1_comparison(0, 0, 1, 0) = 0.0;
        e_coeff1_comparison(1, 0, 0, 0) = 0.0;
        e_coeff1_comparison(1, 0, 1, 0) = 0.5;

        e_coeff1_comparison(0, 0, 0, 1) = -e_coeff1_comparison(0, 0, 0, 0);
        e_coeff1_comparison(0, 0, 1, 1) = -e_coeff1_comparison(0, 0, 1, 0);
        e_coeff1_comparison(1, 0, 0, 1) = -e_coeff1_comparison(1, 0, 0, 0);
        e_coeff1_comparison(1, 0, 1, 1) = -e_coeff1_comparison(1, 0, 1, 0);

        e_coeff1_comparison(0, 1, 0, 1) = -e_coeff1_comparison(0, 0, 0, 0);
        e_coeff1_comparison(0, 1, 1, 1) = -e_coeff1_comparison(0, 0, 1, 0);
        e_coeff1_comparison(1, 1, 0, 1) = -e_coeff1_comparison(1, 0, 0, 0);
        e_coeff1_comparison(1, 1, 1, 1) = -e_coeff1_comparison(1, 0, 1, 0);

        e_coeff1_comparison(0, 1, 0, 0) = e_coeff1_comparison(0, 0, 0, 0);
        e_coeff1_comparison(0, 1, 1, 0) = e_coeff1_comparison(0, 0, 1, 0);
        e_coeff1_comparison(1, 1, 0, 0) = e_coeff1_comparison(1, 0, 0, 0);
        e_coeff1_comparison(1, 1, 1, 0) = e_coeff1_comparison(1, 0, 1, 0);

        // f_coeff1_comparison
        f_coeff1_comparison = e_coeff1_comparison;
        f_coeff1_comparison(0, 1, 0, 0) = -f_coeff1_comparison(0, 0, 0, 0);
        f_coeff1_comparison(0, 1, 1, 0) = -f_coeff1_comparison(0, 0, 1, 0);
        f_coeff1_comparison(1, 1, 0, 0) = -f_coeff1_comparison(1, 0, 0, 0);
        f_coeff1_comparison(1, 1, 1, 0) = -f_coeff1_comparison(1, 0, 1, 0);

        f_coeff1_comparison(0, 1, 0, 1) = f_coeff1_comparison(0, 0, 0, 0);
        f_coeff1_comparison(0, 1, 1, 1) = f_coeff1_comparison(0, 0, 1, 0);
        f_coeff1_comparison(1, 1, 0, 1) = f_coeff1_comparison(1, 0, 0, 0);
        f_coeff1_comparison(1, 1, 1, 1) = f_coeff1_comparison(1, 0, 1, 0);

        // e_coeff2_comparison
        e_coeff2_comparison(0, 0, 0, 0) = 0.375;
        e_coeff2_comparison(0, 0, 1, 0) = 0.375;
        e_coeff2_comparison(1, 0, 0, 0) =-0.375;
        e_coeff2_comparison(1, 0, 1, 0) =-0.375;

        e_coeff2_comparison(0, 0, 0, 1) = 0.125;
        e_coeff2_comparison(0, 0, 1, 1) = 0.125;
        e_coeff2_comparison(1, 0, 0, 1) =-0.125;
        e_coeff2_comparison(1, 0, 1, 1) =-0.125;

        e_coeff2_comparison(0, 1, 0, 0) = e_coeff2_comparison(0, 0, 0, 1);
        e_coeff2_comparison(0, 1, 1, 0) = e_coeff2_comparison(0, 0, 1, 1);
        e_coeff2_comparison(1, 1, 0, 0) = e_coeff2_comparison(1, 0, 0, 1);
        e_coeff2_comparison(1, 1, 1, 0) = e_coeff2_comparison(1, 0, 1, 1);

        e_coeff2_comparison(0, 1, 0, 1) = e_coeff2_comparison(0, 0, 0, 0);
        e_coeff2_comparison(0, 1, 1, 1) = e_coeff2_comparison(0, 0, 1, 0);
        e_coeff2_comparison(1, 1, 0, 1) = e_coeff2_comparison(1, 0, 0, 0);
        e_coeff2_comparison(1, 1, 1, 1) = e_coeff2_comparison(1, 0, 1, 0);

        // f_coeff2_comparison
        f_coeff2_comparison(0, 0, 0, 0) = 0.75;
        f_coeff2_comparison(0, 0, 1, 0) = 0.0;
        f_coeff2_comparison(1, 0, 0, 0) = 0.0;
        f_coeff2_comparison(1, 0, 1, 0) = 0.75;

        f_coeff2_comparison(0, 0, 0, 1) = 0.25;
        f_coeff2_comparison(0, 0, 1, 1) = 0.0;
        f_coeff2_comparison(1, 0, 0, 1) = 0.0;
        f_coeff2_comparison(1, 0, 1, 1) = 0.25;

        f_coeff2_comparison(0, 1, 0, 0) = f_coeff2_comparison(0, 0, 0, 1);
        f_coeff2_comparison(0, 1, 1, 0) = f_coeff2_comparison(0, 0, 1, 1);
        f_coeff2_comparison(1, 1, 0, 0) = f_coeff2_comparison(1, 0, 0, 1);
        f_coeff2_comparison(1, 1, 1, 0) = f_coeff2_comparison(1, 0, 1, 1);

        f_coeff2_comparison(0, 1, 0, 1) = f_coeff2_comparison(0, 0, 0, 0);
        f_coeff2_comparison(0, 1, 1, 1) = f_coeff2_comparison(0, 0, 1, 0);
        f_coeff2_comparison(1, 1, 0, 1) = f_coeff2_comparison(1, 0, 0, 0);
        f_coeff2_comparison(1, 1, 1, 1) = f_coeff2_comparison(1, 0, 1, 0);

        // e_coeff3_comparison
        e_coeff3_comparison(0, 0, 0, 0) = 0.375;
        e_coeff3_comparison(0, 0, 1, 0) = 0.125;
        e_coeff3_comparison(1, 0, 0, 0) = 0.125;
        e_coeff3_comparison(1, 0, 1, 0) = 0.375;

        e_coeff3_comparison(0, 0, 0, 1) = e_coeff3_comparison(0, 0, 0, 0);
        e_coeff3_comparison(0, 0, 1, 1) = e_coeff3_comparison(0, 0, 1, 0);
        e_coeff3_comparison(1, 0, 0, 1) = e_coeff3_comparison(1, 0, 0, 0);
        e_coeff3_comparison(1, 0, 1, 1) = e_coeff3_comparison(1, 0, 1, 0);

        e_coeff3_comparison(0, 1, 0, 1) = -e_coeff3_comparison(0, 0, 0, 0);
        e_coeff3_comparison(0, 1, 1, 1) = -e_coeff3_comparison(0, 0, 1, 0);
        e_coeff3_comparison(1, 1, 0, 1) = -e_coeff3_comparison(1, 0, 0, 0);
        e_coeff3_comparison(1, 1, 1, 1) = -e_coeff3_comparison(1, 0, 1, 0);

        e_coeff3_comparison(0, 1, 0, 0) = -e_coeff3_comparison(0, 0, 0, 0);
        e_coeff3_comparison(0, 1, 1, 0) = -e_coeff3_comparison(0, 0, 1, 0);
        e_coeff3_comparison(1, 1, 0, 0) = -e_coeff3_comparison(1, 0, 0, 0);
        e_coeff3_comparison(1, 1, 1, 0) = -e_coeff3_comparison(1, 0, 1, 0);

        // f_coeff2_comparison
        f_coeff3_comparison(0, 0, 0, 0) = 0.75;
        f_coeff3_comparison(0, 0, 1, 0) = 0.25;
        f_coeff3_comparison(1, 0, 0, 0) = 0.25;
        f_coeff3_comparison(1, 0, 1, 0) = 0.75;

        f_coeff3_comparison(0, 0, 0, 1) = 0.0;
        f_coeff3_comparison(0, 0, 1, 1) = 0.0;
        f_coeff3_comparison(1, 0, 0, 1) = 0.0;
        f_coeff3_comparison(1, 0, 1, 1) = 0.0;

        f_coeff3_comparison(0, 1, 0, 0) = f_coeff3_comparison(0, 0, 0, 1);
        f_coeff3_comparison(0, 1, 1, 0) = f_coeff3_comparison(0, 0, 1, 1);
        f_coeff3_comparison(1, 1, 0, 0) = f_coeff3_comparison(1, 0, 0, 1);
        f_coeff3_comparison(1, 1, 1, 0) = f_coeff3_comparison(1, 0, 1, 1);

        f_coeff3_comparison(0, 1, 0, 1) = f_coeff3_comparison(0, 0, 0, 0);
        f_coeff3_comparison(0, 1, 1, 1) = f_coeff3_comparison(0, 0, 1, 0);
        f_coeff3_comparison(1, 1, 0, 1) = f_coeff3_comparison(1, 0, 0, 0);
        f_coeff3_comparison(1, 1, 1, 1) = f_coeff3_comparison(1, 0, 1, 0);

        REQUIRE(bool(e_coeff0 == e_coeff0_comparison));
        REQUIRE(bool(f_coeff0 == f_coeff0_comparison));
        REQUIRE(bool(e_coeff1 == e_coeff1_comparison));
        REQUIRE(bool(f_coeff1 == f_coeff1_comparison));
        REQUIRE(bool(e_coeff2 == e_coeff2_comparison));
        REQUIRE(bool(f_coeff2 == f_coeff2_comparison));
        REQUIRE(bool(e_coeff3 == e_coeff3_comparison));
        REQUIRE(bool(f_coeff3 == f_coeff3_comparison));
    // }

    // SECTION("PerformSStep")
    // {
        multi_array<double, 2> stau_comparison({2, 2});
        stau_comparison = lr_sol.S;
        stau_comparison(0, 0) += tau * 1.5 * inv_sqrt_e;
        stau_comparison(0, 1) -= tau * inv_sqrt_e;
        stau_comparison(1, 0) -= tau * inv_sqrt_e;
        stau_comparison(1, 1) += tau * 0.5 * inv_sqrt_e;

        PerformSStep(sigma1, sigma2, lr_sol, blas, test_system, grid, partition, w_x_dep, tau);
        REQUIRE(bool(lr_sol.S == stau_comparison));
        // }
}