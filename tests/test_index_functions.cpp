#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"

TEST_CASE("index_functions", "[index_functions]")
{
    vector<string> test_names = {"S1, S2, S3"};
    mysys test_system(test_names);
    myreact test_react0(
        {-1, -1, 0}, {0, 1}, [](std::vector<double> y)
        { return y[0] * y[1]; },
        test_system);
    myreact test_react1(
        {0, -1, -1}, {1, 2}, [](std::vector<double> y)
        { return y[1] * y[2]; },
        test_system);
    myreact test_react2(
        {1, 0, 0}, {1}, [](std::vector<double> y)
        { return 1.0 / (1.0 + y[1]); },
        test_system);
    myreact test_react3(
        {0, 1, 0}, {0}, [](std::vector<double> y)
        { return 1.0 / (1.0 + y[0]); },
        test_system);

    vector<Index> sigma1(test_system.mu()), sigma2(test_system.mu());
    const Index m1 = 2;
    const Index m2 = 1;
    multi_array<Index, 1> n1({m1});
    multi_array<Index, 1> k1({m1});
    multi_array<double, 1> liml1({m1});
    multi_array<Index, 1> n2({m2});
    multi_array<Index, 1> k2({m2});
    multi_array<double, 1> liml2({m2});
    n1(0) = 4;
    n1(1) = 3;
    k1(0) = 1;
    k1(1) = 1;
    liml1(0) = 0.0;
    liml1(1) = 0.0;
    n2(0) = 2;
    k2(0) = 1;
    liml2(0) = 0.0;
    grid_info grid(m1, m2, 1, n1, n2, k1, k2, liml1, liml2);
    CalculateShiftAmount(sigma1, sigma2, test_system, grid);

    SECTION("VecIndexToState")
    {
        Index comb_index = 0;
        Index stride = 1;
        multi_array<Index, 1> interval({5});
        multi_array<double, 2> lim({5, 2});
        vector<double> state_vec(5);
        vector<double> comparison_vec(5);
        for (Index i = 0; i < 5; i++)
        {
            interval(i) = 5 * (i + 1) + 1;
            comb_index += 5 * i * stride;
            stride *= interval(i);
            lim(i, 0) = 5.0 * (i + 1.0);
            lim(i, 1) = 10.0 * (i + 1.0);
            comparison_vec[i] = 5.0 + 10.0 * i;
        }
        CombIndexToState(state_vec, comb_index, interval, lim);
        REQUIRE(bool(state_vec == comparison_vec));
    }

    SECTION("VecIndexToCombIndex")
    {
        multi_array<Index, 1> vec_index({10});
        multi_array<Index, 1> interval({10});
        Index comb_index;
        Index comparison_index = 592088944020;
        for (Index i = 0; i < 10; i++)
        {
            vec_index(i) = i;
            interval(i) = 20 - i;
        }
        comb_index = VecIndexToCombIndex(vec_index, interval);
        REQUIRE(bool(comb_index == comparison_index));
    }

    SECTION("CombIndexToVecIndex1")
    {
        Index comb_index = 23084307895;
        multi_array<Index, 1> interval({10});
        multi_array<Index, 1> vec_index({10});
        multi_array<Index, 1> comparison_vec({10});
        for (Index i = 0; i < 10; i++)
        {
            interval(i) = 11;
            comparison_vec(i) = i;
        }
        CombIndexToVecIndex(vec_index, comb_index, interval);
        REQUIRE(bool(vec_index == comparison_vec));
    }

    SECTION("CombIndexToVecIndex2")
    {
        Index comb_index = 79;
        vector<Index> interval;
        vector<Index> vec_index(4);
        vector<Index> comparison_vec;
        interval = {4, 2, 3, 5};
        comparison_vec = {3, 1, 0, 3};
        CombIndexToVecIndex(vec_index, comb_index, interval);
        REQUIRE(bool(vec_index == comparison_vec));
    }

    SECTION("CalculateShiftAmount")
    {
        vector<Index> sigma1_comparison, sigma2_comparison;

        sigma1_comparison = {-5, -4, 1, 4};
        sigma2_comparison = {0, -1, 0, 0};
        REQUIRE(bool(sigma1 == sigma1_comparison));
        REQUIRE(bool(sigma2 == sigma2_comparison));
    }

    SECTION("ShiftMultiArrayRows1")
    {
        Index n_rows = grid.dx1;
        Index n_cols = 1;
        multi_array<double, 2> input_array({n_rows, n_cols}), output_array({n_rows, n_cols});
        multi_array<double, 2> comparison_array({n_rows, n_cols});

        set_zero(input_array);
        for (Index i = 0; i < n_rows; i++)
            input_array(i, 0) = (double) i + 1;

        set_zero(comparison_array);
        comparison_array(0, 0) = 6.0;
        comparison_array(1, 0) = 7.0;
        comparison_array(2, 0) = 8.0;
        comparison_array(3, 0) = 0.0; // !
        comparison_array(4, 0) = 10.0;
        comparison_array(5, 0) = 11.0;
        comparison_array(6, 0) = 12.0;

        ShiftMultiArrayRows<1>(output_array, input_array, sigma1[0], test_system.reactions[0]->nu, grid);

        REQUIRE(bool(output_array == comparison_array));
    }

    SECTION("ShiftMultiArrayRows2")
    {
        Index n_rows = grid.dx1;
        Index n_cols = 2;
        multi_array<double, 2> input_array({n_rows, n_cols}), output_array({n_rows, n_cols});
        multi_array<double, 2> comparison_array({n_rows, n_cols});

        set_zero(input_array);
        for (Index i = 0; i < n_rows; i++)
        {
            input_array(i, 0) = (double) i + 1;
            input_array(i, 1) = (double) i + 2;
        }

        set_zero(comparison_array);
        comparison_array(5, 0) = 1.0;
        comparison_array(6, 0) = 2.0;
        comparison_array(7, 0) = 3.0;
        comparison_array(8, 0) = 0.0; // !
        comparison_array(9, 0) = 5.0;
        comparison_array(10, 0) = 6.0;
        comparison_array(11, 0) = 7.0;

        comparison_array(5, 1) = 2.0;
        comparison_array(6, 1) = 3.0;
        comparison_array(7, 1) = 4.0;
        comparison_array(8, 1) = 0.0; // !
        comparison_array(9, 1) = 6.0;
        comparison_array(10, 1) = 7.0;
        comparison_array(11, 1) = 8.0;

        ShiftMultiArrayRows<1>(output_array, input_array, -sigma1[0], test_system.reactions[0]->minus_nu, grid);

        REQUIRE(bool(output_array == comparison_array));
    }

    SECTION("ShiftMultiArrayRows3")
    {
        Index n_rows = grid.dx2;
        Index n_cols = 1;
        multi_array<double, 2> input_array({n_rows, n_cols}), output_array({n_rows, n_cols});
        multi_array<double, 2> comparison_array({n_rows, n_cols});

        set_zero(input_array);
        for (Index i = 0; i < n_rows; i++)
        {
            input_array(i, 0) = (double) i + 1;
        }

        set_zero(comparison_array);
        comparison_array(0, 0) = 2.0;

        ShiftMultiArrayRows<2>(output_array, input_array, sigma2[1], test_system.reactions[1]->nu, grid);

        REQUIRE(bool(output_array == comparison_array));
    }
}