#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"

TEST_CASE("index_functions", "[index_functions]")
{
    SECTION("VecIndexToState")
    {
        multi_array<Index, 1> vec_index({5});
        multi_array<Index, 1> interval({5});
        multi_array<double, 1> limit({5});
        multi_array<double, 1> state_vec({5});
        multi_array<double, 1> comparison_vec({5});
        for (Index i = 0; i < 5; i++)
        {
            vec_index(i) = 5 * i;
            interval(i) = 5 * (i + 1) + 1;
            limit(i) = 10.0 * (i + 1);
            comparison_vec(i) = 10.0 * i;
        }
        state_vec = VecIndexToState(vec_index, interval, limit);
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
        vec_index = CombIndexToVecIndex(comb_index, interval);
        REQUIRE(bool(vec_index == comparison_vec));
    }

    SECTION("CombIndexToVecIndex2")
    {
        Index comb_index = 79;
        vector<Index> interval;
        vector<Index> vec_index;
        vector<Index> comparison_vec;
        interval = {4, 2, 3, 5};
        comparison_vec = {3, 1, 0, 3};
        vec_index = CombIndexToVecIndex(comb_index, interval);
        REQUIRE(bool(vec_index == comparison_vec));
    }

    SECTION("CalculateShiftAmount")
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
        vector<Index> sigma1, sigma2;
        vector<Index> sigma1_comparison, sigma2_comparison;
        
        const Index m1 = 2;
        const Index m2 = 1;
        multi_array<Index, 1> n1({m1});
        multi_array<Index, 1> k1({m1});
        multi_array<Index, 1> n2({m2});
        multi_array<Index, 1> k2({m2});
        n1(0) = 4;
        n1(1) = 3;
        k1(0) = 1;
        k1(1) = 2;
        n2(0) = 2;
        k2(0) = 1;
        grid_info<m1, m2> grid(2, n1, n2, k1, k2);

        sigma1_comparison = {-9, -8, 1, 8};
        sigma2_comparison = {0, -1, 0, 0};
        CalculateShiftAmount<m1, m2>(sigma1, sigma2, test_system, grid);
        REQUIRE(bool(sigma1 == sigma1_comparison));
        REQUIRE(bool(sigma2 == sigma2_comparison));
    }

    SECTION("ShiftMultiArrayRows1")
    {
        Index n_rows = 12;
        Index n_cols = 10;
        multi_array<double, 2> input_array({n_rows, n_cols}), output_array({n_rows, n_cols});
        multi_array<double, 2> comparison_array({n_rows, n_cols});

        set_zero(input_array);
        for (Index i = 0; i < n_cols; i++)
            input_array(i, i) = 1.0;

        set_zero(comparison_array);
        for (Index i = 0; i < 5; i++)
            comparison_array(i, 0) = 1.0;
        for (Index i = 0; i < n_rows - 5; i++)
            comparison_array(i + 5, i) = 1.0;

        ShiftMultiArrayRows(output_array, input_array, 5);
        REQUIRE(bool(output_array == comparison_array));
    }

    SECTION("ShiftMultiArrayRows2")
    {
        Index n_rows = 12;
        Index n_cols = 10;
        multi_array<double, 2> input_array({n_rows, n_cols}), output_array({n_rows, n_cols});
        multi_array<double, 2> comparison_array({n_rows, n_cols});

        set_zero(input_array);
        for (Index i = 0; i < n_cols; i++)
            input_array(i, i) = 1.0;

        set_zero(comparison_array);
        for (Index i = 5; i < n_cols; i++)
            comparison_array(i - 5, i) = 1.0;

        ShiftMultiArrayRows(output_array, input_array, -5);
        REQUIRE(bool(output_array == comparison_array));
    }
}