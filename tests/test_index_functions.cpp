#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"

TEST_CASE("IncrVecIndex", "[IncrVecIndex]")
{
    std::array<Index, 5> input = {0, 1, 1, 1, 3};
    std::array<Index, 5> interval = {1, 2, 3, 4, 5};

    std::array<Index, 5> output1_ref = {0, 0, 2, 1, 3};
    std::array<Index, 5> output2_ref = {0, 0, 0, 2, 3};
    std::array<Index, 5> output3_ref = {0, 0, 0, 0, 4};

    IndexFunction::IncrVecIndex(std::begin(interval), std::begin(input), std::end(input));
    REQUIRE(bool(input == output1_ref));

    IndexFunction::IncrVecIndex(std::begin(interval), std::begin(input), std::end(input));
    IndexFunction::IncrVecIndex(std::begin(interval), std::begin(input), std::end(input));
    REQUIRE(bool(input == output2_ref));

    input = {0, 1, 2, 3, 3};
    IndexFunction::IncrVecIndex(std::begin(interval), std::begin(input), std::end(input));
    REQUIRE(bool(input == output3_ref));
}

TEST_CASE("VecIndexToCombIndex", "[VecIndexToCombIndex]")
{
    std::vector<Index> vec_index(10);
    std::vector<Index> interval(10);
    Index comb_index;
    Index comparison_index = 592088944020;
    for (Index i = 0; i < 10; ++i)
    {
        vec_index[i] = i;
        interval[i] = 20 - i;
    }
    comb_index = IndexFunction::VecIndexToCombIndex(std::begin(vec_index), std::end(vec_index), std::begin(interval));
    REQUIRE(bool(comb_index == comparison_index));
}

TEST_CASE("CombIndexToVecIndex", "[CombIndexToVecIndex]")
{
    Index comb_index = 23084307895;
    std::vector<Index> interval(10);
    std::vector<Index> vec_index(10);
    std::vector<Index> comparison_vec(10);
    for (Index i = 0; i < 10; ++i)
    {
        interval[i] = 11;
        comparison_vec[i] = i;
    }
    IndexFunction::CombIndexToVecIndex(comb_index, std::begin(interval), std::begin(vec_index), std::end(vec_index));
    REQUIRE(bool(vec_index == comparison_vec));

    comb_index = 79;
    vec_index.resize(4);
    comparison_vec.resize(4);
    interval = {4, 2, 3, 5};
    comparison_vec = {3, 1, 0, 3};
    IndexFunction::CombIndexToVecIndex(comb_index, std::begin(interval), std::begin(vec_index), std::end(vec_index));
    REQUIRE(bool(vec_index == comparison_vec));
}

TEST_CASE("VecIndexToDepCombIndex", "[VecIndexToDepCombIndex]")
{
    std::vector<Index> vec_index = {6, 3, 2, 4, 11};
    std::vector<Index> n_dep = {6, 5, 13};
    std::vector<Index> idx_dep = {1, 3, 4};
    Index comb_index = IndexFunction::VecIndexToDepCombIndex(std::begin(vec_index), std::begin(n_dep), std::begin(idx_dep), std::end(idx_dep));

    Index comparison_comb_index = 357;
    REQUIRE(bool(comb_index == comparison_comb_index));
}

#ifdef __OPENMP__
TEST_CASE("SetVecIndex", "[SetVecIndex]")
{
    omp_set_dynamic(0);
    omp_set_num_threads(4);
    Index dx = 399;
    vector<Index> interval = {7, 19, 3};

#pragma omp parallel
    {
        std::vector<Index> vec_index(3);
        std::vector<Index> comparison_vec_index(3);
        Index chunk_size = IndexFunction::SetVecIndex(std::begin(vec_index), std::end(vec_index), std::begin(interval), dx);

        REQUIRE((chunk_size == 100));
        switch (omp_get_num_threads())
        {
            case 0:
                comparison_vec_index = {0, 0, 0};
                REQUIRE(bool(vec_index == comparison_vec_index));
            case 1:
                comparison_vec_index = {2, 14, 0};
                REQUIRE(bool(vec_index == comparison_vec_index));
            case 2:
                comparison_vec_index = {4, 9, 1};
                REQUIRE(bool(vec_index == comparison_vec_index));
            case 3:
                comparison_vec_index = {6, 4, 2};
                REQUIRE(bool(vec_index == comparison_vec_index));
        }
    }
}
#endif

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
        multi_array<double, 1> liml({5});
        multi_array<Index, 1> binsize({5});
        vector<double> state_vec(5);
        vector<double> comparison_vec(5);
        for (Index i = 0; i < 5; i++)
        {
            interval(i) = 5 * (i + 1) + 1;
            binsize(i) = 1;
            comb_index += 5 * i * stride;
            stride *= interval(i);
            liml(i) = 5.0 * (i + 1.0);
            comparison_vec[i] = 5.0 + 10.0 * i;
        }
        CombIndexToState(state_vec, comb_index, interval, liml, binsize);
        REQUIRE(bool(state_vec == comparison_vec));
    }

    SECTION("VecIndexToCombIndex")
    {
        vector<Index> vec_index(10);
        multi_array<Index, 1> interval({10});
        Index comb_index;
        Index comparison_index = 592088944020;
        for (Index i = 0; i < 10; i++)
        {
            vec_index[i] = i;
            interval(i) = 20 - i;
        }
        comb_index = VecIndexToCombIndex(vec_index, interval);
        REQUIRE(bool(comb_index == comparison_index));
    }

    SECTION("CombIndexToVecIndex1")
    {
        Index comb_index = 23084307895;
        multi_array<Index, 1> interval({10});
        vector<Index> vec_index(10);
        vector<Index> comparison_vec(10);
        for (Index i = 0; i < 10; i++)
        {
            interval(i) = 11;
            comparison_vec[i] = i;
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