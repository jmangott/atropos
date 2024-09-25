#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <generic/index_functions.hpp>

#include "index_functions.hpp"

TEST_CASE("IncrVecIndex", "[IncrVecIndex]")
{
    std::array<Index, 5> input = {0, 1, 1, 1, 3};
    std::array<Index, 5> interval = {1, 2, 3, 4, 5};

    std::array<Index, 5> output1_ref = {0, 0, 2, 1, 3};
    std::array<Index, 5> output2_ref = {0, 0, 0, 2, 3};
    std::array<Index, 5> output3_ref = {0, 0, 0, 0, 4};

    IndexFunction::IncrVecIndex(std::begin(interval), std::begin(input),
                                std::end(input));
    REQUIRE(bool(input == output1_ref));

    IndexFunction::IncrVecIndex(std::begin(interval), std::begin(input),
                                std::end(input));
    IndexFunction::IncrVecIndex(std::begin(interval), std::begin(input),
                                std::end(input));
    REQUIRE(bool(input == output2_ref));

    input = {0, 1, 2, 3, 3};
    IndexFunction::IncrVecIndex(std::begin(interval), std::begin(input),
                                std::end(input));
    REQUIRE(bool(input == output3_ref));
}

TEST_CASE("VecIndexToCombIndex", "[VecIndexToCombIndex]")
{
    std::vector<Index> vec_index(10);
    std::vector<Index> interval(10);
    Index comb_index;
    Index comparison_index = 592088944020;
    for (Index i = 0; i < 10; ++i) {
        vec_index[i] = i;
        interval[i] = 20 - i;
    }
    comb_index = IndexFunction::VecIndexToCombIndex(
        std::begin(vec_index), std::end(vec_index), std::begin(interval));
    REQUIRE(bool(comb_index == comparison_index));
}

TEST_CASE("CombIndexToVecIndex", "[CombIndexToVecIndex]")
{
    Index comb_index = 23084307895;
    std::vector<Index> interval(10);
    std::vector<Index> vec_index(10);
    std::vector<Index> comparison_vec(10);
    for (Index i = 0; i < 10; ++i) {
        interval[i] = 11;
        comparison_vec[i] = i;
    }
    IndexFunction::CombIndexToVecIndex(comb_index, std::begin(interval),
                                       std::begin(vec_index), std::end(vec_index));
    REQUIRE(bool(vec_index == comparison_vec));

    comb_index = 79;
    vec_index.resize(4);
    comparison_vec.resize(4);
    interval = {4, 2, 3, 5};
    comparison_vec = {3, 1, 0, 3};
    IndexFunction::CombIndexToVecIndex(comb_index, std::begin(interval),
                                       std::begin(vec_index), std::end(vec_index));
    REQUIRE(bool(vec_index == comparison_vec));
}

TEST_CASE("VecIndexToDepCombIndex", "[VecIndexToDepCombIndex]")
{
    std::vector<Index> vec_index = {6, 3, 2, 4, 11};
    std::vector<Index> n_dep = {6, 5, 13};
    std::vector<Index> idx_dep = {1, 3, 4};
    Index comb_index =
        IndexFunction::VecIndexToDepCombIndex(std::begin(vec_index), std::begin(n_dep),
                                              std::begin(idx_dep), std::end(idx_dep));

    Index comparison_comb_index = 357;
    REQUIRE(bool(comb_index == comparison_comb_index));
}

#ifdef __OPENMP__
TEST_CASE("SetVecIndex", "[SetVecIndex]")
{
    omp_set_dynamic(0);
    omp_set_num_threads(4);
    Index dx = 399;
    std::vector<Index> interval = {7, 19, 3};

#pragma omp parallel
    {
        std::vector<Index> vec_index(3);
        std::vector<Index> comparison_vec_index(3);
        Index chunk_size = IndexFunction::SetVecIndex(
            std::begin(vec_index), std::end(vec_index), std::begin(interval), dx);

        REQUIRE((chunk_size == 100));
        switch (omp_get_num_threads()) {
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
