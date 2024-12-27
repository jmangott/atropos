#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "index_functions.hpp"

TEST_CASE("VecIndexToDepCombIndex", "[VecIndexToDepCombIndex]")
{
    std::vector<Index> vec_index = {6, 3, 2, 4, 11};
    std::vector<Index> n_dep = {6, 5, 13};
    std::vector<Index> idx_dep = {1, 3, 4};
    Index comb_index = VecIndexToDepCombIndex(std::begin(vec_index), std::begin(n_dep),
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
        Index chunk_size = SetVecIndex(std::begin(vec_index), std::end(vec_index),
                                       std::begin(interval), dx);

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
