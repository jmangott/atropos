#ifndef INDEX_FUNCTIONS_HPP
#define INDEX_FUNCTIONS_HPP

#include <iterator>
#include <vector>

#include <generic/storage.hpp>

#ifdef __OPENMP__
#include <omp.h>
#endif

#include "grid_class.hpp"

template <class InputIt, class InputItInt, class InputItDep>
Index VecIndexToDepCombIndex(InputIt first, InputItInt first_int, InputItDep first_dep,
                             InputItDep last_dep)
{
    Index comb_index = 0;
    Index stride = 1;
    for (; first_dep != last_dep; ++first_int, ++first_dep) {
        comb_index += *std::next(first, *first_dep) * stride;
        stride *= *first_int;
    }
    return comb_index;
}

#ifdef __OPENMP__
template <class InputIt, class InputItInt>
Index SetVecIndex(InputIt first, InputIt last, InputItInt first_int, const Index dx)
{
    Index chunk_size, start_index;
    int num_threads = omp_get_num_threads();
    int thread_num = omp_get_thread_num();
    chunk_size = (Index)std::ceil((double)dx / num_threads);
    start_index = thread_num * chunk_size;
    Ensign::IndexFunction::comb_index_to_vec_index(start_index, first_int, first, last);
    return chunk_size;
}
#endif

#endif
