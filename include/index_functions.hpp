#ifndef INDEX_FUNCTIONS_HPP
#define INDEX_FUNCTIONS_HPP

#include <iterator>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <generic/storage.hpp>

#ifdef __OPENMP__
#include <omp.h>
#endif

#include "grid_class.hpp"
#include "reaction_class.hpp"

// TODO: write tests

namespace IndexFunction
{
    template<class InputIt, class InputItInt>
    Index VecIndexToCombIndex(InputIt first, InputIt last, InputItInt first_int)
    {
        Index comb_index = 0;
        Index stride = 1;
        for (; first != last; ++first, ++first_int)
        {
            comb_index += *first * stride;
            stride *= *first_int;
        }
        return comb_index;
    }

    template <class InputIt, class OutputIt>
    void CombIndexToVecIndex(InputIt first, OutputIt d_first, OutputIt d_last)
    {
        assert(d_first != d_last);
        Index comb_index;
        for (; d_first != std::next(d_last, -1); ++first, ++d_first)
        {
            *d_first = comb_index % *first;
            comb_index = Index(comb_index / *first);
        }
        *(std::next(d_last, -1)) = comb_index;
    }

    template <class InputIt, class InputItInt, class InputItDep>
    Index VecIndexToDepCombIndex(InputIt first, InputItInt first_int, InputItDep first_dep, InputIt last_dep)
    {
        Index comb_index = 0;
        Index stride = 1;
        for (; first_dep != last_dep; ++first_int, ++first_dep)
        {
            comb_index += *std::next(first, *first_dep) * stride;
            stride *= *first_int;
        }
        return comb_index;
    }

    template <class InputIt, class OutputIt>
    void IncrVecIndex(InputIt first, OutputIt d_first, OutputIt d_last)
    {
        assert(d_first != d_last);
        for (; d_first != std::next(d_last, -1); ++first, ++d_first)
        {
            ++(*d_first);
            if (*d_first < *first)
                return;
            *d_first = typename std::iterator_traits<OutputIt>::value_type(0);
        }
        ++(*(std::next(d_last, -1)));
    }

    #ifdef __OPENMP__
    template <InputIt, InputItInt>
    Index SetVecIndex(InputIt first, InputIt last, InputItInt first_int, Index dx)
    {
        Index chunk_size, start_index;
        int num_threads = omp_get_num_threads();
        int thread_num = omp_get_thread_num();
        chunk_size = (Index)std::ceil((double)dx / num_threads);
        start_index = thread_num * chunk_size;
        CombIndexToVecIndex(first_int, first, last);
        return chunk_size;
    }
    #endif
}

//-----------------------------------------------------------------------------

Index VecIndexToCombIndex(std::vector<Index> vec_index, multi_array<Index, 1> interval);

inline void CombIndexToVecIndex(std::vector<Index> &vec_index, Index comb_index, const std::vector<Index> &interval)
{
    Index dim = interval.size();
    for (Index i = 0; i < (dim - 1); i++)
    {
        vec_index[i] = comb_index % interval[i];
        comb_index = Index (comb_index / interval[i]);
    }
    if (dim > 0) vec_index[dim - 1] = comb_index;
}

inline void CombIndexToVecIndex(std::vector<Index> &vec_index, Index comb_index, const multi_array<Index, 1> &interval)
{
    Index dim = interval.shape()[0];
    for (Index i = 0; i < (dim - 1); i++)
    {
        vec_index[i] = comb_index % interval(i);
        comb_index = Index (comb_index / interval(i));
    }
    if (dim > 0) vec_index[dim - 1] = comb_index;
}

inline void CombIndexToState(std::vector<double> &state, Index comb_index, const multi_array<Index, 1> &interval, const multi_array<double, 1> &liml, const multi_array<Index, 1> &binsize)
{
    Index dim = interval.shape()[0];
    for (Index i = 0; i < (dim - 1); i++)
    {
        state[i] = liml(i) + binsize(i) * (comb_index % interval(i));
        comb_index = Index (comb_index / interval(i));
    }
    if (dim > 0) state[dim - 1] = liml(dim - 1) + binsize(dim - 1) * comb_index;
}


// Convert a combined index for the participating species in reaction mu to a general combined index, where the remaining species have population number 0
Index DepCombIndexToCombIndex(Index comb_index_dep, std::vector<Index> n_dep, multi_array<Index, 1> n, std::vector<Index> dep_vec);


inline void IncrVecIndex(std::vector<Index> &vec_index, const multi_array<Index, 1> &interval, Index dim)
{
    for (Index k = 0; k < dim - 1; k++)
    {
        vec_index[k]++;
        if (vec_index[k] < interval(k))
            return;
        vec_index[k] = 0;
    }
    if (dim > 0) vec_index[dim - 1]++;
}


#ifdef __OPENMP__
inline Index SetVecIndex(std::vector<Index> &vec_index, const multi_array<Index, 1> &interval, Index dx)
{
    Index chunk_size, start_index;
    int num_threads = omp_get_num_threads();
    int thread_num = omp_get_thread_num();
    chunk_size = (Index) std::ceil((double) dx / num_threads);
    start_index = thread_num * chunk_size;
    CombIndexToVecIndex(vec_index, start_index, interval);
    return chunk_size;
}
#endif

inline Index VecIndexToDepCombIndex(std::vector<Index> &vec_index, const std::vector<Index> &n_dep, const std::vector<Index> &dep_vec)
{
    Index comb_index = 0;
    Index stride = 1;
    for (vector<Index>::size_type i = 0; i < dep_vec.size(); i++)
    {
        comb_index += vec_index[dep_vec[i]] * stride;
        stride *= n_dep[i];
    }
    return comb_index;
}

// Calculate for all reactions the number of indices by which arrays have to be shifted in order to calculate coefficients
void CalculateShiftAmount(std::vector<Index> &sigma1, std::vector<Index> &sigma2, mysys reaction_system, grid_info grid);

// Calculate `output_array`, where rows of `input_array` are shifted by `shift`
// NOTE: for positive values of `shift` rows are shifted to larger row indices
template <int id>
void ShiftMultiArrayRows(multi_array<double, 2> &output_array, const multi_array<double, 2> &input_array, Index shift, vector<int> nu, grid_info grid)
{
    if ((output_array.shape()[0] != input_array.shape()[0]) ||
        (output_array.shape()[1] != input_array.shape()[1]))
    {
        std::cerr << "ERROR: Dimensions of output_array and input_array must be the same!" << endl;
        std::abort();
    }

    Index n_rows = output_array.shape()[0];
    Index n_cols = output_array.shape()[1];

    grid_info *grid_alt;
    Index inc;

    if constexpr (id == 1)
    {
        grid_alt = new grid_info(grid.m1, grid.m2, grid.r, grid.n1, grid.n2, grid.binsize1, grid.binsize2, grid.liml1, grid.liml2);
        inc = 0;
    }
    else if constexpr (id == 2)
    {
        grid_alt = new grid_info(grid.m2, grid.m1, grid.r, grid.n2, grid.n1, grid.binsize2, grid.binsize1, grid.liml2, grid.liml1);
        inc = grid.m1;
    }
    else
    {
        std::cerr << "ERROR: `id` must be 1 (shift in partition 1) or 2 (shift in partition 2)!" << endl;
        std::abort();
    }

    Index k_inc;

    // NOTE: Ensign stores matrices in column-major order

    vector<Index> vec_index(grid_alt->m1);

    for (Index j = 0; j < n_cols; j++)
    {
        std::fill(vec_index.begin(), vec_index.end(), 0);
        
#ifdef __OPENMP__
#pragma omp parallel firstprivate(vec_index) private(k_inc)
#endif
        {
#ifdef __OPENMP__
            Index chunk_size = SetVecIndex(vec_index, grid_alt->n1, grid_alt->dx1);
#endif

#ifdef __OPENMP__
#pragma omp for schedule(static, chunk_size)
#endif
            for (Index i = 0; i < n_rows; i++)
            {
                if ((shift < 0 && i - shift < n_rows) || (shift >= 0 && i - shift >= 0))
                {
                    output_array(i, j) = input_array(i - shift, j);
                }
                else
                {
                    output_array(i, j) = 0.0;
                    IncrVecIndex(vec_index, grid_alt->n1, grid_alt->m1);
                    continue;
                }

                for (int k = 0; k < grid_alt->m1; k++)
                {
                    k_inc = k + inc;
                    if (((nu[k_inc] > 0) &&
                        (vec_index[k] - nu[k_inc] < 0)) ||
                        ((nu[k_inc] < 0) &&
                        (vec_index[k] - nu[k_inc] >= grid_alt->n1(k))))
                    {
                        output_array(i, j) = 0.0;
                        break;
                    }
                }

                IncrVecIndex(vec_index, grid_alt->n1, grid_alt->m1);
            }
        }
    }
    delete grid_alt;
}

#endif