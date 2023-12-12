#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cassert>

#include <generic/storage.hpp>
#include <lr/lr.hpp>

#include "index_functions.hpp"

// TODO: Write tests for these functions

namespace Matrix
{
    template <class InputIt, class OutputIt>
    void RemoveElement(InputIt first, InputIt last, OutputIt d_first, const Index idx)
    {
        std::copy(first, first + idx, d_first);
        std::copy(first + idx + 1, last, d_first + idx);
    }

    template <class T, size_t d>
    void Matricize(const multi_array<T, d> &input, multi_array<T, 2> &output, const Index dim)
    {
        std::array<Index, d> shape{input.shape()};
        std::array<Index, d - 1> cols_shape, vec_index_cols;
        RemoveElement(std::begin(shape), std::end(shape), std::begin(cols_shape), dim);
        std::vector<Index> vec_index(d, 0);
        Index i, j;

        assert(shape[dim] == output.shape()[1] && prod(cols_shape) == output.shape()[0]);

        for (auto const &el : input)
        {
            i = vec_index[dim];
            RemoveElement(std::begin(vec_index), std::end(vec_index), std::begin(vec_index_cols), dim);
            j = IndexFunction::VecIndexToCombIndex(std::begin(vec_index_cols), std::end(vec_index_cols), std::begin(cols_shape));
            output(j, i) = el;
            IndexFunction::IncrVecIndex(std::begin(shape), std::begin(vec_index), std::end(vec_index));
        }
    }

    template <class T, size_t d>
    void Tensorize(const multi_array<T, 2> &input, multi_array<T, d> &output, const Index dim) 
    {
        std::array<Index, d> shape{output.shape()};
        std::array<Index, d - 1> cols_shape, vec_index_cols;
        RemoveElement(std::begin(shape), std::end(shape), std::begin(cols_shape), dim);
        std::vector<Index> vec_index(d, 0);
        Index i, j;

        assert(shape[dim] == input.shape()[1] && prod(cols_shape) == input.shape()[0]);

        for (auto &el : output)
        {
            i = vec_index[dim];
            RemoveElement(std::begin(vec_index), std::end(vec_index), std::begin(vec_index_cols), dim);
            j = IndexFunction::VecIndexToCombIndex(std::begin(vec_index_cols), std::end(vec_index_cols), std::begin(cols_shape));
            el = input(j, i);
            IndexFunction::IncrVecIndex(std::begin(shape), std::begin(vec_index), std::end(vec_index));
        }
    }

    template <class T>
    multi_array<T, 2> Orthogonalize(multi_array<T, 2> &input, const Index n_basisfunctions, std::function<T(T *, T *)> inner_product, const blas_ops &blas)
    {
        Index rows = input.shape()[0];
        Index cols = input.shape()[1];

        assert(n_basisfunctions <= cols);

        std::default_random_engine generator(1234);
        std::normal_distribution<double> distribution(0.0, 1.0);

        for (Index k = n_basisfunctions; k < cols; ++k)
        {
#ifdef __OPENMP__
#pragma omp parallel for
#endif
            for (Index i = 0; i < rows; ++i)
            {
                input(i, k) = distribution(generator);
            }
        }

        multi_array<T, 2> R({cols, cols});

        gram_schmidt gs(&blas);
        gs(input, R, inner_product);

        for (Index j = n_basisfunctions; j < cols; ++j)
        {
            for (Index i = 0; i < cols; ++i)
            {
                R(i, j) = T(0.0);
            }
        }
        return R;
    }

    template <Index inv>
    void ShiftRows(multi_array<double, 2> &output_array, const multi_array<double, 2> &input_array, const grid_parms grid, const Index mu)
    {
        assert(output_array.shape() == input_array.shape());

        Index shift = inv * grid.shift[mu];
        Index n_rows = output_array.shape()[0];
        Index n_cols = output_array.shape()[1];

        // NOTE: Ensign stores matrices in column-major order

        vector<Index> vec_index(grid.d);

        for (Index j = 0; j < n_cols; ++j)
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
                for (Index i = 0; i < n_rows; ++i)
                {
                    if ((shift < 0 && i - shift < n_rows) || (shift >= 0 && i - shift >= 0))
                    {
                        output_array(i, j) = input_array(i - shift, j);
                    }
                    else
                    {
                        output_array(i, j) = 0.0;
                        IndexFunction::IncrVecIndex(std::begin(grid.n), std::begin(vec_index), std::end(vec_index));
                        continue;
                    }

                    for (int k = 0; k < grid.d; ++k)
                    {
                        if (
                            ((inv * grid.nu(mu, k) > 0) && (vec_index[k] - inv * grid.nu(mu, k) < 0)) ||
                            ((inv * grid.nu(mu, k) < 0) && (vec_index[k] - inv * grid.nu(mu, k) >= grid.n[k])))
                        {
                            output_array(i, j) = 0.0;
                            break;
                        }
                    }
                    IndexFunction::IncrVecIndex(std::begin(grid.n), std::begin(vec_index), std::end(vec_index));
                }
            }
        }
    }
}

#endif