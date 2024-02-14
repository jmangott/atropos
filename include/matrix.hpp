#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <algorithm>
#include <cassert>

#include <generic/storage.hpp>
#include <lr/lr.hpp>

#include "index_functions.hpp"
#include "timer_class.hpp"

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
        set_zero(output_array);

        Index shift = inv * grid.shift[mu];
        Index n_rows = output_array.shape()[0];
        Index n_cols = output_array.shape()[1];

        Index min_i = std::max((Index) 0, shift);
        Index max_i = std::min(n_rows, n_rows + shift);

        // NOTE: Ensign stores matrices in column-major order
        // TODO: parallelize loops
        for (Index j = 0; j < n_cols; ++j)
        {
            for (Index i = min_i; i < max_i; ++i)
            {
                output_array(i, j) = input_array(i - shift, j);
            }
        }

        std::vector<Index> vec_index(grid.d);
        IndexFunction::CombIndexToVecIndex(min_i, std::begin(grid.n), std::begin(vec_index), std::end(vec_index));
        for (Index i = min_i; i < max_i; ++i)
        {
            for (Index k = 0; k < grid.d; ++k)
            {
                if ((vec_index[k] - inv * grid.nu(mu, k) < 0) || (vec_index[k] - inv * grid.nu(mu, k) >= grid.n[k]))
                {
                    for (Index j = 0; j < n_cols; ++j)
                    {
                        output_array(i, j) = 0.0;
                    }
                    break;
                }
            }
            IndexFunction::IncrVecIndex(std::begin(grid.n), std::begin(vec_index), std::end(vec_index));
        }
    }
}

#endif