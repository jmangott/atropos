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
    void Matricize(const multi_array<T, d> &input, multi_array<T, 2> &output, Index dim)
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
            j = VecIndexToCombIndex(std::begin(vec_index_cols), std::end(vec_index_cols), std::begin(cols_shape));
            output(j, i) = el;
            IncrVecIndex(std::begin(shape), std::begin(vec_index), std::end(vec_index));
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
            j = VecIndexToCombIndex(std::begin(vec_index_cols), std::end(vec_index_cols), std::begin(cols_shape));
            el = input(j, i);
            IncrVecIndex(std::begin(shape), std::begin(vec_index), std::end(vec_index));
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
}

#endif