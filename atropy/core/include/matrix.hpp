#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <algorithm>
#include <cassert>

#include <generic/index.hpp>
#include <generic/storage.hpp>
#include <generic/timer.hpp>
#include <lr/lr.hpp>

#include "index_functions.hpp"

namespace Matrix {
template <Index inv>
void ShiftRows(Ensign::multi_array<double, 2>& output_array,
               const Ensign::multi_array<double, 2>& input_array, const grid_parms grid,
               const Index mu)
{
    assert(output_array.shape() == input_array.shape());
    Ensign::Matrix::set_zero(output_array);

    Index shift = inv * grid.shift[mu];
    Index n_rows = output_array.shape()[0];
    Index n_cols = output_array.shape()[1];

    Index min_i = std::max((Index)0, shift);
    Index max_i = std::min(n_rows, n_rows + shift);

    // NOTE: Ensign stores matrices in column-major order
#ifdef __OPENMP__
#pragma omp parallel for
#endif
    for (Index j = 0; j < n_cols; ++j) {
        for (Index i = min_i; i < max_i; ++i) {
            output_array(i, j) = input_array(i - shift, j);
        }
    }

    std::vector<Index> vec_index(grid.d);
    Ensign::IndexFunction::comb_index_to_vec_index(
        min_i, std::begin(grid.n), std::begin(vec_index), std::end(vec_index));
    for (Index i = min_i; i < max_i; ++i) {
        for (Index k = 0; k < grid.d; ++k) {
            if ((vec_index[k] - inv * grid.nu(mu, k) < 0) ||
                (vec_index[k] - inv * grid.nu(mu, k) >= grid.n[k])) {
// TODO: improve parallelization
#ifdef __OPENMP__
#pragma omp parallel for
#endif
                for (Index j = 0; j < n_cols; ++j) {
                    output_array(i, j) = 0.0;
                }
                break;
            }
        }
        Ensign::IndexFunction::incr_vec_index(std::begin(grid.n), std::begin(vec_index),
                                              std::end(vec_index));
    }
}
} // namespace Matrix

#endif
