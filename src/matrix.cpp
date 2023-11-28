#include "matrix.hpp"

void Matrix::ShiftMultiArrayRows(multi_array<double, 2> &output_array, const multi_array<double, 2> &input_array, const grid_parms grid, const Index mu)
{
    assert(output_array.shape() == input_array.shape());

    Index n_rows = output_array.shape()[0];
    Index n_cols = output_array.shape()[1];

    // NOTE: Ensign stores matrices in column-major order

    vector<Index> vec_index(grid.d);

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
                if ((grid.shift[mu] < 0 && i - grid.shift[mu] < n_rows) || (grid.shift[mu] >= 0 && i - grid.shift[mu] >= 0))
                {
                    output_array(i, j) = input_array(i - grid.shift[mu], j);
                }
                else
                {
                    output_array(i, j) = 0.0;
                    IncrVecIndex(std::begin(grid.n), std::begin(vec_index), std::end(vec_index));
                    continue;
                }

                for (int k = 0; k < grid.d; k++)
                {
                    if (
                        ((grid.nu(mu, k) > 0) && (vec_index[k] - grid.nu(mu, k) < 0)) ||
                        ((grid.nu(mu, k) < 0) && (vec_index[k] - grid.nu(mu, k) >= grid.n[k])))
                    {
                        output_array(i, j) = 0.0;
                        break;
                    }
                }
                IncrVecIndex(std::begin(grid.n), std::begin(vec_index), std::end(vec_index));
            }
        }
    }
}