#include "grid_class.hpp"
#include "index_functions.hpp"

using std::vector;

multi_array<double, 1> VecIndexToState(multi_array<Index, 1> vec_index, multi_array<Index, 1> interval, multi_array<double, 1> limit)
{
    Index dim = vec_index.shape()[0];
    multi_array<double, 1> state_vec({dim});
    for (Index i = 0; i < dim; i++)
    {
        state_vec(i) = vec_index(i) * limit(i) / (interval(i) - 1.0);
    }
    return state_vec;
}


Index VecIndexToCombIndex(multi_array<Index, 1> vec_index, multi_array<Index, 1> interval)
{
    Index comb_index = 0;
    Index stride = 1;
    for (Index i = 0; i < interval.shape()[0]; i++)
    {
        comb_index += vec_index(i) * stride;
        stride *= interval(i);
    }
    return comb_index;
}


multi_array<Index, 1> CombIndexToVecIndex(Index comb_index, multi_array<Index, 1> interval)
{
    Index dim = interval.shape()[0];
    multi_array<Index, 1> vec_index({dim});
    for (Index i = 0; i < dim; i++)
    {
        if (i == (dim - 1))
            vec_index(i) = comb_index;
        else
        {
            vec_index(i) = comb_index % interval(i);
            comb_index = int(comb_index / interval(i));
        }
    }
    return vec_index;
}


vector<Index> CombIndexToVecIndex(Index comb_index, vector<Index> interval)
{
    Index dim = interval.size();
    vector<Index> vec_index(dim, 0);
    for (Index i = 0; i < dim; i++)
    {
        if (i == (dim - 1))
            vec_index[i] = comb_index;
        else
        {
            vec_index[i] = comb_index % interval[i];
            comb_index = int(comb_index / interval[i]);
        }
    }
    return vec_index;
}


void ShiftMultiArrayRows(multi_array<double, 2> &output_array, multi_array<double, 2> &input_array, int shift)
{
    if ((output_array.shape()[0] != input_array.shape()[0]) ||
        (output_array.shape()[1] != input_array.shape()[1]))
    {
        std::cerr << "ERROR: Dimensions of output_array and input_array must be the same!";
        std::abort();
    }

    Index n_rows = output_array.shape()[0];
    Index n_cols = output_array.shape()[1];

    // NOTE: Ensign stores matrices in column-major order
    for (Index j = 0; j < n_cols; j++)
    {
        for (Index i = 0; i < n_rows; i++)
        {
            if ((shift < 0) && (i - shift >= n_rows))
            {
                output_array(i, j) = input_array(n_rows - 1, j);
            }
            else if ((shift > 0) && (i - shift < 0))
            {
                output_array(i, j) = input_array(0, j);
            }
            else
            {
                output_array(i, j) = input_array(i - shift, j);
            }
        }
    }
}