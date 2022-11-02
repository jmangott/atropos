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


void CalculateShiftAmount(vector<Index> &sigma1, vector<Index> &sigma2, mysys reaction_system, multi_array<Index, 1> n_xx1, multi_array<Index, 1> n_xx2, multi_array<Index, 1> k_xx1, multi_array<Index, 1> k_xx2)
{
    Index stride1, stride2;
    Index sigma1_sum, sigma2_sum;
    Index m1 = n_xx1.shape()[0];
    Index m2 = n_xx2.shape()[0];

    // NOTE: when the partition requires a permutation of the original order of species,
    // then also the nu vectors and similar quantities have to be permuted

    for (auto &it : reaction_system.reactions)
    {
        stride1 = 1;
        stride2 = 1;
        sigma1_sum = 0;
        sigma2_sum = 0;
        for (Index i = 0; i < m1; i++)
        {
            sigma1_sum += it->nu[i] * stride1 * k_xx1(i);
            stride1 *= n_xx1(i);
        }
        for (Index i = 0; i < m2; i++)
        {
            sigma2_sum += it->nu[i + m1] * stride2 * k_xx2(i);
            stride2 *= n_xx2(i);
        }
        sigma1.push_back(sigma1_sum);
        sigma2.push_back(sigma2_sum);
    }
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
            if ((shift < 0) && (i < -shift))
            {
                output_array(i, j) = input_array(0, j);
            }
            else if ((shift > 0) && (i >= (n_rows - shift)))
            {
                output_array(i, j) = input_array(n_rows - 1, j);
            }
            else
            {
                output_array(i, j) = input_array(i + shift, j);
            }
        }
    }
}