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

Index VecIndexToCombIndex(vector<Index> vec_index, multi_array<Index, 1> interval)
{
    Index comb_index = 0;
    Index stride = 1;
    for (Index i = 0; i < interval.shape()[0]; i++)
    {
        comb_index += vec_index[i] * stride;
        stride *= interval(i);
    }
    return comb_index;
}


Index VecIndexToCombIndex(vector<Index> vec_index, vector<Index> interval)
{
    Index comb_index = 0;
    Index stride = 1;
    for (vector<Index>::size_type i = 0; i < interval.size(); i++)
    {
        comb_index += vec_index[i] * stride;
        stride *= interval[i];
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


Index DepCombIndexToCombIndex(Index comb_index_dep, vector<Index> n_dep, multi_array<Index, 1> n, vector<Index> dep_vec)
{
    // cout << "deptest0";
    vector<Index> vec_index(n.shape()[0], 0);
    vector<Index> vec_index_dep;
    Index comb_index;
    // cout << "1";
    vec_index_dep = CombIndexToVecIndex(comb_index_dep, n_dep);
    // Convert vec_index_dep to a vector with n.shape()[0]
    // cout << "2";
    for (vector<Index>::size_type i = 0; i < dep_vec.size(); i++)
        vec_index[dep_vec[i]] = vec_index_dep[i];
    // cout << "3" << endl;
    comb_index = VecIndexToCombIndex(vec_index, n);
    return comb_index;
}


Index CombIndexToDepCombIndex(Index comb_index, vector<Index> n_dep, multi_array<Index, 1> n, vector<Index> dep_vec)
{
    Index comb_index_dep;
    multi_array<Index, 1> vec_index({n.shape()[0]});
    vector<Index> vec_index_dep(n_dep.size());
    vec_index = CombIndexToVecIndex(comb_index, n);
    for (vector<Index>::size_type i = 0; i < dep_vec.size(); i++)
    {
        vec_index_dep[i] = vec_index(dep_vec[i]);
    }
    comb_index_dep = VecIndexToCombIndex(vec_index_dep, n_dep);
    return comb_index_dep;
}


Index DepVecIndexRemCombIndexToCombIndex(vector<Index> vec_index_dep, Index comb_index_rem, vector<Index> n_rem, multi_array<Index, 1> n, vector<Index> dep_vec)
{
    vector<Index> vec_index_rem;
    Index comb_index;
    vec_index_rem = CombIndexToVecIndex(comb_index_rem, n_rem);

    // vec_index_c_rem contains now the real population number
    for (vector<Index>::size_type i = 0; i < dep_vec.size(); i++)
        vec_index_rem[dep_vec[i]] = vec_index_dep[i];
    comb_index = VecIndexToCombIndex(vec_index_rem, n);
    return comb_index;
}


void CalculateShiftAmount(std::vector<Index> &sigma1, std::vector<Index> &sigma2, mysys reaction_system, grid_info grid)
{
    Index stride1, stride2;
    Index sigma1_sum, sigma2_sum;

    for (auto &it : reaction_system.reactions)
    {
        stride1 = 1;
        stride2 = 1;
        sigma1_sum = 0;
        sigma2_sum = 0;
        for (Index i = 0; i < grid.m1; i++)
        {
            sigma1_sum += it->nu[i] * stride1 / grid.k1(i);
            stride1 *= grid.n1(i);
        }
        for (Index i = 0; i < grid.m2; i++)
        {
            sigma2_sum += it->nu[i + grid.m1] * stride2 / grid.k2(i);
            stride2 *= grid.n2(i);
        }
        sigma1.push_back(sigma1_sum);
        sigma2.push_back(sigma2_sum);
    }
}


void ShiftMultiArrayRows(int id, multi_array<double, 2> &output_array, const multi_array<double, 2> &input_array, int shift, vector<int> nu, grid_info grid, mysys reaction_system)
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

    if (id == 1)
    {
        grid_alt = new grid_info(grid.m1, grid.m2, grid.r, grid.n1, grid.n2, grid.k1, grid.k2);
        inc = 0;
    }
    else if (id == 2)
    {
        grid_alt = new grid_info(grid.m2, grid.m1, grid.r, grid.n2, grid.n1, grid.k2, grid.k1);
        inc = grid.m1;
    }
    else
    {
        std::cerr << "ERROR: `id` must be 1 (shift in partition 1) or 2 (shift in partition 2)!" << endl;
        std::abort();
    }

    multi_array<Index, 1> vec_index({grid_alt->m1});
    Index k_inc;

    // NOTE: Ensign stores matrices in column-major order
    for (Index j = 0; j < n_cols; j++)
    {
        for (Index i = 0; i < n_rows; i++)
        {
            if ((shift < 0 && i - shift < n_rows) || (shift >= 0 && i - shift >= 0))
            {
                output_array(i, j) = input_array(i - shift, j);
            }
            else
            {
                output_array(i, j) = 0;
                continue;
            }

            vec_index = CombIndexToVecIndex(i, grid_alt->n1);
            for (int k = 0; k < grid_alt->m1; k++)
            {
                k_inc = k + inc;
                if ((nu[k_inc] / grid_alt->k1(k) > 0) &&
                    (vec_index(k) - nu[k_inc] / grid_alt->k1(k) < 0))
                {
                    output_array(i, j) = 0;
                    break;
                }
                else if ((nu[k_inc] / grid_alt->k1(k) < 0) &&
                         (vec_index(k) - nu[k_inc] / grid_alt->k1(k) >= grid_alt->n1(k)))
                {
                    output_array(i, j) = 0;
                    break;
                }
            }
        }
    }
    delete grid_alt;
}