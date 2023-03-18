#include "grid_class.hpp"
#include "index_functions.hpp"

using std::vector;


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


Index DepCombIndexToCombIndex(Index comb_index_dep, vector<Index> n_dep, multi_array<Index, 1> n, vector<Index> dep_vec)
{
    vector<Index> vec_index(n.shape()[0], 0);
    vector<Index> vec_index_dep(n_dep.size());
    Index comb_index;
    CombIndexToVecIndex(vec_index_dep, comb_index_dep, n_dep);
    // Convert vec_index_dep to a vector with n.shape()[0]
    for (vector<Index>::size_type i = 0; i < dep_vec.size(); i++)
        vec_index[dep_vec[i]] = vec_index_dep[i];
    comb_index = VecIndexToCombIndex(vec_index, n);
    return comb_index;
}


Index DepVecIndexRemCombIndexToCombIndex(vector<Index> vec_index_dep, Index comb_index_rem, vector<Index> n_rem, multi_array<Index, 1> n, vector<Index> dep_vec)
{
    vector<Index> vec_index_rem(n_rem.size());
    Index comb_index;
    CombIndexToVecIndex(vec_index_rem, comb_index_rem, n_rem);

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

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        stride1 = 1;
        stride2 = 1;
        sigma1_sum = 0;
        sigma2_sum = 0;
        for (Index i = 0; i < grid.m1; i++)
        {
            sigma1_sum += reaction_system.reactions[mu]->nu[i] * stride1;
            stride1 *= grid.n1(i);
        }
        for (Index i = 0; i < grid.m2; i++)
        {
            sigma2_sum += reaction_system.reactions[mu]->nu[i + grid.m1] * stride2;
            stride2 *= grid.n2(i);
        }
        sigma1[mu] = sigma1_sum;
        sigma2[mu] = sigma2_sum;
    }
}