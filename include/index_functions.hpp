#ifndef INDEX_FUNCTIONS_HPP
#define INDEX_FUNCTIONS_HPP

#include <iostream>
#include <stdexcept>
#include <vector>

#include <generic/storage.hpp>

#include "grid_class.hpp"
#include "reaction_class.hpp"

// TODO: Rewrite all index functions in the same style as `CombIndexToState`

// Convert index vector to state vector
multi_array<double, 1> VecIndexToState(multi_array<Index, 1> vec_index, multi_array<Index, 1> interval, multi_array<double,1 > liml, multi_array<double, 1> limr);


// Convert index vector associated with the population numbers to combined index
Index VecIndexToCombIndex(multi_array<Index, 1> vec_index, multi_array<Index, 1> interval);


Index VecIndexToCombIndex(std::vector<Index> vec_index, multi_array<Index, 1> interval);


Index VecIndexToCombIndex(std::vector<Index> vec_index, std::vector<Index> interval);


// inline void CombIndexToVecIndex(multi_array<Index, 1> &vec_index, Index comb_index, const multi_array<Index, 1> &interval)
// {
//     Index dim = interval.shape()[0];
//     for (Index i = 0; i < dim; i++)
//     {
//         if (i == (dim - 1))
//             vec_index(i) = comb_index;
//         else
//         {
//             vec_index(i) = comb_index % interval(i);
//             comb_index = int(comb_index / interval(i));
//         }
//     }
// }


inline void CombIndexToVecIndex(multi_array<Index, 1> &vec_index, Index comb_index, const multi_array<Index, 1> &interval)
{
    Index dim = interval.shape()[0];
    for (Index i = 0; i < (dim - 1); i++)
    {
        vec_index(i) = comb_index % interval(i);
        comb_index = int(comb_index / interval(i));
    }
    if (dim > 0) vec_index(dim - 1) = comb_index;
}


// inline void CombIndexToVecIndex(std::vector<Index> &vec_index, Index comb_index, const std::vector<Index> &interval)
// {
//     Index dim = interval.size();
//     for (Index i = 0; i < dim; i++)
//     {
//         if (i == (dim - 1))
//             vec_index[i] = comb_index;
//         else
//         {
//             vec_index[i] = comb_index % interval[i];
//             comb_index = int(comb_index / interval[i]);
//         }
//     }
// }


inline void CombIndexToVecIndex(std::vector<Index> &vec_index, Index comb_index, const std::vector<Index> &interval)
{
    Index dim = interval.size();
    for (Index i = 0; i < (dim - 1); i++)
    {
        vec_index[i] = comb_index % interval[i];
        comb_index = int(comb_index / interval[i]);
    }
    if (dim > 0) vec_index[dim - 1] = comb_index;
}


inline void CombIndexToState(std::vector<double> &state, Index comb_index, const multi_array<Index, 1> &interval, const multi_array<double, 1> &liml, const multi_array<double, 1> &limr)
{
    Index dim = interval.shape()[0];
    for (Index i = 0; i < dim; i++)
    {
        if (i == (dim - 1))
            state[i] = liml(i) + (limr(i) - liml(i)) * comb_index / (interval(i) - 1.0);
        else
        {
            state[i] = liml(i) + (limr(i) - liml(i)) * (comb_index % interval(i)) / (interval(i) - 1.0);
            comb_index = int(comb_index / interval(i));
        }
    }
}


// Convert a combined index for the participating species in reaction mu to a general combined index, where the remaining species have population number 0
Index DepCombIndexToCombIndex(Index comb_index_dep, std::vector<Index> n_dep, multi_array<Index, 1> n, std::vector<Index> dep_vec);


// Convert a general combined index to a combined index for the participating species in reaction mu
inline Index CombIndexToDepCombIndex(Index comb_index, const vector<Index> &n_dep, const multi_array<Index, 1> &n, const vector<Index> &dep_vec)
{
    Index comb_index_dep = 0;
    Index stride = 1;
    multi_array<Index, 1> vec_index(n.shape());
    vector<Index> vec_index_dep(n_dep.size());
    CombIndexToVecIndex(vec_index, comb_index, n);
    for (vector<Index>::size_type i = 0; i < dep_vec.size(); i++)
    {
        comb_index_dep += vec_index(dep_vec[i]) * stride;
        stride *= n_dep[i];
    }
    return comb_index_dep;
}


// Convert a combined index for the participating and the remaining species in reaction mu to a general combined index
Index DepVecIndexRemCombIndexToCombIndex(std::vector<Index> vec_index_dep, Index comb_index_rem, std::vector<Index> n_rem, multi_array<Index, 1> n, std::vector<Index> dep_vec);


// Calculate for all reactions the number of indices by which arrays have to be shifted in order to calculate the coefficients C1, C2, D1, D2
void CalculateShiftAmount(std::vector<Index> &sigma1, std::vector<Index> &sigma2, mysys reaction_system, grid_info grid);


// Calculate `output_array`, where rows of `input_array` are shifted by `shift`
// NOTE: for positive values of `shift` rows are shifted to larger row indices
void ShiftMultiArrayRows(int id, multi_array<double, 2> &output_array, const multi_array<double, 2> &input_array, Index shift, vector<int> nu, grid_info grid);
// void ShiftMultiArrayRows(multi_array<double, 2> &output_array, const multi_array<double, 2> &input_array, int shift);

#endif