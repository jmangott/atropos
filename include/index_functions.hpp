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
multi_array<double, 1> VecIndexToState(multi_array<Index, 1> vec_index, multi_array<Index, 1> interval, multi_array<double, 1> limit);


// Convert index vector associated with the population numbers to combined index
Index VecIndexToCombIndex(multi_array<Index, 1> vec_index, multi_array<Index, 1> interval);


Index VecIndexToCombIndex(std::vector<Index> vec_index, multi_array<Index, 1> interval);


Index VecIndexToCombIndex(std::vector<Index> vec_index, std::vector<Index> interval);


// Convert combined index to index vector associated with the population numbers
multi_array<Index, 1> CombIndexToVecIndex(Index comb_index, multi_array<Index, 1> interval);


// Convert combined index to index vector associated with the population numbers
std::vector<Index> CombIndexToVecIndex(Index comb_index, std::vector<Index> interval);


inline void CombIndexToState(std::vector<double> &state, Index comb_index, const multi_array<Index, 1> &interval, const multi_array<double, 1> &limit)
{
    Index dim = interval.shape()[0];
    for (Index i = 0; i < dim; i++)
    {
        if (i == (dim - 1))
            state[i] = comb_index * limit(i) / (interval(i) - 1.0);
        else
        {
            state[i] = (comb_index % interval(i)) * limit(i) / (interval(i) - 1.0);
            comb_index = int(comb_index / interval(i));
        }
    }
}


// Convert a combined index for the participating species in reaction mu to a general combined index, where the remaining species have population number 0
Index DepCombIndexToCombIndex(Index comb_index_dep, std::vector<Index> n_dep, multi_array<Index, 1> n, std::vector<Index> dep_vec);


// Convert a general combined index to a combined index for the participating species in reaction mu
Index CombIndexToDepCombIndex(Index comb_index, vector<Index> n_dep, multi_array<Index, 1> n, vector<Index> dep_vec);


// Convert a combined index for the participating and the remaining species in reaction mu to a general combined index
Index DepVecIndexRemCombIndexToCombIndex(std::vector<Index> vec_index_dep, Index comb_index_rem, std::vector<Index> n_rem, multi_array<Index, 1> n, std::vector<Index> dep_vec);


// Calculate for all reactions the number of indices by which arrays have to be shifted in order to calculate the coefficients C1, C2, D1, D2
void CalculateShiftAmount(std::vector<Index> &sigma1, std::vector<Index> &sigma2, mysys reaction_system, grid_info grid);


// Calculate `output_array`, where rows of `input_array` are shifted by `shift`
// NOTE: for positive values of `shift` rows are shifted to larger row indices
void ShiftMultiArrayRows(int id, multi_array<double, 2> &output_array, const multi_array<double, 2> &input_array, int shift, vector<int> nu, grid_info grid, mysys reaction_system);
// void ShiftMultiArrayRows(multi_array<double, 2> &output_array, const multi_array<double, 2> &input_array, int shift);

#endif