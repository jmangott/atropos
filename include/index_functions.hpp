#ifndef INDEX_FUNCTIONS_HPP
#define INDEX_FUNCTIONS_HPP

#include <iostream>
#include <stdexcept>
#include <vector>

#include <generic/storage.hpp>

#include "grid_class.hpp"
#include "reaction_class.hpp"

// TODO: Create CombIndexToState (?)


// Convert index vector to state vector
multi_array<double, 1> VecIndexToState(multi_array<Index, 1> vec_index, multi_array<Index, 1> interval, multi_array<double, 1> limit);


// Convert index vector associated with the population numbers to combined index
Index VecIndexToCombIndex(multi_array<Index, 1> vec_index, multi_array<Index, 1> interval);


// Convert combined index to index vector associated with the population numbers
multi_array<Index, 1> CombIndexToVecIndex(Index comb_index, multi_array<Index, 1> interval);


// Convert combined index to index vector associated with the population numbers
std::vector<Index> CombIndexToVecIndex(Index comb_index, std::vector<Index> interval);


// Calculate for all reactions the number of indices by which arrays have to be shifted in order to calculate the coefficients C1, C2, D1, D2
template<Index m1, Index m2>
void CalculateShiftAmount(std::vector<Index> &sigma1, std::vector<Index> &sigma2, mysys reaction_system, grid_info<m1, m2> grid)
{
    Index stride1, stride2;
    Index sigma1_sum, sigma2_sum;

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
            sigma1_sum += it->nu[i] * stride1 * grid.k1(i);
            stride1 *= grid.n1(i);
        }
        for (Index i = 0; i < m2; i++)
        {
            sigma2_sum += it->nu[i + m1] * stride2 * grid.k2(i);
            stride2 *= grid.n2(i);
        }
        sigma1.push_back(sigma1_sum);
        sigma2.push_back(sigma2_sum);
    }
}


// Calculate `output_array`, where rows of `input_array` are shifted by `shift`
void ShiftMultiArrayRows(multi_array<double, 2> &output_array, multi_array<double, 2> &input_array, int shift);

#endif