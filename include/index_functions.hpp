#ifndef INDEX_FUNCTIONS_HPP
#define INDEX_FUNCTIONS_HPP

#include <iostream>
#include <stdexcept>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "reaction_class.hpp"

// Convert index vector to state vector
multi_array<double, 1> VecIndexToState(multi_array<Index, 1> vec_index, multi_array<Index, 1> interval, multi_array<double, 1> limit);


// Convert index vector associated with the population numbers to combined index
Index VecIndexToCombIndex(multi_array<Index, 1> vec_index, multi_array<Index, 1> interval);


// Convert combined index to index vector associated with the population numbers
multi_array<Index, 1> CombIndexToVecIndex(Index comb_index, multi_array<Index, 1> interval);


// Convert combined index to index vector associated with the population numbers
std::vector<Index> CombIndexToVecIndex(Index comb_index, std::vector<Index> interval);


// Calculate for all reactions the number of indices by which arrays have to be shifted in order to calculate the coefficients C1, C2, D1, D2
void CalculateShiftAmount(vector<Index> &sigma1, std::vector<Index> &sigma2, mysys mysystem, multi_array<Index, 1> n_xx1, multi_array<Index, 1> n_xx2, multi_array<Index, 1> k_xx1, multi_array<Index, 1> k_xx2);

// Calculate `output_array`, where rows of `input_array` are shifted by `shift`
void ShiftMultiArrayCols(multi_array<double, 2> &output_array, multi_array<double, 2> &input_array, int shift);

#endif