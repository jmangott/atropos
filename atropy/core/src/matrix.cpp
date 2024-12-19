#include "matrix.hpp"

// Shift operator
template <>
void Matrix::ShiftRows<1>(Ensign::multi_array<double, 2>& output_array,
                          const Ensign::multi_array<double, 2>& input_array,
                          const grid_parms grid, const Index mu);

// Inverse shift operator
template <>
void Matrix::ShiftRows<-1>(Ensign::multi_array<double, 2>& output_array,
                           const Ensign::multi_array<double, 2>& input_array,
                           const grid_parms grid, const Index mu);
