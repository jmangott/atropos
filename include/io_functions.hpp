#ifndef IO_FUNCTIONS_HPP
#define IO_FUNCTIONS_HPP

#include <fstream>
#include <iostream>
#include <sstream>

#include <generic/netcdf.hpp>
#include <generic/storage.hpp>

#include <lr/lr.hpp>

#include "grid_class.hpp"


// Read in multi_array from a .csv file
void ReadInMultiArray(multi_array<double, 2> &output_array, string filename);


// Write out multi_array to a .csv file
void WriteOutMultiArray(const multi_array<double, 2> &input_array, string filename);

#endif