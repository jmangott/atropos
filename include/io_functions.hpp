#ifndef IO_FUNCTIONS_HPP
#define IO_FUNCTIONS_HPP

#include <fstream>
#include <iostream>
#include <sstream>

#include <generic/storage.hpp>


// Read in multi_array from a .csv file
void ReadInMultiArray(multi_array<double, 2> &output_array, string filename);

#endif