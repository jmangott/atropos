#ifndef IO_FUNCTIONS_HPP
#define IO_FUNCTIONS_HPP

#include <fstream>
#include <iostream>
#include <sstream>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

void ReadInMultiArray(multi_array<double, 2> &output_array, string filename);

#endif