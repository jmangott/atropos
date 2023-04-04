#ifndef IO_FUNCTIONS_HPP
#define IO_FUNCTIONS_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <netcdf.h>

#include <generic/netcdf.hpp>
#include <generic/storage.hpp>

#include <lr/lr.hpp>

#include "grid_class.hpp"
#include "reaction_class.hpp"


// Read multi_array from a .csv file
void ReadCSV(multi_array<double, 2> &output_array, std::string filename);


// Write multi_array to a .csv file
void WriteCSV(const multi_array<double, 2> &input_array, std::string filename);


// Write `lr_sol` to a .netCDF file
void WriteNC(std::string fn, const lr2<double> &lr_sol, vector<string> names, grid_info grid, double *t, double *dt);

// Read a .netCDF file to multi_arrays `xx1` and `xx2`
void ReadNC(std::string fn, multi_array<double, 2> &xx1, multi_array<double, 2> &xx2);

#endif