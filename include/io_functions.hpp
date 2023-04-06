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

// Write `lr_sol` to a .netCDF file
void WriteNC(std::string fn, const lr2<double> &lr_sol, vector<string> names, grid_info grid, double *t, double *dt);


// Read a .netCDF file to multi_arrays `xx1`, `xx2` and `ss`
void ReadNC(std::string fn, multi_array<double, 2> &xx1, multi_array<double, 2> &xx2, multi_array<double, 2> &ss, Index &n_basisfunctions);

#endif