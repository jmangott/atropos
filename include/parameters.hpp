#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <algorithm>
#include <vector>

#include <generic/storage.hpp>
#include <lr/lr.hpp>

/////////////////////////////////////////////
//////////////// PARAMETERS /////////////////
/////////////////////////////////////////////

constexpr Index kR = 5;                   // rank
constexpr Index kD = 2;                   // number of species
constexpr Index kN = 51;                  // number of grid points for one species
constexpr Index kK = 1;                   // grid point density
constexpr Index kNsteps = 1000;           // number of time steps
constexpr double kTstar = 1.0;            // final time

// Derived quantities
constexpr Index kM1 = kD / 2;             // number of species in partition 1
constexpr Index kM2 = kD - kM1;           // number of species in partition 2
constexpr double kTau = kTstar / kNsteps; // time step size


/////////////////////////////////////////////
//////// ARRAYS FOR DATUM GENERATION ////////
/////////////////////////////////////////////

// Number of grid points for each species
multi_array<Index, 1> n_xx1({kM1});
multi_array<Index, 1> n_xx2({kM2});

// Inverse of grid spacing (i. e. density of grid points)
multi_array<Index, 1> k_xx1({kM1});
multi_array<Index, 1> k_xx2({kM2});

// Limits for the population number
multi_array<double, 1> lim_xx1({kM1});
multi_array<double, 1> lim_xx2({kM2});

// Grid spacing
multi_array<double, 1> h_xx1({kM1});
multi_array<double, 1> h_xx2({kM2});

// Total number of grid points
Index dxx1_mult;
Index dxx2_mult;

#endif