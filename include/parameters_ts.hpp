#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <algorithm>
#include <vector>

#include <generic/storage.hpp>
#include <lr/lr.hpp>

/////////////////////////////////////////////
//////////////// PARAMETERS /////////////////
/////////////////////////////////////////////

constexpr Index kR = 5;                     // rank
constexpr Index kD = 2;                     // number of species
constexpr Index kN = 51;                    // number of grid points for one species
constexpr Index kK = 1;                     // grid point density
constexpr double kTstar = 20;              // final time
double kTau = 0.001;                         // time step size

// Derived quantities
constexpr Index kM1 = kD / 2;               // number of species in partition 1
constexpr Index kM2 = kD - kM1;             // number of species in partition 2
Index kNsteps = ceil(kTstar/kTau);          // number of time steps

#endif