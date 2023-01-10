#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <algorithm>
#include <vector>

#include <generic/storage.hpp>
#include <lr/lr.hpp>

/////////////////////////////////////////////
//////////////// PARAMETERS /////////////////
/////////////////////////////////////////////

// // Toggle switch
// constexpr Index kR = 5;                           // rank
// constexpr Index kD = 2;                           // number of species
// constexpr Index kM1 = 1;                          // number of species in partition 1
// constexpr Index kM2 = 1;                          // number of species in partition 2

// Lambda phage
constexpr Index kR = 4;                           // rank
constexpr Index kD = 5;                           // number of species
constexpr Index kM1 = 2;                          // number of species in partition 1
constexpr Index kM2 = 3;                          // number of species in partition 2

// constexpr Index kN = 51;                       // number of grid points for one species
multi_array<Index, 1> kN1({kM1});
multi_array<Index, 1> kN2({kM2});
// constexpr Index kK = 1;                        // grid point density
multi_array<Index, 1> kK1({kM1});
multi_array<Index, 1> kK2({kM2});
constexpr double kTstar = 10;                     // final time
double kTau = 0.01;                               // time step size
Index kNsteps = ceil(kTstar / kTau);              // number of time steps
constexpr Index kSnapshot = 100;                  // number of time steps between snapshots

#endif