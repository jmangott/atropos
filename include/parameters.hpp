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

// // Toggle switch
// vector<Index> kN1{51};
// vector<Index> kN2{51};
// vector<Index> kK1{1};
// vector<Index> kK2{1};

// Lambda phage
vector<Index> kN1 {16, 41};
vector<Index> kN2 {11, 11, 11};
vector<Index> kK1 {1, 1};
vector<Index> kK2 {1, 1, 1};

constexpr double kTstar = 10;                     // final time
double kTau = 0.01;                               // time step size
Index kNsteps = ceil(kTstar / kTau);              // number of time steps
constexpr Index kSnapshot = 100;                  // number of time steps between snapshots

#endif