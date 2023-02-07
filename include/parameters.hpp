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

// std::vector<Index> kN1{51};
// std::vector<Index> kN2{51};
// std::vector<Index> kK1{1};
// std::vector<Index> kK2{1};
// std::vector<double> kLiml1{0.0};
// std::vector<double> kLiml2{0.0};

// // Lambda phage
// constexpr Index kR = 8;                           // rank
// constexpr Index kD = 5;                           // number of species
// constexpr Index kM1 = 2;                          // number of species in partition 1
// constexpr Index kM2 = 3;                          // number of species in partition 2

// std::vector<Index> kN1 {16, 41};
// std::vector<Index> kN2 {11, 11, 11};
// std::vector<Index> kK1 {1, 1};
// std::vector<Index> kK2 {1, 1, 1};
// std::vector<double> kLiml1 {0.0, 0.0};
// std::vector<double> kLiml2 {0.0, 0.0, 0.0};

// TGFb6
constexpr Index kR = 5;                            // rank
constexpr Index kD = 8;                            // number of species
constexpr Index kM1 = 6;                           // number of species in partition 1
constexpr Index kM2 = 2;                           // number of species in partition 2

std::vector<Index> kN1 {5, 5, 21, 21, 26, 21};
std::vector<Index> kN2 {151, 151};
std::vector<Index> kK1 {1, 1, 1, 1, 1, 1};
std::vector<Index> kK2 {1, 1};
std::vector<double> kLiml1 {333.0, 0.0, 20.0, 0.0, 1890.0, 0.0};
std::vector<double> kLiml2 {470.0, 20.0};

constexpr double kTstar = 2.0;                    // final time
double kTau = 0.01;                               // time step size
Index kNsteps = ceil(kTstar / kTau);              // number of time steps
constexpr Index kSnapshot = 10;                   // number of time steps between snapshots

#endif