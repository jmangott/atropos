#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <algorithm>
#include <string>
#include <vector>

#include <generic/storage.hpp>
#include <lr/lr.hpp>

/////////////////////////////////////////////
/////////////////// MODEL ///////////////////
/////////////////////////////////////////////

// #include "models/reactions_ts.hpp"
// #include "models/reactions_lp.hpp"
// #include "models/reactions_tgfb6.hpp"
#include "models/reactions_tyson.hpp"


/////////////////////////////////////////////
//////////////// PARAMETERS /////////////////
/////////////////////////////////////////////

// // Toggle switch
// constexpr Index kR = 5;
// constexpr Index kD = 2;
// constexpr Index kM1 = 1;
// constexpr Index kM2 = 1;

// std::vector<Index> kN1{51};
// std::vector<Index> kN2{51};
// std::vector<Index> kK1{1};
// std::vector<Index> kK2{1};
// std::vector<double> kLiml1{0.0};
// std::vector<double> kLiml2{0.0};

// // Toggle switch new
// constexpr Index kR = 5;
// constexpr Index kD = 2;
// constexpr Index kM1 = 1;
// constexpr Index kM2 = 1;

// std::vector<Index> kN1 {301};
// std::vector<Index> kN2 {301};
// std::vector<Index> kK1 {1};
// std::vector<Index> kK2 {1};
// std::vector<double> kLiml1 {0.0};
// std::vector<double> kLiml2 {0.0};

// // Lambda phage
// constexpr Index kR = 8;
// constexpr Index kD = 5;
// constexpr Index kM1 = 2;
// constexpr Index kM2 = 3;

// std::vector<Index> kN1 {16, 41};
// std::vector<Index> kN2 {11, 11, 11};
// std::vector<Index> kK1 {1, 1};
// std::vector<Index> kK2 {1, 1, 1};
// std::vector<double> kLiml1 {0.0, 0.0};
// std::vector<double> kLiml2 {0.0, 0.0, 0.0};

// // Lambda phage new
// constexpr Index kR = 4;
// constexpr Index kD = 5;
// constexpr Index kM1 = 2;
// constexpr Index kM2 = 3;

// // std::vector<Index> kN1 {6, 151};
// // std::vector<Index> kN2 {6, 11, 21};
// std::vector<Index> kN1 {6, 41};
// std::vector<Index> kN2 {6, 11, 21};
// std::vector<Index> kK1 {1, 1};
// std::vector<Index> kK2 {1, 1, 1};
// // std::vector<double> kLiml1 {0.0, 0.0};
// std::vector<double> kLiml1 {0.0, 70.0};
// std::vector<double> kLiml2 {0.0, 0.0, 0.0};

// // TGFb6
// constexpr Index kR = 4;
// constexpr Index kD = 8;
// constexpr Index kM1 = 4;
// constexpr Index kM2 = 4;

// std::vector<Index> kN1 {5, 5, 151, 151};
// std::vector<Index> kN2 {26, 21, 21, 21};
// std::vector<Index> kK1 {1, 1, 1, 1};
// std::vector<Index> kK2 {1, 1, 1, 1};
// std::vector<double> kLiml1 {333.0, 0.0, 470.0, 20.0};
// std::vector<double> kLiml2 {1890.0, 0.0, 20.0, 0.0};

// Tyson
constexpr Index kR = 6;                        // rank
constexpr Index kD = 5;                        // number of species
constexpr Index kM1 = 2;                       // number of species in partition 1
constexpr Index kM2 = 3;                       // number of species in partition 2

std::vector<Index> kN1 {1521, 1901};           // number of grid points for species in partition 1
std::vector<Index> kN2 {16, 2041, 1251};       // number of grid points for species in partition 2
std::vector<Index> kK1 {1, 1};
std::vector<Index> kK2 {1, 1, 1};
std::vector<double> kLiml1 {0.0, 0.0};         // left population number limit for species in partition 1
std::vector<double> kLiml2 {0.0, 4000.0, 0.0}; // left population number limit for species in partition 2

constexpr double kTstar = 0.01;                // final time
double kTau = 0.01;                            // time step size
Index kNsteps = ceil(kTstar / kTau);           // number of time steps
constexpr Index kSnapshot = 1;                 // number of time steps between snapshots

constexpr Index kNBasisFunctions = 1;          // number of basis functions for the initial condition,
                                               // must be > 0 and <= kR!

constexpr Index kNSubsteps = 100;              // number of explicit Euler substeps for one 
                                               // integration step for second order method

std::string kFilename = "tgfb6";               // name of the output folder

constexpr bool kPrintDiagnostics = true;       // if `true`, diagnositcs (e.g. memory requirements, 
                                               // maximum propensity value) value will be displayed

#endif