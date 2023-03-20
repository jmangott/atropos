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

#include "models/reactions_tyson.hpp"


/////////////////////////////////////////////
//////////////// PARAMETERS /////////////////
/////////////////////////////////////////////

constexpr Index kR = 6;                        // rank
constexpr Index kD = 5;                        // total number of species
constexpr Index kM1 = 3;                       // number of species in partition 1
constexpr Index kM2 = 2;                       // number of species in partition 2

std::vector<Index> kN1 {21, 61, 11};           // number of grid points for species in partition 1
std::vector<Index> kN2 {71, 61};               // number of grid points for species in partition 2
std::vector<Index> kBinsize1 {1, 1, 1};        // binsize for species in partition 1
std::vector<Index> kBinsize2 {1, 1};           // binsize for species in partition 2
std::vector<double> kLiml1 {0.0, 0.0, 0.0};    // left population number limit for species in partition 1
std::vector<double> kLiml2 {0.0, 0.0};         // left population number limit for species in partition 2

constexpr double kTstar = 40.00;               // final time
double kTau = 0.00001;                         // time step size
Index kNsteps = ceil(kTstar / kTau);           // number of time steps
constexpr Index kSnapshot = 100;               // number of time steps between snapshots

constexpr Index kNBasisFunctions = 1;          // number of basis functions for the initial condition,
                                               // must be > 0 and <= kR!

constexpr Index kNSubsteps = 1000;             // number of explicit Euler substeps for one 
                                               // integration step for second order method

std::string kFilename = "tyson_small";         // name of the output folder

constexpr bool kPrintDiagnostics = true;       // if `true`, diagnostics (e.g. memory requirements, 
                                               // maximum propensity value) will be displayed

#endif