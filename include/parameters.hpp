#include "models/reactions_tyson.hpp"

constexpr Index kR = 6;  // rank
constexpr Index kM1 = 3;  // number of species in partition 1
constexpr Index kM2 = 2;  // number of species in partition 2

constexpr Index kN1[kM1] {21, 61, 11};  // number of grid points for species in partition 1
constexpr Index kN2[kM2] {71, 61};  // number of grid points for species in partition 2
constexpr Index kBinsize1[kM1] {1, 1, 1};  // binsize for species in partition 1
constexpr Index kBinsize2[kM2] {1, 1};  // binsize for species in partition 2
constexpr double kLiml1[kM1] {0, 0, 0};  // left population number limit for species in partition 1
constexpr double kLiml2[kM2] {0, 0};  // left population number limit for species in partition 2

constexpr double kTstar = 1.0;  // final time
double kTau = 0.01;  // time step size
constexpr Index kSnapshot = 100;  // number of time steps between snapshots
constexpr bool kSecondOrder = true;  // flag for activating the second order method
constexpr Index kNSubsteps = 10;  // number of explicit Euler substeps for one integration step for second order method
constexpr char kFilename[] = "tyson1";  // name of the output folder