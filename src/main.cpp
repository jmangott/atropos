#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "reactions.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;

int main()
{
    constexpr int r = 10; // desired rank
    constexpr size_t d = 2;
 
    std::fill_n(std::back_inserter(xx), d, 1.0);
    nn = {"A1", "B2"};

    double tstar = 20; // final time
    double tau = 0.00625; // time step size
    Index nsteps = tstar / tau;

    array<Index, d> n_xx = {100, 100}; // number of grid points for each species population
    array<double, d> lim_xx = {10.0, 10.0}; // limits for the population number


    int m = 1;
    double x_max = 10.0;
    int n = 100; 
    // multi_array<double, 2> p({m, d - m});

    // Set up the probability for t = 0


    return 0;
}