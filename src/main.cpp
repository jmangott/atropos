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

    Index nsteps = 1000; // # time steps
    double tstar = 1.0; // final time
    double tau = tstar / nsteps; // time step size

    array<Index, d> n_xx = {50, 50}; // number of grid points for each species population
    array<double, d> lim_xx = {10.0, 10.0}; // limits for the population number


    int m1 = d / 2;
    int m2 = d - m1;
    double x_max = 10.0;
    int n = 100;

    // Initial datum generation
    array<double, 2> h_xx;
    for (int ii = 0; ii < 2; ii++)
    {
        h_xx[ii] = lim_xx[ii] / (n_xx[ii] - 1);
    }



    // NOTE: this quantity should be never formed in the real implementation!
    multi_array<double, 2> cc({d, d});
    multi_array<double, 2> p({m1, m2});

    cout << mysystem.reactions[0]->propensity() << endl;

    // Set up the probability for t = 0

    return 0;
}