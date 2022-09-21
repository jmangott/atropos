#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "reactions.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;
using std::vector;

int main()
{
    constexpr int r = 5; // rank
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

    // cout << mysystem.reactions[0]->propensity() << endl;

    // Set up S, X1 and X2h for t = 0
    multi_array<double, 1> ss({r});
    multi_array<double, 2> xx1({n_xx[0], r});
    multi_array<double, 2> xx2h({n_xx[1], r});

    string line;

    // Read in S
    ifstream s_input_file("../input/s.csv");
    int ii = 0;
    while (getline (s_input_file, line))
    {
        ss(ii) = stod(line);
        ii++;
    }
    s_input_file.close();

    // Read in X1
    ii = 0;
    int jj = 0;
    stringstream ssline;
    string element;

    ifstream xx1_input_file("../input/u.csv");
    while (getline (xx1_input_file, line))
    {
        ssline.str(line);
        while (getline(ssline, element, ','))
        {
            xx1(ii, jj) = stod(element);
            jj++;
        }
        ii++;
        jj = 0;
        ssline.clear();
    }
    xx1_input_file.close();

    // Read in X2
    ii = 0;
    jj = 0;

    ifstream xx2h_input_file("../input/vh.csv");
    while (getline (xx2h_input_file, line))
    {
        ssline.str(line);
        while (getline(ssline, element, ','))
        {
            xx2h(ii, jj) = stod(element);
            jj++;
        }
        ii++;
        jj = 0;
        ssline.clear();
    }
    xx2h_input_file.close();

    return 0;
}