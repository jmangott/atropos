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

    Index dxx1_mult = n_xx[0];
    Index dxx2_mult = n_xx[1];

    vector<const double *> x1, x2;

    // Check if reaction class works
    // cout << mysystem.reactions[0]->propensity() << endl;

    // For coefficients
    multi_array<double, 2> C1v({r, r});
    multi_array<double, 2> C1w({r, r});

    multi_array<double, 2> C2v({r, r});
    multi_array<double, 2> C2w({r, r});

    multi_array<complex<double>, 2> C2vc({r, r});
    multi_array<complex<double>, 2> C2wc({r, r});

    multi_array<double, 1> we_v({dxx2_mult});
    multi_array<double, 1> we_w({dxx2_mult});

    // Temporary objects for multiplication
    multi_array<double, 2> tmpX({dxx1_mult, r});

    // Set up S, X1 and X2h for t = 0
    multi_array<double, 1> ss({r});
    multi_array<double, 2> xx1({dxx1_mult, r});
    multi_array<double, 2> xx2h({dxx2_mult, r});

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

    x1.push_back(xx1.begin());
    x2.push_back(xx2h.begin());
    
    // Set up the low-rank structure and the inner products
    lr2<double> lr_sol(r, {dxx1_mult, dxx2_mult});

    // TODO: check if inner product is set up correctly
    std::function<double(double *, double *)> ip_xx1 = inner_product_from_const_weight(h_xx[0], dxx1_mult);
    std::function<double(double *, double *)> ip_xx2 = inner_product_from_const_weight(h_xx[1], dxx2_mult);

    blas_ops blas;
    initialize(lr_sol, x1, x2, ip_xx1, ip_xx2, blas);

    /////////////////////////////////////////////
    ////////////////// K-STEP ///////////////////
    /////////////////////////////////////////////

    tmpX = lr_sol.X;
    blas.matmul(tmpX, lr_sol.S, lr_sol.X);

    // TODO:
    // Insert here function for calculating the second argument of coeff(lr_sol.V, ..., we_v, etc.)
    //                                                                              ^

    coeff(lr_sol.V, lr_sol.V, we_v, C1v, blas);
    coeff(lr_sol.V, lr_sol.V, we_w, C1w, blas);

    return 0;
}