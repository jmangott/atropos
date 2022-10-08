#include <algorithm>
#include <fstream>
#include <iostream>
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

// TODO: improve speed by using a lookup table?
vector<double> index_to_state(Index idx, vector<Index> interval)
{
    int stride = 1;
    size_t dim = interval.size();
    vector<double> state_vec(dim, 0);
    for (size_t k = 0; k < dim; k++)
    {
        if (k == (dim - 1))
        {
            state_vec[k] = (int(idx / stride));
        }
        else
        {
            state_vec[k] = (int((idx % interval[k]) / stride));
            stride *= interval[k];
        }
    }
    return state_vec;
}

int main()
{
    /////////////////////////////////////////////
    ////// DECLARATION AND INITIALIZATION ///////
    /////////////////////////////////////////////

    constexpr int r = 5; // rank
    constexpr int d = 2; // # species
 
    std::fill_n(std::back_inserter(xx), d, 1.0);
    nn = {"A1", "B2"};

    // Index nsteps = 1000; // # time steps
    // double tstar = 1.0; // final time
    // double tau = tstar / nsteps; // time step size

    constexpr size_t m1 = d / 2;
    constexpr size_t m2 = d - m1;

    // Number of grid points for each species in partition 1 and 2
    vector<Index> n_xx1(m1, 51);
    vector<Index> n_xx2(m2, 51);

    // Limits for the population number in partition 1 and 2
    vector<int> lim_xx1(m1, 50);
    vector<int> lim_xx2(m2 ,50);

    // Initial datum generation
    vector<double> h_xx1(m1, 0);
    vector<double> h_xx2(m2, 0);
    vector<int> k_xx1(m1, 0);
    vector<int> k_xx2(m2, 0);
    for (size_t ii = 0; ii < m1; ii++)
    {
        h_xx1[ii] = double(lim_xx1[ii]) / (n_xx1[ii] - 1);
        k_xx1[ii] = int((n_xx1[ii] - 1) / lim_xx1[ii]);
    }

    for (size_t ii = 0; ii < m2; ii++)
    {
        h_xx2[ii] = lim_xx2[ii] / (n_xx2[ii] - 1);
        k_xx2[ii] = int((n_xx2[ii] - 1) / lim_xx2[ii]);
    }

    Index dxx1_mult = 1;
    Index dxx2_mult = 1;
    for (auto ele : n_xx1) dxx1_mult *= ele;
    for (auto ele : n_xx2) dxx2_mult *= ele;

    vector<const double *> x1, x2;

    // For coefficients
    multi_array<double, 2> c2({r, r});
    multi_array<double, 2> d2({r, r});

    multi_array<double, 1> w_x2({dxx2_mult});

    // Temporary objects for multiplication
    multi_array<double, 2> tmp_x({dxx1_mult, r});

    // Set up S, X1 and X2h for t = 0
    multi_array<double, 2> ss({r, r});
    multi_array<double, 2> xx1({dxx1_mult, r});
    multi_array<double, 2> xx2({dxx2_mult, r});


    /////////////////////////////////////////////
    /////////// READ IN S, X1 and X2 ////////////
    /////////////////////////////////////////////

    // TODO: Use .netcdf instead of .csv

    string line;

    // Read in S
    ifstream s_input_file("../input/s.csv");
    int ii = 0;
    int jj = 0;
    stringstream ssline;
    string element;

    while (getline(s_input_file, line))
    {
        ssline.str(line);
        while (getline(ssline, element, ','))
        {
            ss(ii, jj) = stod(element);
            jj++;
        }
        ii++;
        jj = 0;
        ssline.clear();
    }
    s_input_file.close();

    // Read in X1
    ii = 0;
    jj = 0;

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

    ifstream xx2_input_file("../input/vh.csv");
    while (getline (xx2_input_file, line))
    {
        ssline.str(line);
        while (getline(ssline, element, ','))
        {
            xx2(ii, jj) = stod(element);
            jj++;
        }
        ii++;
        jj = 0;
        ssline.clear();
    }
    xx2_input_file.close();

    double *it1 = xx1.begin();
    double *it2 = xx2.begin();
    for (size_t i = 0; i < r; i++)
    {
        x1.push_back(it1);
        x2.push_back(it2);
        it1 += m1 * n_xx1[0];
        it2 += m2 * n_xx2[0];
    }
    
    // Set up the low-rank structure and the inner products
    lr2<double> lr_sol(r, {dxx1_mult, dxx2_mult});

    std::function<double(double *, double *)> ip_xx1 = inner_product_from_const_weight(h_xx1[0], dxx1_mult);
    std::function<double(double *, double *)> ip_xx2 = inner_product_from_const_weight(h_xx2[0], dxx2_mult);

    blas_ops blas;
    initialize(lr_sol, x1, x2, ip_xx1, ip_xx2, blas);
    lr_sol.S = ss;


    /////////////////////////////////////////////
    ////////////////// K-STEP ///////////////////
    /////////////////////////////////////////////

    tmp_x = lr_sol.X;
    blas.matmul(tmp_x, lr_sol.S, lr_sol.X); // lr_sol.X contains now K
    
    vector<double> sigma1, sigma2;
    double sigma2_sum;
    
    // NOTE: when the partition requires a permutation of the original order of species,
    // then also the nu vectors and similar quantities have to be permuted

    for(auto it : mysystem.reactions) {
        sigma2_sum = 0;
        for (size_t i = 0; i < m2; i++)
        {
            sigma2_sum =+ it->nu[i] * pow(n_xx2[i], i);
        }
        sigma2.push_back(sigma2_sum);
    }


    // Calculate C^2 for given mu and alpha_tilde1

    // TODO: rewrite the whole section as function, e.g. calcC2(mu, alpha_tilde1)
    // (or independent of alpha_tilde1)

    int mu = 2;
    int alpha_tilde1 = 30;

    // Calculate the shifted X1

    // TODO: write this shifting procedure as a function
    vector<const double *> x2_shift = x2;
    for (size_t i = 0; i < r; i++)
    {
        x2_shift[i] += int(k_xx2[0] * sigma2[mu]);
    }

    multi_array<double, 2> xx2_shift({dxx2_mult, r});
    for (int j = 0; j < r; j++) 
    {
        for (Index i = 0; i < dxx2_mult; i++)
        {
            if ((k_xx2[0] * sigma2[mu] < 0) && 
                (i < -k_xx2[0] * sigma2[mu]))
            {
                xx2_shift(i, j) = x2[j][0];
            }
            else if ((k_xx2[0] * sigma2[mu] > 0) && 
                     (i >= (dxx2_mult - k_xx2[0] * sigma2[mu])))
            {
                xx2_shift(i, j) = x2[j][dxx2_mult - 1];
            }
            else
            {
                xx2_shift(i, j) = x2_shift[j][i];
            }
        }
    }

    // for (int i = 0; i < dxx2_mult; i++)
    // {
    //     for (int j = 0; j < r; j++)
    //     {
    //         cout << xx2_shift(i, j) << " ";
    //     }
    //     cout << endl;
    // }

    vector<size_t> dep_vec, dep_vec1, dep_vec2;
    dep_vec = mysystem.reactions[0]->depends_on;
    for (auto & ele : dep_vec)
    {
        (ele < m1) ? dep_vec1.push_back(ele) : dep_vec2.push_back(ele); 
    }

    // int stride = 1;
    vector<double> a_vec(dxx2_mult, 0);
    vector<double> n_tot(d, 0);
    vector<double> n1(d, 0);
    vector<double> n2(d, 0);

    n1 = index_to_state(alpha_tilde1, n_xx1);
    for (size_t i = 0; i < n1.size(); i++) n_tot[i] = n1[i];

    // for (size_t k = 0; k < m1; k++)
    // {
    //     if (k == (m1 - 1))
    //     {
    //         n_tot[k] = int(alpha_tilde1 / stride);
    //     }
    //     else
    //     {
    //         n_tot[k] = int((alpha_tilde1 % n_xx1[k]) / stride);
    //         stride *= n_xx1[k];
    //     }
    // }

    for (int alpha2 = 0; alpha2 < dxx2_mult; alpha2++)
    {
        // stride = 1;
        // for (size_t k = 0; k < m2; k++)
        // {
        //     if (k == (m2 - 1))
        //     {
        //         n_tot[m1 + k] = int(i / stride);
        //     }
        //     else
        //     {
        //         n_tot[m1 + k] = int((i % n_xx2[k]) / stride);
        //         stride *= n_xx2[k];
        //     }
        // }
        n2 = index_to_state(alpha2, n_xx2);
        for (size_t i = 0; i < n2.size(); i++) n_tot[m1 + i] = n1[i];
        a_vec[alpha2] = mysystem.reactions[mu]->propensity(n_tot);
    }

    // Set up w_x2
    double d_xx2 = 1;
    for (auto & ele : h_xx2) d_xx2 *= ele;
    for (size_t i = 0; i < a_vec.size(); i++)
    {
        w_x2(i) = a_vec[i] * d_xx2;
    }

    // Calculate coefficient c2
    coeff(xx2_shift, lr_sol.V, w_x2, c2, blas);


    // Calculate coefficient d2
    coeff(lr_sol.V, lr_sol.V, w_x2, d2, blas);




    return 0;
}