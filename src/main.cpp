#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
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
using std::map;
using std::stringstream;
using std::string;
using std::vector;

// TODO: Improve speed by using a lookup table (?)

// Convert combined index to vector containing the population numbers
vector<double> IndexToState(Index idx, vector<Index> interval)
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


// Calculate the integration weight for coefficients C2 and D2 (depending on `alpha1` and `mu`)
multi_array<double, 1> CalculateWeightX2(int d, vector<Index> n_xx1, vector<Index> n_xx2, Index dxx2_mult, vector<double> h_xx2, vector<double> n1, int mu)
{
    size_t m1 = n_xx1.size();
    size_t m2 = n_xx2.size();
    vector<double> a_vec(dxx2_mult, 0);
    vector<double> n_tot(d, 0);
    vector<double> n2(m2, 0);
    multi_array<double, 1> w_x2({dxx2_mult});

    for (int alpha2 = 0; alpha2 < dxx2_mult; alpha2++)
    {
        n2 = IndexToState(alpha2, n_xx2);
        for (size_t i = 0; i < n2.size(); i++)
            n_tot[m1 + i] = n1[i];
        a_vec[alpha2] = mysystem.reactions[mu]->propensity(n_tot);
    }

    double h_xx2_mult = 1;
    for (auto &ele : h_xx2)
        h_xx2_mult *= ele;
    for (int i = 0; i < int(a_vec.size()); i++)
        w_x2(i) = a_vec[i] * h_xx2_mult;

    return w_x2;
}


// Calculate array_shift, where rows of input_pointer are shifted by index_shift
void ShiftMultiArrayCols(multi_array<double, 2> &array_shift, vector<int> k_vec, vector<const double *> &input_pointer, int index_shift)
{
    int n_rows = array_shift.shape()[0];
    int n_cols = array_shift.shape()[1];
    vector<const double *> pointer_shift = input_pointer;
    for (int i = 0; i < n_cols; i++)
        pointer_shift[i] += int(k_vec[0] * index_shift);

    for (int j = 0; j < n_cols; j++)
    {
        for (int i = 0; i < n_rows; i++)
        {
            if ((k_vec[0] * index_shift < 0) &&
                (i < -k_vec[0] * index_shift))
            {
                array_shift(i, j) = input_pointer[j][0];
            }
            else if ((k_vec[0] * index_shift > 0) &&
                     (i >= (n_rows - k_vec[0] * index_shift)))
            {
                array_shift(i, j) = input_pointer[j][n_rows - 1];
            }
            else
            {
                array_shift(i, j) = pointer_shift[j][i];
            }
        }
    }
}


// TODO: Use .netcdf instead of .csv

// Read in a .csv file and store it as a multi_array<double, 2>
void ReadInMultiArray(multi_array<double, 2> &output_array, string filename)
{
    int ii = 0;
    int jj = 0;
    string line;
    ifstream input_file(filename);
    stringstream ss_line;
    string element;

    while (getline(input_file, line))
    {
        ss_line.str(line);
        while (getline(ss_line, element, ','))
        {
            output_array(ii, jj) = stod(element);
            jj++;
        }
        ii++;
        jj = 0;
        ss_line.clear();
    }
    input_file.close();
}


// Calculate coefficients C2 and D2 for all values of`dep_vec` for a given reaction `mu`
void CalculateCoefficientsX2(multi_array<double, 2> &c2, multi_array<double, 2> &d2, int r, int d, vector<Index> n_xx1, vector<Index> n_xx2, vector<Index> n_xx1_coeff, int dxx2_mult, vector<double> h_xx2, vector<int> k_xx2, vector<const double *> &x2, lr2<double> lr_sol, blas_ops blas, int index_shift, int mu, int alpha_tilde1, vector<size_t> dep_vec)
{
    multi_array<double, 1> w_x2({dxx2_mult});
    vector<double> n1(n_xx1.size(), 0);
    vector<double> n1_coeff;

    // Calculate the shifted X2
    multi_array<double, 2> xx2_shift({dxx2_mult, r});
    ShiftMultiArrayCols(xx2_shift, k_xx2, x2, index_shift);

    n1_coeff = IndexToState(alpha_tilde1, n_xx1_coeff);

    // Convert n1_coeff to a vector with size kM1
    for (size_t j = 0; j < dep_vec.size(); j++)
        n1[dep_vec[j]] = n1_coeff[j];

    w_x2 = CalculateWeightX2(d, n_xx1, n_xx2, dxx2_mult, h_xx2, n1, mu);
    coeff(xx2_shift, lr_sol.V, w_x2, c2, blas);
    coeff(lr_sol.V, lr_sol.V, w_x2, d2, blas);
}


int main()
{
    /////////////////////////////////////////////
    ////// DECLARATION AND INITIALIZATION ///////
    /////////////////////////////////////////////

    constexpr int kR = 5; // rank
    constexpr int kD = 2; // # species
 
    std::fill_n(std::back_inserter(xx), kD, 1.0);
    nn = {"A1", "B2"};

    // Index nsteps = 1000; // # time steps
    // double tstar = 1.0; // final time
    // double tau = tstar / nsteps; // time step size

    constexpr size_t kM1 = kD / 2;
    constexpr size_t kM2 = kD - kM1;

    // Number of grid points for each species in partition 1 and 2
    vector<Index> n_xx1(kM1, 51);
    vector<Index> n_xx2(kM2, 51);

    // Limits for the population number in partition 1 and 2
    vector<int> lim_xx1(kM1, 50);
    vector<int> lim_xx2(kM2 ,50);

    // Initial datum generation
    vector<double> h_xx1(kM1, 0);
    vector<double> h_xx2(kM2, 0);
    vector<int> k_xx1(kM1, 0);
    vector<int> k_xx2(kM2, 0);
    for (size_t ii = 0; ii < kM1; ii++)
    {
        h_xx1[ii] = double(lim_xx1[ii]) / (n_xx1[ii] - 1);
        k_xx1[ii] = int((n_xx1[ii] - 1) / lim_xx1[ii]);
    }

    for (size_t ii = 0; ii < kM2; ii++)
    {
        h_xx2[ii] = lim_xx2[ii] / (n_xx2[ii] - 1);
        k_xx2[ii] = int((n_xx2[ii] - 1) / lim_xx2[ii]);
    }

    Index dxx1_mult = 1;
    Index dxx2_mult = 1;
    for (auto ele : n_xx1)
        dxx1_mult *= ele;
    for (auto ele : n_xx2)
        dxx2_mult *= ele;

    vector<const double *> x1, x2;

    // For coefficients
    multi_array<double, 2> c2({kR, kR});
    multi_array<double, 2> d2({kR, kR});

    multi_array<double, 1> w_x2({dxx2_mult});

    // Temporary objects for multiplication
    multi_array<double, 2> tmp_x({dxx1_mult, kR});

    // Set up S, X1 and X2h for t = 0
    multi_array<double, 2> ss({kR, kR});
    multi_array<double, 2> xx1({dxx1_mult, kR});
    multi_array<double, 2> xx2({dxx2_mult, kR});

    // Read in S, X1 and X2
    ReadInMultiArray(ss, "../input/s.csv");
    ReadInMultiArray(xx1, "../input/u.csv");
    ReadInMultiArray(xx2, "../input/vh.csv");
    
    // Point to the beginning of every column of X1 and X2
    double *it1 = xx1.begin();
    double *it2 = xx2.begin();
    for (size_t i = 0; i < kR; i++)
    {
        x1.push_back(it1);
        x2.push_back(it2);
        it1 += kM1 * n_xx1[0];
        it2 += kM2 * n_xx2[0];
    }
    
    // Set up the low-rank structure and the inner products
    lr2<double> lr_sol(kR, {dxx1_mult, dxx2_mult});

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
    
    vector<int> sigma1, sigma2;
    int sigma2_sum;
    
    // NOTE: when the partition requires a permutation of the original order of species,
    // then also the nu vectors and similar quantities have to be permuted

    for(auto it : mysystem.reactions) {
        sigma2_sum = 0;
        for (size_t i = 0; i < kM2; i++)
            sigma2_sum =+ it->nu[i] * pow(n_xx2[i], i);
        sigma2.push_back(sigma2_sum);
    }

    // Calculate C2 and D2 for given mu and alpha_tilde1
    int mu = 2;
    int alpha_tilde1 = 30;

    vector<size_t> dep_vec, dep_vec1, dep_vec2;
    dep_vec = mysystem.reactions[0]->depends_on;
    for (auto & ele : dep_vec)
        (ele < kM1) ? dep_vec1.push_back(ele) : dep_vec2.push_back(ele);

    int dxx1_coeff_mult = 1;
    vector<Index> n_xx1_coeff;

    for (auto ele : dep_vec1)
    {
        dxx1_coeff_mult *= ele;
        n_xx1_coeff.push_back(n_xx1[ele]);
    }

    CalculateCoefficientsX2(c2, d2, kR, kD, n_xx1, n_xx2, n_xx1_coeff, dxx2_mult, h_xx2, k_xx2, x2, lr_sol, blas, sigma2[mu], mu, alpha_tilde1, dep_vec1);

    for (int i = 0; i < kR; i++)
    {
        for (int j = 0; j < kR; j++)
        {
            cout << c2(i, j) << " ";
        }
        cout << endl;
    }

    return 0;
}