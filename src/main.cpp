#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
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


// Convert index vector to state vector
vector<double> IndexToState(vector<int> index_vec, vector<Index> interval, vector<double> limit)
{
    size_t dim = index_vec.size();
    vector<double> state_vec(dim, 0.0);
    for (size_t i = 0; i < dim; i++)
    {
        state_vec[i] = index_vec[i] * limit[i] / (interval[i] - 1.0);
    }
    return state_vec;
}


// Convert index vector associated with the population numbers to combined index
Index VecIndexToCombIndex(vector<double> state_vec, vector<Index> interval, vector<double> limit)
{
    Index combined_index = 0;
    Index state_index;

    for (size_t i = 0; i < interval.size(); i++)
    {
        state_index = int(state_vec[i] * (interval[i] - 1.0) / limit[i]);
        combined_index += state_index * pow(interval[i], i);
    }
    return combined_index;
}


// Convert combined index to index vector associated with the population numbers
vector<double> CombIndexToVecIndex(Index combined_index, vector<Index> interval, vector<double> limit)
{
    int stride = 1;
    size_t dim = interval.size();
    double state_index;
    vector<double> state_vec(dim, 0.0);
    for (size_t k = 0; k < dim; k++)
    {
        if (k == (dim - 1))
        {
            state_index = (int(combined_index / stride));
        }
        else
        {
            state_index = (int((combined_index % interval[k]) / stride));
            stride *= interval[k];
        }
        if ((std::abs(interval[k] - 1.0) < 1.0e-10) &&
            (std::abs(interval[k] - 1.0 - limit[k]) < 1.0e-10))
            state_vec[k] = 0;
        else
            state_vec[k] = state_index * limit[k] / (interval[k] - 1.0);
    }
    return state_vec;
}


// Calculate the integration weight for coefficients C2 and D2 (depending on `alpha1` and `mu`)
multi_array<double, 1> CalculateWeightX2(vector<Index> n_xx1, vector<Index> n_xx2, vector<double> lim_xx2, Index dxx2_mult, vector<double> h_xx2, vector<double> n1, int mu)
{
    size_t m1 = n_xx1.size();
    size_t m2 = n_xx2.size();
    size_t d = m1 + m2;
    vector<double> a_vec(dxx2_mult, 0);
    vector<double> n_tot(d, 0);
    vector<double> n2(m2, 0);
    multi_array<double, 1> w_x2({dxx2_mult});

    for (size_t i = 0; i < m1; i++)
    {
        n_tot[i] = n1[i];
    }

    // Calculate propensity vector
    for (int alpha2 = 0; alpha2 < dxx2_mult; alpha2++)
    {
        n2 = CombIndexToVecIndex(alpha2, n_xx2, lim_xx2);
        for (size_t i = 0; i < m2; i++)
            n_tot[m1 + i] = n2[i];
        a_vec[alpha2] = mysystem.reactions[mu]->propensity(n_tot);
    }

    // Calculate integration weight
    double h_xx2_mult = 1;
    for (auto &ele : h_xx2)
        h_xx2_mult *= ele;
    for (int i = 0; i < int(a_vec.size()); i++)
        w_x2(i) = a_vec[i] * h_xx2_mult;

    return w_x2;
}


// Calculate `output_array`, where rows of `input_array` are shifted by `shift`
void ShiftMultiArrayCols(multi_array<double, 2> &output_array, multi_array<double, 2> &input_array, int shift)
{
    if ((output_array.shape()[0] != input_array.shape()[0]) ||
        (output_array.shape()[1] != input_array.shape()[1]))
    {
        std::cerr << "ERROR: Dimensions of output_array and input_array must be the same!";
        std::abort();
    }

    int n_rows = output_array.shape()[0];
    int n_cols = output_array.shape()[1];

    for (int j = 0; j < n_cols; j++)
    {
        for (int i = 0; i < n_rows; i++)
        {
            if ((shift < 0) && (i < -shift))
            {
                output_array(i, j) = input_array(0, j);
            }
            else if ((shift > 0) && (i >= (n_rows - shift)))
            {
                output_array(i, j) = input_array(n_rows - 1, j);
            }
            else
            {
                output_array(i, j) = input_array(i + shift, j);
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


// Calculate coefficients C2 and D2 for all values of `dep_vec` for a given reaction `mu`
void CalculateCoefficientsX2(multi_array<double, 2> &c2, multi_array<double, 2> &d2, vector<Index> n_xx1, vector<Index> n_xx2, vector<double> lim_xx1, vector<double> h_xx2, lr2<double> &lr_sol, blas_ops blas, int shift, int mu, vector<double> n1)
{
    Index dxx2_mult = lr_sol.V.shape()[0];
    multi_array<double, 1> w_x2({dxx2_mult});

    // Calculate the shifted X2
    multi_array<double, 2> xx2_shift(lr_sol.V.shape());

    ShiftMultiArrayCols(xx2_shift, lr_sol.V, shift);

    w_x2 = CalculateWeightX2(n_xx1, n_xx2, lim_xx1, dxx2_mult, h_xx2, n1, mu);
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

    Index nsteps = 1000; // # time steps
    double tstar = 1.0; // final time
    double tau = tstar / nsteps; // time step size

    constexpr size_t kM1 = kD / 2;
    constexpr size_t kM2 = kD - kM1;

    // Number of grid points for each species in partition 1 and 2
    vector<Index> n_xx1(kM1, 51);
    vector<Index> n_xx2(kM2, 51);

    // Limits for the population number in partition 1 and 2
    vector<double> lim_xx1(kM1, 50);
    vector<double> lim_xx2(kM2 ,50);

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
    int sigma1_sum, sigma2_sum;
    
    // NOTE: when the partition requires a permutation of the original order of species,
    // then also the nu vectors and similar quantities have to be permuted

    for(auto it : mysystem.reactions) {
        sigma1_sum = 0;
        sigma2_sum = 0;
        for (size_t i = 0; i < kM1; i++)
            sigma1_sum += it->nu[i] * pow(n_xx1[i], i);
        for (size_t i = 0; i < kM2; i++)
            sigma2_sum += it->nu[i] * pow(n_xx2[i], i);
        sigma1.push_back(sigma1_sum);
        sigma2.push_back(sigma2_sum);
    }

    vector<size_t> dep_vec, dep_vec1, dep_vec2;
    vector<Index> n_xx1_coeff, n_xx1_reduced;
    vector<double> n1_coeff, n1_reduced, n1_zero(kM1, 0);
    vector<double> lim_xx1_coeff, lim_xx1_reduced;
    Index dxx1_coeff_mult, dxx1_reduced_mult;
    dep_vec = mysystem.reactions[0]->depends_on;
    int alpha1;

    multi_array<double, 2> prod_c2K({dxx1_mult, kR});
    multi_array<double, 2> prod_c2K_shift({dxx1_mult, kR});
    multi_array<double, 2> prod_d2K({dxx1_mult, kR});

    for (size_t mu = 0; mu < mysystem.mu(); mu++)
    {
        dep_vec = mysystem.reactions[mu]->depends_on;
        dep_vec1.clear();
        dep_vec2.clear();
        dxx1_coeff_mult = 1;
        dxx1_reduced_mult = 1;
        n_xx1_coeff.clear();
        n_xx1_reduced = n_xx1;
        lim_xx1_coeff.clear();
        lim_xx1_reduced = lim_xx1;
        std::fill(n1_zero.begin(), n1_zero.end(), 0.0);
        for (auto &ele : dep_vec)
        {
            if (ele < kM1)
            {
                dep_vec1.push_back(ele);
                dxx1_coeff_mult *= n_xx1[ele];
                n_xx1_coeff.push_back(n_xx1[ele]);
                n_xx1_reduced[ele] = 1;
                lim_xx1_coeff.push_back(lim_xx1[ele]);
                lim_xx1_reduced[ele] = 0.0;
            }
            else
            {
                dep_vec2.push_back(ele);
            }
        }

        for (auto &ele : n_xx1_reduced)
        {
            dxx1_reduced_mult *= ele;
        }

        // Loop through all species in partition 1 on which the propensity for reaction mu depends
        for (int i = 0; i < dxx1_coeff_mult; i++)
        {
            n1_coeff = CombIndexToVecIndex(i, n_xx1_coeff, lim_xx1_coeff);
            
            // Convert n1_coeff to a vector with size kM1
            for (size_t j = 0; j < dep_vec1.size(); j++)
                n1_zero[dep_vec1[j]] = n1_coeff[j];

            // Loop through the remaining species in partition 1
            for (int k = 0; k < dxx1_reduced_mult; k++)
            {
                n1_reduced = CombIndexToVecIndex(k, n_xx1_reduced, lim_xx1_reduced);
                
                // n1_reduced contains now the real population number
                for (size_t l = 0; l < dep_vec1.size(); l++)
                    n1_reduced[dep_vec1[l]] = n1_coeff[l];
                alpha1 = VecIndexToCombIndex(n1_reduced, n_xx1, lim_xx1);

                CalculateCoefficientsX2(c2, d2, n_xx1, n_xx2, lim_xx1, h_xx2, lr_sol, blas, k_xx2[0] * sigma2[mu], mu, n1_zero);

                // Calculate matrix-vector multiplication of C2*K and D2*K
                set_zero(prod_c2K);
                set_zero(prod_d2K);

                for (int j = 0; j < kR; j++)
                {
                    prod_c2K(alpha1, j) = 0.0;
                    for (int l = 0; l < kR; l++)
                    {
                        prod_c2K(alpha1, j) += tau * c2(j, l) * lr_sol.X(alpha1, l);
                        prod_d2K(alpha1, j) += tau * d2(j, l) * lr_sol.X(alpha1, l);
                    }
                }
                // Shift prod_c2K
                ShiftMultiArrayCols(prod_c2K_shift, prod_c2K, -k_xx1[0] * sigma1[mu]);

                prod_c2K_shift -= prod_d2K;
                lr_sol.X += prod_c2K_shift;
            }
        }
    }

    return 0;
}