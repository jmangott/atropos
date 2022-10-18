#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "index_functions.hpp"
#include "io_functions.hpp"
#include "reactions.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;

// Calculate the integration weight for coefficients C2 and D2 (depending on `alpha1` and `mu`)
multi_array<double, 1> CalculateWeightX2(multi_array<Index, 1> n_xx1, multi_array<Index, 1> n_xx2, multi_array<double, 1> lim_xx1, multi_array<double, 1> lim_xx2, Index dxx2_mult, multi_array<double, 1> h_xx2, multi_array<Index, 1> vec_index1, Index mu)
{
    Index m1 = n_xx1.shape()[0];
    Index m2 = n_xx2.shape()[0];
    Index d = m1 + m2;
    multi_array<double, 1> a_vec({dxx2_mult});
    multi_array<double, 1> state1({m1});
    multi_array<double, 1> state2({m2});
    vector<double> state_tot(d, 0.0);
    multi_array<Index, 1> vec_index2({m2});
    multi_array<double, 1> w_x2({dxx2_mult});

    state1 = VecIndexToState(vec_index1, n_xx1, lim_xx1);

    for (Index i = 0; i < m1; i++)
    {
        state_tot[i] = state1(i);
    }

    // Calculate propensity vector
    for (Index alpha2 = 0; alpha2 < dxx2_mult; alpha2++)
    {
        vec_index2 = CombIndexToVecIndex(alpha2, n_xx2);
        state2 = VecIndexToState(vec_index2, n_xx2, lim_xx2);

        // `state_tot` is the concatenation of `state1` and `state2`
        for (Index i = 0; i < m2; i++)
            state_tot[m1 + i] = state2(i);
        a_vec(alpha2) = mysystem.reactions[mu]->propensity(state_tot);
    }

    // Calculate integration weight
    double h_xx2_mult = 1;
    for (auto const &ele : h_xx2)
        h_xx2_mult *= ele;
    for (Index i = 0; i < dxx2_mult; i++)
        w_x2(i) = a_vec(i) * h_xx2_mult;

    return w_x2;
}


// Calculate coefficients C2 and D2 for all values of `dep_vec` for a given reaction `mu`
void CalculateCoefficientsX2(multi_array<double, 2> &c2, multi_array<double, 2> &d2, multi_array<Index, 1> n_xx1, multi_array<Index, 1> n_xx2, multi_array<double, 1> lim_xx1, multi_array<double, 1> lim_xx2, multi_array<double, 1> h_xx2, lr2<double> &lr_sol, blas_ops blas, Index shift, multi_array<Index, 1> vec_index1, Index mu)
{
    Index dxx2_mult = lr_sol.V.shape()[0];
    multi_array<double, 1> w_x2({dxx2_mult});

    // Calculate the shifted X2
    multi_array<double, 2> xx2_shift(lr_sol.V.shape());

    ShiftMultiArrayCols(xx2_shift, lr_sol.V, shift);

    w_x2 = CalculateWeightX2(n_xx1, n_xx2, lim_xx1, lim_xx2, dxx2_mult, h_xx2, vec_index1, mu);
    coeff(xx2_shift, lr_sol.V, w_x2, c2, blas);
    coeff(lr_sol.V, lr_sol.V, w_x2, d2, blas);
}


int main()
{
    /////////////////////////////////////////////
    ////// DECLARATION AND INITIALIZATION ///////
    /////////////////////////////////////////////

    constexpr Index kR = 5; // rank
    constexpr Index kD = 2; // # species

    constexpr Index kN = 51;
    constexpr Index kK = 1;
    constexpr Index kZero = 0;
 
    nn = {"S1", "S2"};

    Index nsteps = 1000; // # time steps
    double tstar = 1.0; // final time
    double tau = tstar / nsteps; // time step size

    constexpr Index kM1 = kD / 2;
    constexpr Index kM2 = kD - kM1;

    // Number of grid points for each species in partition 1 and 2
    multi_array<Index, 1> n_xx1({kM1});
    multi_array<Index, 1> n_xx2({kM2});

    for (auto &ele : n_xx1)
        ele = kN;
    for (auto &ele : n_xx2)
        ele = kN;

    // Limits for the population number in partition 1 and 2
    multi_array<double, 1> lim_xx1({kM1});
    multi_array<double, 1> lim_xx2({kM2});

    // Initial datum generation
    multi_array<double, 1> h_xx1({kM1});
    multi_array<double, 1> h_xx2({kM2});
    multi_array<Index, 1> k_xx1({kM1});
    multi_array<Index, 1> k_xx2({kM2});

    for (auto &ele : k_xx1)
        ele = kK;
    for (auto &ele : k_xx2)
        ele = kK;

    for (Index ii = 0; ii < kM1; ii++)
    {
        if ((n_xx1(ii) - 1) % k_xx1(ii) != 0)
        {
            std::cerr << "ERROR: (n_xx1 - 1) must be a multiple of k_xx1!";
            std::abort();
        }
        else
        {
            lim_xx1(ii) = (n_xx1(ii) - 1) / k_xx1(ii);
            h_xx1(ii) = 1.0 / k_xx1(ii);
        }
    }

    for (Index ii = 0; ii < kM2; ii++)
    {
        if ((n_xx2(ii) - 1) % k_xx2(ii) != 0)
        {
            std::cerr << "ERROR: (n_xx2 - 1) must be a multiple of k_xx2!";
            std::abort();
        }
        else
        {
            lim_xx2(ii) = (n_xx2(ii) - 1) / k_xx2(ii);
            h_xx2(ii) = 1.0 / k_xx2(ii);
        }
    }

    Index dxx1_mult = 1;
    Index dxx2_mult = 1;
    for (auto const ele : n_xx1)
        dxx1_mult *= ele;
    for (auto const ele : n_xx2)
        dxx2_mult *= ele;

    vector<const double *> x1, x2;

    // For coefficients
    multi_array<double, 2> c2({kR, kR});
    multi_array<double, 2> d2({kR, kR});

    multi_array<double, 1> w_x2({dxx2_mult});

    // Temporary objects for multiplication
    multi_array<double, 2> tmp_x({dxx1_mult, kR});

    // Set up S, X1 and X2h for t = 0
    // multi_array<double, 2> ss({kR, kR});
    multi_array<double, 2> xx1({dxx1_mult, kR});
    multi_array<double, 2> xx2({dxx2_mult, kR});

    // Read in S, X1 and X2
    // ReadInMultiArray(ss, "../input/s.csv");
    ReadInMultiArray(xx1, "../input/u.csv");
    ReadInMultiArray(xx2, "../input/vh.csv");
    
    // Point to the beginning of every column of X1 and X2
    double *it1 = xx1.begin();
    double *it2 = xx2.begin();
    for (int i = 0; i < kR; i++)
    {
        x1.push_back(it1);
        x2.push_back(it2);
        it1 += kM1 * n_xx1(0);
        it2 += kM2 * n_xx2(0);
    }
    
    // Set up the low-rank structure and the inner products
    lr2<double> lr_sol(kR, {dxx1_mult, dxx2_mult});

    std::function<double(double *, double *)> ip_xx1 = inner_product_from_const_weight(h_xx1(0), dxx1_mult);
    std::function<double(double *, double *)> ip_xx2 = inner_product_from_const_weight(h_xx2(0), dxx2_mult);

    blas_ops blas;
    initialize(lr_sol, x1, x2, ip_xx1, ip_xx2, blas);
    // lr_sol.S = ss;

    // Calculate the shift amount for all reactions (this has to be done only once)
    vector<Index> sigma1, sigma2;
    CalculateShiftAmount(sigma1, sigma2, mysystem, n_xx1, n_xx2);


    /////////////////////////////////////////////
    ////////////////// K-STEP ///////////////////
    /////////////////////////////////////////////

    tmp_x = lr_sol.X;
    blas.matmul(tmp_x, lr_sol.S, lr_sol.X); // lr_sol.X contains now K

    vector<Index> dep_vec, dep_vec1, dep_vec2;
    vector<Index> n_xx1_coeff;
    multi_array<Index, 1> n_xx1_reduced({kM1});
    vector<Index> vec_index1_coeff;
    multi_array<Index, 1> vec_index1_reduced({kM1});
    multi_array<Index, 1> vec_index1_zero({kM1});
    vector<double> lim_xx1_coeff;
    multi_array<double, 1> lim_xx1_reduced({kM1});
    Index dxx1_coeff_mult, dxx1_reduced_mult;
    dep_vec = mysystem.reactions[0]->depends_on;
    Index alpha1;

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

        for (auto &ele : vec_index1_zero)
            ele = kZero;

        for (auto const &ele : dep_vec)
        {
            if (ele < kM1)
            {
                dep_vec1.push_back(ele);
                dxx1_coeff_mult *= n_xx1(ele);
                n_xx1_coeff.push_back(n_xx1(ele));
                n_xx1_reduced(ele) = 1;
                lim_xx1_coeff.push_back(lim_xx1(ele));
                lim_xx1_reduced(ele) = 0.0;
            }
            else
            {
                dep_vec2.push_back(ele);
            }
        }

        for (auto const &ele : n_xx1_reduced)
        {
            dxx1_reduced_mult *= ele;
        }

        // Loop through all species in partition 1 on which the propensity for reaction mu depends
        for (Index i = 0; i < dxx1_coeff_mult; i++)
        {
            vec_index1_coeff = CombIndexToVecIndex(i, n_xx1_coeff);
            
            // Convert vec_index1_coeff to a vector with size kM1
            for (vector<Index>::size_type j = 0; j < dep_vec1.size(); j++)
                vec_index1_zero(dep_vec1[j]) = vec_index1_coeff[j];

            // Loop through the remaining species in partition 1
            for (Index k = 0; k < dxx1_reduced_mult; k++)
            {
                vec_index1_reduced = CombIndexToVecIndex(k, n_xx1_reduced);
                
                // vec_index1_reduced contains now the real population number
                for (vector<Index>::size_type l = 0; l < dep_vec1.size(); l++)
                    vec_index1_reduced(dep_vec1[l]) = vec_index1_coeff[l];
                alpha1 = VecIndexToCombIndex(vec_index1_reduced, n_xx1);

                CalculateCoefficientsX2(c2, d2, n_xx1, n_xx2, lim_xx1, lim_xx2, h_xx2, lr_sol, blas, k_xx2(0) * sigma2[mu], vec_index1_zero, mu);

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
                ShiftMultiArrayCols(prod_c2K_shift, prod_c2K, -k_xx1(0) * sigma1[mu]);

                prod_c2K_shift -= prod_d2K;
                lr_sol.X += prod_c2K_shift;
            }
        }
    }

    return 0;
}