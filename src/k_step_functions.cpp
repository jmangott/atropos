#include "k_step_functions.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;


multi_array<double, 1> CalculateWeightX2(multi_array<Index, 1> n_xx1, multi_array<Index, 1> n_xx2, multi_array<double, 1> lim_xx1, multi_array<double, 1> lim_xx2, Index dxx2_mult, multi_array<double, 1> h_xx2, multi_array<Index, 1> vec_index1, mysys reaction_system, Index mu)
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
        a_vec(alpha2) = reaction_system.reactions[mu]->propensity(state_tot);
    }

    // Calculate integration weight
    double h_xx2_mult = 1;
    for (auto const &ele : h_xx2)
        h_xx2_mult *= ele;
    for (Index i = 0; i < dxx2_mult; i++)
        w_x2(i) = a_vec(i) * h_xx2_mult;

    return w_x2;
}


void CalculateCoefficientsX2(multi_array<double, 2> &c2, multi_array<double, 2> &d2, multi_array<Index, 1> n_xx1, multi_array<Index, 1> n_xx2, multi_array<double, 1> lim_xx1, multi_array<double, 1> lim_xx2, multi_array<double, 1> h_xx2, lr2<double> &lr_sol, blas_ops blas, Index shift, multi_array<Index, 1> vec_index1, mysys reaction_system, Index mu)
{
    Index dxx2_mult = lr_sol.V.shape()[0];
    multi_array<double, 1> w_x2({dxx2_mult});

    multi_array<double, 2> xx2_shift(lr_sol.V.shape());

    // TODO: check shift sign
    // Calculate the shifted X2 (the value -shift is due to *inverse* shift operation)
    // NOTE: shift^-1(V^T) != (shift^-1(V))^T
    // CASE 1: (shift^-1(V))^T
    ShiftMultiArrayRows(xx2_shift, lr_sol.V, -shift);

    // CASE 2: shift^-1(V^T)
    // NOTE: I think in reality this case makes no sense, as one would have to assume (instead of constant basis functions for population numbers on the edge of the considered box) that the basis functions f_i themselves are approximately equal for large and small indices i, but the order of the basis functions is arbitrary
    // multi_array<double, 2> xx2_trans(lr_sol.V.shape());
    // xx2_trans = lr_sol.V;
    // transpose_inplace(xx2_trans);
    // ShiftMultiArrayRows(xx2_shift, xx2_trans, -shift);
    // transpose_inplace(xx2_shift);

    w_x2 = CalculateWeightX2(n_xx1, n_xx2, lim_xx1, lim_xx2, dxx2_mult, h_xx2, vec_index1, reaction_system, mu);

    coeff(xx2_shift, lr_sol.V, w_x2, c2, blas);
    coeff(lr_sol.V, lr_sol.V, w_x2, d2, blas);
}


void PerformKStep(multi_array<Index, 1> n_xx1, multi_array<Index, 1> n_xx2, multi_array<double, 1> h_xx2, multi_array<double, 1> lim_xx1, multi_array<double, 1> lim_xx2, vector<Index> sigma1, vector<Index> sigma2, lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, double tau)
{
    Index dxx1_mult = lr_sol.X.shape()[0];
    Index dxx2_mult = lr_sol.V.shape()[0];
    Index r = lr_sol.V.shape()[1];
    Index m1 = n_xx1.shape()[0];

    // For coefficients
    multi_array<double, 2> c2({r, r});
    multi_array<double, 2> d2({r, r});

    multi_array<double, 1> w_x2({dxx2_mult});

    vector<Index> dep_vec, dep_vec1, dep_vec2;
    vector<Index> n_xx1_coeff;
    multi_array<Index, 1> n_xx1_reduced({m1});
    vector<Index> vec_index1_coeff;
    multi_array<Index, 1> vec_index1_reduced({m1});
    multi_array<Index, 1> vec_index1_zero({m1});
    vector<double> lim_xx1_coeff;
    multi_array<double, 1> lim_xx1_reduced({m1});
    Index dxx1_coeff_mult, dxx1_reduced_mult;
    dep_vec = reaction_system.reactions[0]->depends_on;
    Index alpha1;

    multi_array<double, 2> prod_KC2({dxx1_mult, r});
    multi_array<double, 2> prod_KC2_shift({dxx1_mult, r});
    multi_array<double, 2> prod_KD2({dxx1_mult, r});
    multi_array<double, 2> k_dot({dxx1_mult, r});

    for (vector<myreact *>::size_type mu = 0; mu < reaction_system.mu(); mu++)
    {
        dep_vec = reaction_system.reactions[mu]->depends_on;
        dep_vec1.clear();
        dep_vec2.clear();
        dxx1_coeff_mult = 1;
        dxx1_reduced_mult = 1;
        n_xx1_coeff.clear();
        n_xx1_reduced = n_xx1;
        lim_xx1_coeff.clear();
        lim_xx1_reduced = lim_xx1;

        for (auto &ele : vec_index1_zero)
            ele = 0;

        for (auto const &ele : dep_vec)
        {
            if (ele < m1)
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
            
            // Convert vec_index1_coeff to a vector with size m1
            for (vector<Index>::size_type j = 0; j < dep_vec1.size(); j++)
                vec_index1_zero(dep_vec1[j]) = vec_index1_coeff[j];
            
            CalculateCoefficientsX2(c2, d2, n_xx1, n_xx2, lim_xx1, lim_xx2, h_xx2, lr_sol, blas, sigma2[mu], vec_index1_zero, reaction_system, mu);

            // Loop through the remaining species in partition 1
            for (Index k = 0; k < dxx1_reduced_mult; k++)
            {
                vec_index1_reduced = CombIndexToVecIndex(k, n_xx1_reduced);
                
                // vec_index1_reduced contains now the real population number
                for (vector<Index>::size_type l = 0; l < dep_vec1.size(); l++)
                    vec_index1_reduced(dep_vec1[l]) = vec_index1_coeff[l];
                alpha1 = VecIndexToCombIndex(vec_index1_reduced, n_xx1);

                // Calculate matrix-vector multiplication of C2*K and D2*K
                for (Index j = 0; j < r; j++)
                {
                    prod_KC2(alpha1, j) = 0.0;
                    prod_KD2(alpha1, j) = 0.0;
                    for (Index l = 0; l < r; l++)
                    {
                        prod_KC2(alpha1, j) += tau * lr_sol.X(alpha1, l) * c2(l, j);
                        prod_KD2(alpha1, j) += tau * lr_sol.X(alpha1, l) * d2(l, j);
                    }
                }
            }
        }
        // Shift prod_c2K
        // TODO: check shift signature
        ShiftMultiArrayRows(prod_KC2_shift, prod_KC2, +sigma1[mu]); // -sigma1[mu]

        // Calculate k_dot = shift(C2 * K) - D2 * K
        prod_KC2_shift -= prod_KD2;
        k_dot += prod_KC2_shift;
    }
    lr_sol.X += k_dot;
}
