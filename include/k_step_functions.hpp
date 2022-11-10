#ifndef K_STEP_FUNCTIONS_HPP
#define K_STEP_FUNCTIONS_HPP

#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"
#include "k_step_functions.hpp"
#include "reaction_class.hpp"


// Calculate the integration weight for coefficients C2 and D2 (depending on `alpha1` and `mu`)
template <Index m1, Index m2>
multi_array<double, 1> CalculateWeightX2(multi_array<Index, 1> vec_index1, mysys reaction_system, grid_info<m1, m2> grid, Index mu)
{
    Index d = m1 + m2;
    multi_array<double, 1> a_vec({grid.dx2});
    multi_array<double, 1> state1({m1});
    multi_array<double, 1> state2({m2});
    std::vector<double> state_tot(d, 0.0);
    multi_array<Index, 1> vec_index2({m2});
    multi_array<double, 1> w_x2({grid.dx2});

    state1 = VecIndexToState(vec_index1, grid.n1, grid.lim1);

    for (Index i = 0; i < m1; i++)
    {
        state_tot[i] = state1(i);
    }

    // Calculate propensity vector
    for (Index alpha2 = 0; alpha2 < grid.dx2; alpha2++)
    {
        vec_index2 = CombIndexToVecIndex(alpha2, grid.n2);
        state2 = VecIndexToState(vec_index2, grid.n2, grid.lim2);

        // `state_tot` is the concatenation of `state1` and `state2`
        for (Index i = 0; i < m2; i++)
            state_tot[m1 + i] = state2(i);
        a_vec(alpha2) = reaction_system.reactions[mu]->propensity(state_tot);
    }

    // Calculate integration weight
    double h_xx2_mult = 1;
    for (Index i = 0; i < m2; i++)
        h_xx2_mult *= grid.h2(i);
    for (Index i = 0; i < grid.dx2; i++)
        w_x2(i) = a_vec(i) * h_xx2_mult;

    return w_x2;
}


// Calculate coefficients C2 and D2 for all values of `dep_vec` for a given reaction `mu`
template <Index m1, Index m2>
void CalculateCoefficientsX2(multi_array<double, 2> &c2, multi_array<double, 2> &d2, lr2<double> &lr_sol, blas_ops blas, Index shift, multi_array<Index, 1> vec_index1, mysys reaction_system, grid_info<m1, m2> grid, Index mu)
{
    multi_array<double, 1> w_x2({grid.dx2});
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

    w_x2 = CalculateWeightX2<m1, m2>(vec_index1, reaction_system, grid, mu);

    coeff(xx2_shift, lr_sol.V, w_x2, c2, blas);
    coeff(lr_sol.V, lr_sol.V, w_x2, d2, blas);
}


// Perform K-Step with time step size `tau`
template <Index m1, Index m2>
void PerformKStep(std::vector<Index> sigma1, vector<Index> sigma2, lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info<m1, m2> grid, double tau)
{
    // For coefficients
    multi_array<double, 2> c2({grid.r, grid.r});
    multi_array<double, 2> d2({grid.r, grid.r});

    multi_array<double, 1> w_x2({grid.dx2});

    std::vector<Index> dep_vec, dep_vec1, dep_vec2;
    std::vector<Index> n_xx1_coeff;
    multi_array<Index, 1> n_xx1_reduced({m1});
    std::vector<Index> vec_index1_coeff;
    multi_array<Index, 1> vec_index1_reduced({m1});
    multi_array<Index, 1> vec_index1_zero({m1});
    std::vector<double> lim_xx1_coeff;
    multi_array<double, 1> lim_xx1_reduced({m1});
    Index dxx1_coeff_mult, dxx1_reduced_mult;
    dep_vec = reaction_system.reactions[0]->depends_on;
    Index alpha1;

    multi_array<double, 2> prod_KC2({grid.dx1, grid.r});
    multi_array<double, 2> prod_KC2_shift({grid.dx1, grid.r});
    multi_array<double, 2> prod_KD2({grid.dx2, grid.r});
    multi_array<double, 2> k_dot({grid.dx1, grid.r});
    set_zero(k_dot);

    for (std::vector<myreact *>::size_type mu = 0; mu < reaction_system.mu(); mu++)
    {
        dep_vec = reaction_system.reactions[mu]->depends_on;
        dep_vec1.clear();
        dep_vec2.clear();
        dxx1_coeff_mult = 1;
        dxx1_reduced_mult = 1;
        n_xx1_coeff.clear();
        n_xx1_reduced = grid.n1;
        lim_xx1_coeff.clear();
        lim_xx1_reduced = grid.lim1;
        set_zero(prod_KC2);
        set_zero(prod_KD2);

        for (Index i = 0; i < m1; i++)
            vec_index1_zero(i) = 0;

        for (auto const &ele : dep_vec)
        {
            if (ele < m1)
            {
                dep_vec1.push_back(ele);
                dxx1_coeff_mult *= grid.n1(ele);
                n_xx1_coeff.push_back(grid.n1(ele));
                n_xx1_reduced(ele) = 1;
                lim_xx1_coeff.push_back(grid.lim1(ele));
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
            for (std::vector<Index>::size_type j = 0; j < dep_vec1.size(); j++)
                vec_index1_zero(dep_vec1[j]) = vec_index1_coeff[j];

            CalculateCoefficientsX2<m1, m2>(c2, d2, lr_sol, blas, sigma2[mu], vec_index1_zero, reaction_system, grid, mu);

            // Loop through the remaining species in partition 1
            for (Index k = 0; k < dxx1_reduced_mult; k++)
            {
                vec_index1_reduced = CombIndexToVecIndex(k, n_xx1_reduced);

                // vec_index1_reduced contains now the real population number
                for (std::vector<Index>::size_type l = 0; l < dep_vec1.size(); l++)
                    vec_index1_reduced(dep_vec1[l]) = vec_index1_coeff[l];
                alpha1 = VecIndexToCombIndex(vec_index1_reduced, grid.n1);

                // Calculate matrix-vector multiplication of C2*K and D2*K
                for (Index j = 0; j < grid.r; j++)
                {
                    for (Index l = 0; l < grid.r; l++)
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

#endif