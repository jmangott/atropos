#include "s_step_functions.hpp"

using std::vector;

void CalculateCoefficientsB(multi_array<double, 3> &b_coeff_vec_shift, multi_array<double, 3> &b_coeff_vec, const lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, Index mu, vector<Index> sigma2)
{
    Index d = grid.m1 + grid.m2;
    multi_array<double, 2> b_coeff({grid.r, grid.r});
    multi_array<double, 2> b_coeff_shift({grid.r, grid.r});
    multi_array<double, 1> a_vec({grid.dx2});
    multi_array<double, 1> state1({grid.m1});
    multi_array<double, 1> state2({grid.m2});
    vector<double> state_tot(d, 0.0);
    multi_array<Index, 1> vec_index1({grid.m1});
    multi_array<Index, 1> vec_index2({grid.m2});
    multi_array<double, 1> w_x2({grid.dx2});
    multi_array<double, 2> xx2_shift(lr_sol.V.shape());

    // Calculate the shifted X2
    ShiftMultiArrayRows(xx2_shift, lr_sol.V, -sigma2[mu]);

    for (Index alpha1 = 0; alpha1 < grid.dx1; alpha1++)
    {
        vec_index1 = CombIndexToVecIndex(alpha1, grid.n1);
        w_x2 = CalculateWeightX(2, vec_index1, reaction_system, grid, mu);

        coeff(xx2_shift, lr_sol.V, w_x2, b_coeff_shift, blas);
        coeff(lr_sol.V, lr_sol.V, w_x2, b_coeff, blas);

        for (Index i = 0; i < grid.r; i++)
        {
            for (Index j = 0; j < grid.r; j++)
            {
                b_coeff_vec_shift(alpha1, i, j) = b_coeff_shift(i, j);
                b_coeff_vec(alpha1, i, j) = b_coeff(i, j);
            }
        }
    }
}


void CalculateCoefficientsS(multi_array<double, 4> &e_coeff_tot, multi_array<double, 4> &f_coeff_tot, const multi_array<double, 3> &b_coeff_vec_shift, const multi_array<double, 3> &b_coeff_vec, const lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, Index mu, vector<Index> sigma1)
{
    multi_array<double, 2> e_coeff({grid.r, grid.r});
    multi_array<double, 2> f_coeff({grid.r, grid.r});
    multi_array<double, 1> w_x1({grid.dx1});
    multi_array<double, 1> w_x1_shift({grid.dx1});
    multi_array<double, 2> xx1_shift(lr_sol.X.shape());

    // Calculate the shifted X1
    ShiftMultiArrayRows(xx1_shift, lr_sol.X, -sigma1[mu]);

    // For integration weight
    double h_xx1_mult = 1;
    for (Index i = 0; i < grid.m1; i++)
        h_xx1_mult *= grid.h1(i);

    for (Index j = 0; j < grid.r; j++)
    {
        for (Index l = 0; l < grid.r; l++)
        {
            // Calculate integration weights
            for (Index alpha1 = 0; alpha1 < grid.dx1; alpha1++)
            {
                w_x1_shift(alpha1) = b_coeff_vec_shift(alpha1, j, l) * h_xx1_mult;
                w_x1(alpha1) = b_coeff_vec(alpha1, j, l) * h_xx1_mult;
            }

            coeff(xx1_shift, lr_sol.X, w_x1_shift, e_coeff, blas);
            coeff(lr_sol.X, lr_sol.X, w_x1, f_coeff, blas);

            for (Index i = 0; i < grid.r; i++)
            {
                for (Index k = 0; k < grid.r; k++)
                {
                    e_coeff_tot(i, j, k, l) = e_coeff(i, k);
                    f_coeff_tot(i, j, k, l) = f_coeff(i, k);
                }
            }
        }
    }
}

// Perform S-Step with time step size `tau`
void PerformSStep(vector<Index> sigma1, std::vector<Index> sigma2, lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, double tau)
{
    multi_array<double, 3> b_coeff_vec_shift({grid.dx1, grid.r, grid.r});
    multi_array<double, 3> b_coeff_vec({grid.dx1, grid.r, grid.r});
    multi_array<double, 4> e_coeff({grid.r, grid.r, grid.r, grid.r});
    multi_array<double, 4> f_coeff({grid.r, grid.r, grid.r, grid.r});
    multi_array<double, 2> s_dot(lr_sol.S.shape());
    set_zero(s_dot);

    for (std::vector<myreact *>::size_type mu = 0; mu < reaction_system.mu(); mu++)
    {
        CalculateCoefficientsB(b_coeff_vec_shift, b_coeff_vec, lr_sol, blas, reaction_system, grid, mu, sigma2);
        CalculateCoefficientsS(e_coeff, f_coeff, b_coeff_vec_shift, b_coeff_vec, lr_sol, blas, reaction_system, grid, mu, sigma1);

        for (Index i = 0; i < grid.r; i++)
        {
            for (Index j = 0; j < grid.r; j++)
            {
                for (Index k = 0; k < grid.r; k++)
                {
                    for (Index l = 0; l < grid.r; l++)
                    {
                        s_dot(i, j) += tau * lr_sol.S(k, l) * e_coeff(i, j, k, l);
                        s_dot(i, j) -= tau * lr_sol.S(k, l) * f_coeff(i, j, k, l);
                    }
                }
            }
        }
    }
    lr_sol.S -= s_dot;
}