#include "s_step_functions.hpp"

using std::vector;

void CalculateCoefficientsB(multi_array<double, 3> &b_coeff_vec_shift, multi_array<double, 3> &b_coeff_vec, const lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, partition_info<1> partition1, partition_info<2> partition2, Index mu, vector<Index> sigma2, const vector<multi_array<double, 2>> &w_x_dep)
{
    multi_array<double, 2> b_coeff({grid.r, grid.r});
    multi_array<double, 2> b_coeff_shift({grid.r, grid.r});
    multi_array<double, 1> w_x2({grid.dx2});
    multi_array<double, 2> xx2_shift(lr_sol.V.shape());

    // Calculate the shifted X2
    ShiftMultiArrayRows(2, xx2_shift, lr_sol.V, -sigma2[mu], reaction_system.reactions[mu]->minus_nu, grid);

    Index alpha2_dep;

    for (Index alpha1_dep = 0; alpha1_dep < partition1.dx_dep(mu); alpha1_dep++)
    {
        get_time::start("bweightvec");
        for (Index alpha2 = 0; alpha2 < grid.dx2; alpha2++)
        {
            alpha2_dep = CombIndexToDepCombIndex(alpha2, partition2.n_dep[mu], grid.n2, partition2.dep_vec[mu]);
            w_x2(alpha2) = w_x_dep[mu](alpha1_dep, alpha2_dep) * grid.h2_mult;
        }
        get_time::stop("bweightvec");
        get_time::start("bcoeff_s");
        coeff(xx2_shift, lr_sol.V, w_x2, b_coeff_shift, blas);
        coeff(lr_sol.V, lr_sol.V, w_x2, b_coeff, blas);
        get_time::stop("bcoeff_s");

        for (Index i = 0; i < grid.r; i++)
        {
            for (Index j = 0; j < grid.r; j++)
            {
                b_coeff_vec_shift(alpha1_dep, i, j) = b_coeff_shift(i, j);
                b_coeff_vec(alpha1_dep, i, j) = b_coeff(i, j);
            }
        }
    }
}


void CalculateCoefficientsS(multi_array<double, 4> &e_coeff_tot, multi_array<double, 4> &f_coeff_tot, const multi_array<double, 3> &b_coeff_vec_shift, const multi_array<double, 3> &b_coeff_vec, const lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, partition_info<1> partition1, Index mu, vector<Index> sigma1)
{
    multi_array<double, 2> e_coeff({grid.r, grid.r});
    multi_array<double, 2> f_coeff({grid.r, grid.r});
    multi_array<double, 1> w_x1({grid.dx1});
    multi_array<double, 1> w_x1_shift({grid.dx1});
    multi_array<double, 2> xx1_shift(lr_sol.X.shape());

    Index alpha1_dep;

    // Calculate the shifted X1
    ShiftMultiArrayRows(1, xx1_shift, lr_sol.X, -sigma1[mu], reaction_system.reactions[mu]->minus_nu, grid);

    for (Index j = 0; j < grid.r; j++)
    {
        for (Index l = 0; l < grid.r; l++)
        {
            // Calculate integration weights
            get_time::start("sweightvec");
            for (Index alpha1 = 0; alpha1 < grid.dx1; alpha1++)
            {
                // get_time::start("sweightvec_index");
                alpha1_dep = CombIndexToDepCombIndex(alpha1, partition1.n_dep[mu], grid.n1, partition1.dep_vec[mu]);
                // get_time::stop("sweightvec_index");
                w_x1_shift(alpha1) = b_coeff_vec_shift(alpha1_dep, j, l) * grid.h1_mult;
                w_x1(alpha1) = b_coeff_vec(alpha1_dep, j, l) * grid.h1_mult;
            }
            get_time::stop("sweightvec");

            get_time::start("scoeff_s");
            coeff(xx1_shift, lr_sol.X, w_x1_shift, e_coeff, blas);
            coeff(lr_sol.X, lr_sol.X, w_x1, f_coeff, blas);
            get_time::stop("scoeff_s");

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
void PerformSStep(vector<Index> sigma1, std::vector<Index> sigma2, lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, partition_info<1> partition1, partition_info<2> partition2, const vector<multi_array<double, 2>> &w_x_dep, double tau)
{
    // multi_array<double, 3> b_coeff_vec_shift;
    // multi_array<double, 3> b_coeff_vec;
    multi_array<double, 4> e_coeff({grid.r, grid.r, grid.r, grid.r});
    multi_array<double, 4> f_coeff({grid.r, grid.r, grid.r, grid.r});
    multi_array<double, 2> s_dot(lr_sol.S.shape());
    set_zero(s_dot);

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        // TODO: Ask Lukas if this is good practice
        multi_array<double, 3> b_coeff_vec_shift({partition1.dx_dep(mu), grid.r, grid.r});
        multi_array<double, 3> b_coeff_vec({partition1.dx_dep(mu), grid.r, grid.r});
        
        // This causes a memory leak
        // b_coeff_vec_shift.resize({partition1.dx_dep(mu), grid.r, grid.r});
        // b_coeff_vec.resize({partition1.dx_dep(mu), grid.r, grid.r});

        get_time::start("bcoeff");
        CalculateCoefficientsB(b_coeff_vec_shift, b_coeff_vec, lr_sol, blas, reaction_system, grid, partition1, partition2, mu, sigma2, w_x_dep);
        get_time::stop("bcoeff");
        get_time::start("scoeff");
        CalculateCoefficientsS(e_coeff, f_coeff, b_coeff_vec_shift, b_coeff_vec, lr_sol, blas, reaction_system, grid, partition1, mu, sigma1);
        get_time::stop("scoeff");

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