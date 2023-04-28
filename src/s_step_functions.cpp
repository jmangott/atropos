#include "s_step_functions.hpp"

using std::vector;

void CalculateCoefficientsB(multi_array<double, 3> &b_coeff_vec_shift, multi_array<double, 3> &b_coeff_vec, const lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, partition_info<1> partition1, partition_info<2> partition2, Index mu, vector<Index> sigma2, const vector<multi_array<double, 2>> &w_x_dep)
{
    multi_array<double, 2> b_coeff({grid.r, grid.r});
    multi_array<double, 2> b_coeff_shift({grid.r, grid.r});
    multi_array<double, 1> w_x2({grid.dx2});
    multi_array<double, 2> xx2_shift(lr_sol.V.shape());

    vector<Index> vec_index(grid.m2);
    multi_array<Index, 1> vec_index_start({grid.m2});
    std::fill(vec_index.begin(), vec_index.end(), 0);
    std::fill(vec_index_start.begin(), vec_index_start.end(), 0);

    // Calculate the shifted X2
    ShiftMultiArrayRows<2>(xx2_shift, lr_sol.V, -sigma2[mu], reaction_system.reactions[mu]->minus_nu, grid);

    for (Index alpha1_dep = 0; alpha1_dep < partition1.dx_dep(mu); alpha1_dep++)
    {
#ifdef __OPENMP__
#pragma omp parallel firstprivate(vec_index_start, vec_index)
#endif
        {
            Index alpha2_dep;
#ifdef __OPENMP__
            SetVecIndexStart(vec_index_start, grid.n2, grid.dx2);
#endif
            std::copy(vec_index_start.begin(), vec_index_start.end(), vec_index.begin());

#ifdef __OPENMP__
#pragma omp for schedule(static)
#endif
            for (Index alpha2 = 0; alpha2 < grid.dx2; alpha2++)
            {
                alpha2_dep = VecIndextoDepCombIndex(vec_index, partition2.n_dep[mu], partition2.dep_vec[mu]);
                w_x2(alpha2) = w_x_dep[mu](alpha1_dep, alpha2_dep) * grid.h2_mult;
                IncrVecIndex(vec_index, grid.n2, grid.m2);
            }
        }

        coeff(xx2_shift, lr_sol.V, w_x2, b_coeff_shift, blas);
        coeff(lr_sol.V, lr_sol.V, w_x2, b_coeff, blas);

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


void CalculateCoefficientsS(multi_array<double, 5> &e_coeff_tot, multi_array<double, 5> &f_coeff_tot, vector<Index> sigma1, vector<Index> sigma2, const lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, partition_info<1> partition1, partition_info<2> partition2, const vector<multi_array<double, 2>> &w_x_dep)
{
    multi_array<double, 2> e_coeff({grid.r, grid.r});
    multi_array<double, 2> f_coeff({grid.r, grid.r});
    multi_array<double, 1> w_x1({grid.dx1});
    multi_array<double, 1> w_x1_shift({grid.dx1});
    multi_array<double, 2> xx1_shift(lr_sol.X.shape());

    vector<Index> vec_index(grid.m1);
    multi_array<Index, 1> vec_index_start({grid.m1});
    std::fill(vec_index.begin(), vec_index.end(), 0);
    std::fill(vec_index_start.begin(), vec_index_start.end(), 0);

    multi_array<double, 3> b_coeff_vec;
    multi_array<double, 3> b_coeff_vec_shift;

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        // Calculate the shifted X1
        ShiftMultiArrayRows<1>(xx1_shift, lr_sol.X, -sigma1[mu], reaction_system.reactions[mu]->minus_nu, grid);
        b_coeff_vec.resize({partition1.dx_dep(mu), grid.r, grid.r});
        b_coeff_vec_shift.resize({partition1.dx_dep(mu), grid.r, grid.r});
        CalculateCoefficientsB(b_coeff_vec_shift, b_coeff_vec, lr_sol, blas, reaction_system, grid, partition1, partition2, mu, sigma2, w_x_dep);

        for (Index j = 0; j < grid.r; j++)
        {
            for (Index l = 0; l < grid.r; l++)
            {
#ifdef __OPENMP__
#pragma omp parallel firstprivate(vec_index_start, vec_index)
#endif
                {
                    Index alpha1_dep;
#ifdef __OPENMP__
                    SetVecIndexStart(vec_index_start, grid.n1, grid.dx1);
#endif
                    std::copy(vec_index_start.begin(), vec_index_start.end(), vec_index.begin());

#ifdef __OPENMP__
#pragma omp for schedule(static)
#endif
                    // Calculate integration weights
                    for (Index alpha1 = 0; alpha1 < grid.dx1; alpha1++)
                    {
                        // get_time::start("sweightvec_index");
                        alpha1_dep = VecIndextoDepCombIndex(vec_index, partition1.n_dep[mu], partition1.dep_vec[mu]);
                        // get_time::stop("sweightvec_index");
                        w_x1_shift(alpha1) = b_coeff_vec_shift(alpha1_dep, j, l) * grid.h1_mult;
                        w_x1(alpha1) = b_coeff_vec(alpha1_dep, j, l) * grid.h1_mult;
                        IncrVecIndex(vec_index, grid.n1, grid.m1);
                    }
                }

                coeff(xx1_shift, lr_sol.X, w_x1_shift, e_coeff, blas);
                coeff(lr_sol.X, lr_sol.X, w_x1, f_coeff, blas);

                for (Index i = 0; i < grid.r; i++)
                {
                    for (Index k = 0; k < grid.r; k++)
                    {
                        e_coeff_tot(mu, i, j, k, l) = e_coeff(i, k);
                        f_coeff_tot(mu, i, j, k, l) = f_coeff(i, k);
                    }
                }
            }
        }
    }
}


// Perform S-Step with time step size `tau`
void PerformSStep(multi_array<double, 2> &s_dot, const multi_array<double, 2> &s, const multi_array<double, 5> &e_coeff, const multi_array<double, 5> &f_coeff, vector<Index> sigma1, std::vector<Index> sigma2, blas_ops blas, mysys reaction_system, grid_info grid, partition_info<1> partition1, partition_info<2> partition2, const vector<multi_array<double, 2>> &w_x_dep, double tau)
{
    // multi_array<double, 2> s_dot(lr_sol.S.shape());
    set_zero(s_dot);

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
#ifdef __OPENMP__
#pragma omp parallel for collapse(4)
#endif
        for (Index i = 0; i < grid.r; i++)
        {
            for (Index j = 0; j < grid.r; j++)
            {
                for (Index k = 0; k < grid.r; k++)
                {
                    for (Index l = 0; l < grid.r; l++)
                    {
                        s_dot(i, j) += tau * s(k, l) * e_coeff(mu, i, j, k, l);
                        s_dot(i, j) -= tau * s(k, l) * f_coeff(mu, i, j, k, l);
                    }
                }
            }
        }
    }
    // lr_sol.S -= s_dot;
}