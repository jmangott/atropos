#include "s_step_functions.hpp"

using std::vector;

void CalculateCoefficientsS(multi_array<double, 5> &e_coeff_tot, multi_array<double, 5> &f_coeff_tot, const vector<multi_array<double, 3>> &c_coeff_dep, const vector<multi_array<double, 3>> &d_coeff_dep, vector<Index> sigma1, vector<Index> sigma2, const lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, partition_info<1> partition1, partition_info<2> partition2, const vector<multi_array<double, 2>> &w_x_dep)
{
    multi_array<double, 2> e_coeff({grid.r, grid.r});
    multi_array<double, 2> f_coeff({grid.r, grid.r});
    multi_array<double, 1> w_x1({grid.dx1});
    multi_array<double, 1> w_x1_shift({grid.dx1});
    multi_array<double, 2> xx1_shift(lr_sol.X.shape());

    vector<Index> vec_index(grid.m1);
    multi_array<Index, 1> vec_index_start({grid.m1});

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        // Calculate the shifted X1
        ShiftMultiArrayRows<1>(xx1_shift, lr_sol.X, -sigma1[mu], reaction_system.reactions[mu]->minus_nu, grid);

        for (Index j = 0; j < grid.r; j++)
        {
            for (Index l = 0; l < grid.r; l++)
            {
                std::fill(vec_index.begin(), vec_index.end(), 0);
                
#ifdef __OPENMP__
#pragma omp parallel firstprivate(vec_index)
#endif
                {
                    Index alpha1_dep;
#ifdef __OPENMP__
                    Index chunk_size = SetVecIndex(vec_index, grid.n1, grid.dx1);
#endif

#ifdef __OPENMP__
#pragma omp for schedule(static, chunk_size)
#endif
                    // Calculate integration weights
                    for (Index alpha1 = 0; alpha1 < grid.dx1; alpha1++)
                    {
                        // get_time::start("sweightvec_index");
                        alpha1_dep = VecIndexToDepCombIndex(vec_index, partition1.n_dep[mu], partition1.dep_vec[mu]);
                        // get_time::stop("sweightvec_index");
                        w_x1_shift(alpha1) = c_coeff_dep[mu](alpha1_dep, j, l) * grid.h1_mult;
                        w_x1(alpha1) = d_coeff_dep[mu](alpha1_dep, j, l) * grid.h1_mult;
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
    set_zero(s_dot);

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
// TODO: This works on the workstation only for OMP_NUM_THREADS=6 (or 5) if collapse(<3). Why!?
#ifdef __OPENMP__
#pragma omp parallel for collapse(2)
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
}