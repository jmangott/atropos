#ifndef KL_STEP_FUNCTIONS_HPP
#define KL_STEP_FUNCTIONS_HPP

#include <array>
#include <memory>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"
#include "partition_class.hpp"
#include "reaction_class.hpp"
#include "timer_class.hpp"


// Calculate coefficients C1 and D1 (`id` = 1) or C2 and D2 (`id` = 2) for all values of `dep_vec` for a given reaction `mu`
template <Index id>
void CalculateCoefficientsKL(std::vector<multi_array<double, 3>> &c_coeff_dep, std::vector<multi_array<double, 3>> &d_coeff_dep, const std::vector<Index> &sigma_c, const lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, partition_info<id> partition, partition_info<id == 1 ? 2 : 1> partition2, const std::vector<multi_array<double, 2>> &w_x_dep)
{
    multi_array<double, 2> c_coeff({grid.r, grid.r});
    multi_array<double, 2> d_coeff({grid.r, grid.r});
    // TODO: use a pointer for `tmp_xx`
    std::array<Index, 2> tmp_xx_dim;
    Index weight_dim;
    (id == 1) ? (tmp_xx_dim = lr_sol.V.shape()) : (tmp_xx_dim = lr_sol.X.shape());
    (id == 1) ? (weight_dim = grid.dx2) : (weight_dim = grid.dx1);
    multi_array<double, 2> tmp_xx(tmp_xx_dim), xx_shift(tmp_xx_dim);
    (id == 1) ? (tmp_xx = lr_sol.V) : (tmp_xx = lr_sol.X);
    multi_array<double, 1> weight({weight_dim});
    Index alpha1_dep;
    vector<Index> vec_index;

    (id == 1) ? vec_index.resize(grid.m2) : vec_index.resize(grid.m1);

    // TODO: write a custom `coeff` routine, such that the conversion to a `weight` vector with length dx1 or dx2 is no longer needed
    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        // get_time::start("shift_kl");
        ShiftMultiArrayRows<id == 1 ? 2 : 1>(xx_shift, tmp_xx, -sigma_c[mu], reaction_system.reactions[mu]->minus_nu, grid);
        // get_time::stop("shift_kl");

        for (Index alpha2_dep = 0; alpha2_dep < partition.dx_dep(mu); alpha2_dep++)
        {
            // get_time::start("weight_kl");
            std::fill(vec_index.begin(), vec_index.end(), 0);
#ifdef __OPENMP__
#pragma omp parallel firstprivate(vec_index) private(alpha1_dep)
#endif
            {
                if constexpr (id == 1)
                {
#ifdef __OPENMP__
#pragma omp for
#endif
                    for (Index alpha1 = 0; alpha1 < grid.dx2; alpha1++)
                    {
                        alpha1_dep = VecIndextoDepCombIndex(vec_index, partition2.n_dep[mu], partition2.dep_vec[mu]);
                        weight(alpha1) = w_x_dep[mu](alpha2_dep, alpha1_dep) * grid.h2_mult;
                        IncrVecIndex(vec_index, grid.n2, grid.m2);
                    }
                }
                else if constexpr (id == 2)
                {
#ifdef __OPENMP__
#pragma omp for
#endif
                    for (Index alpha1 = 0; alpha1 < grid.dx1; alpha1++)
                    {
                        alpha1_dep = VecIndextoDepCombIndex(vec_index, partition2.n_dep[mu], partition2.dep_vec[mu]);
                        weight(alpha1) = w_x_dep[mu](alpha1_dep, alpha2_dep) * grid.h1_mult;
                        IncrVecIndex(vec_index, grid.n1, grid.m1);
                    }
                }
            }
            // get_time::stop("weight_kl");

            // get_time::start("coeff_kl");
            coeff(xx_shift, tmp_xx, weight, c_coeff, blas);
            coeff(tmp_xx, tmp_xx, weight, d_coeff, blas);
            // get_time::stop("coeff_kl");

            // get_time::start("write_coeff_kl");
#ifdef __OPENMP__
#pragma omp parallel for collapse(2)
#endif
            for (Index i = 0; i < grid.r; i++)
            {
                for (Index j = 0; j < grid.r; j++)
                {
                    c_coeff_dep[mu](alpha2_dep, i, j) = c_coeff(i, j);
                    d_coeff_dep[mu](alpha2_dep, i, j) = d_coeff(i, j);
                }
            }
            // get_time::stop("write_coeff_kl");
        }
        // cout << "K (=1) or L (=2): " << id << endl;
        // cout << "mu: " << mu << endl;
        // cout << get_time::sorted_output() << endl;
        // cout << endl;
        // get_time::reset();
    }
}


// Perform K-Step (`id` = 1) or L-Step (`id` = 2) with time step size `tau`
template <Index id>
void PerformKLStep(multi_array<double, 2> &kl_dot, const multi_array<double, 2> &kl, const std::vector<multi_array<double, 3>> &c_coeff_dep, const std::vector<multi_array<double, 3>> &d_coeff_dep, const std::vector<Index> &sigma, blas_ops blas, mysys reaction_system, grid_info grid, partition_info<id> partition, const std::vector<multi_array<double, 2>> &w_x_dep, double tau)
{
    Index dx1 = kl.shape()[0];
    Index r = kl.shape()[1];
    multi_array<double, 2> prod_klc({dx1, r});
    multi_array<double, 2> prod_klc_shift({dx1, r});
    multi_array<double, 2> prod_kld({dx1, r});
    // multi_array<double, 2> kl_dot({dx1, r});
    set_zero(kl_dot);

    Index dim_n, alpha;

    (id == 1) ? dim_n = grid.n1.shape()[0] : dim_n = grid.n2.shape()[0];
    multi_array<Index, 1> n({dim_n});
    (id == 1) ? n = grid.n1 : n = grid.n2;

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        set_zero(prod_klc);
        set_zero(prod_kld);

#ifdef __OPENMP__
#pragma omp parallel for private(alpha)
#endif
        for (Index i = 0; i < dx1; i++)
        {
            alpha = CombIndexToDepCombIndex(i, partition.n_dep[mu], n, partition.dep_vec[mu]);

            // Calculate matrix-vector multiplication of C2*K and D2*K
            for (Index j = 0; j < grid.r; j++)
            {
                for (Index l = 0; l < grid.r; l++)
                {
                    prod_klc(i, j) += tau * kl(i, l) * c_coeff_dep[mu](alpha, j, l);
                    prod_kld(i, j) += tau * kl(i, l) * d_coeff_dep[mu](alpha, j, l);
                }
            }
        }
        // Shift prod_KC
        ShiftMultiArrayRows<id>(prod_klc_shift, prod_klc, sigma[mu], reaction_system.reactions[mu]->nu, grid);

        // Calculate k_dot = shift(C1,2 * K) - D1,2 * K
        kl_dot += prod_klc_shift;
        kl_dot -= prod_kld;
    }
    // kl += kl_dot;
}

#endif