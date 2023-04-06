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

    // TODO: write a custom `coeff` routine, such that the conversion to a `weight` vector with length dx1 or dx2 is no longer needed
#pragma omp parallel for
    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        ShiftMultiArrayRows<id == 1 ? 2 : 1>(xx_shift, tmp_xx, -sigma_c[mu], reaction_system.reactions[mu]->minus_nu, grid);

        for (Index alpha2_dep = 0; alpha2_dep < partition.dx_dep(mu); alpha2_dep++)
        {
            if (id == 1)
            {
                for (Index alpha1 = 0; alpha1 < grid.dx2; alpha1++)
                {
                    alpha1_dep = CombIndexToDepCombIndex(alpha1, partition2.n_dep[mu], grid.n2, partition2.dep_vec[mu]);
                    weight(alpha1) = w_x_dep[mu](alpha2_dep, alpha1_dep) * grid.h2_mult;
                }
            }
            else if (id == 2)
            {
                for (Index alpha1 = 0; alpha1 < grid.dx1; alpha1++)
                {
                    alpha1_dep = CombIndexToDepCombIndex(alpha1, partition2.n_dep[mu], grid.n1, partition2.dep_vec[mu]);
                    weight(alpha1) = w_x_dep[mu](alpha1_dep, alpha2_dep) * grid.h1_mult;
                }
            }
            else
            {
                std::cerr << "ERROR: id must be 1 (=L) or 2 (=K)!" << endl;
                std::abort();
            }

            coeff(xx_shift, tmp_xx, weight, c_coeff, blas);
            coeff(tmp_xx, tmp_xx, weight, d_coeff, blas);

            for (Index i = 0; i < grid.r; i++)
            {
                for (Index j = 0; j < grid.r; j++)
                {
                    c_coeff_dep[mu](alpha2_dep, i, j) = c_coeff(i, j);
                    d_coeff_dep[mu](alpha2_dep, i, j) = d_coeff(i, j);
                }
            }
        }
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

    Index alpha;
    Index dim_n;

    (id == 1) ? dim_n = grid.n1.shape()[0] : dim_n = grid.n2.shape()[0];
    multi_array<Index, 1> n({dim_n});
    (id == 1) ? n = grid.n1 : n = grid.n2;

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        set_zero(prod_klc);
        set_zero(prod_kld);

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