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
void CalculateCoefficientsX(multi_array<double, 2> &c_coeff, multi_array<double, 2> &d_coeff, const lr2<double> &lr_sol, blas_ops blas, const multi_array<double, 2> &xx_shift, Index alpha2_dep, mysys reaction_system, grid_info grid, partition_info<id == 1 ? 2 : 1> partition2, Index mu, std::vector<multi_array<double, 2>> &w_x_dep)
{
    // TODO: use a pointer for `tmp_xx`
    std::array<Index, 2> tmp_xx_dim;
    Index weight_dim;
    (id == 1) ? (tmp_xx_dim = lr_sol.V.shape()) : (tmp_xx_dim = lr_sol.X.shape());
    (id == 1) ? (weight_dim = grid.dx2) : (weight_dim = grid.dx1);
    multi_array<double, 2> tmp_xx(tmp_xx_dim);
    multi_array<double, 1> weight({weight_dim});
    Index alpha1_dep;

    // TODO: write a custom `coeff` routine, such that the conversion to a `weight` vector with length dx1 or dx2 is no longer needed
    if (id == 1)
    {
        tmp_xx = lr_sol.V;
        for (Index alpha1 = 0; alpha1 < grid.dx2; alpha1++)
        {
            alpha1_dep = CombIndexToDepCombIndex(alpha1, partition2.n_dep[mu], grid.n2, partition2.dep_vec[mu]);
            weight(alpha1) = w_x_dep[mu](alpha2_dep, alpha1_dep) * grid.h2_mult;
        }
    }
    else if (id == 2)
    {
        tmp_xx = lr_sol.X;
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
}


// Perform K-Step (`id` = 1) or L-Step (`id` = 2) with time step size `tau`
template <Index id>
void PerformKLStep(std::vector<Index> sigma1, std::vector<Index> sigma2, lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, partition_info<id> partition, partition_info<id == 1 ? 2 : 1> partition2, std::vector<multi_array<double, 2>> &w_x_dep, double tau)
{
    grid_info *grid_alt;
    std::array<Index, 2> tmp_xx_dim, tmp_xx_c_dim;
    if (id == 1)
    {
        tmp_xx_dim = lr_sol.V.shape();
        tmp_xx_c_dim = lr_sol.X.shape();
    }
    else if (id == 2)
    {
        tmp_xx_dim = lr_sol.X.shape();
        tmp_xx_c_dim = lr_sol.V.shape();
    }

    multi_array<double, 2> xx_shift(tmp_xx_dim);
    multi_array<double, 2> tmp_xx(tmp_xx_dim), tmp_xx_c(tmp_xx_c_dim);
    std::vector<Index> sigma, sigma_c;

    int id_c;
    if (id == 1)
    {
        grid_alt = new grid_info(grid.m1, grid.m2, grid.r, grid.n1, grid.n2, grid.k1, grid.k2, grid.liml1, grid.liml2);
        // TODO: `tmp_xx` and `tmp_xx_c` could be replaced by pointers
        tmp_xx = lr_sol.V;
        tmp_xx_c = lr_sol.X;
        sigma_c = sigma2;
        sigma = sigma1;
        id_c = 2;
    }
    else if (id == 2)
    {
        grid_alt = new grid_info(grid.m2, grid.m1, grid.r, grid.n2, grid.n1, grid.k2, grid.k1, grid.liml2, grid.liml1);
        tmp_xx = lr_sol.X;
        tmp_xx_c = lr_sol.V;
        sigma_c = sigma1;
        sigma = sigma2;
        id_c = 1;
    }
    else
    {
        std::cerr << "ERROR: `id` must be 1 (=K) or 2 (=L)!" << std::endl;
        std::abort();
    }

    multi_array<double, 2> c_coeff({grid_alt->r, grid_alt->r});
    multi_array<double, 2> d_coeff({grid_alt->r, grid_alt->r});

    // TODO: replace these quantities by lr_sol.V or lr_sol.X
    multi_array<double, 2> prod_KLC({grid_alt->dx1, grid_alt->r});
    multi_array<double, 2> prod_KLC_shift({grid_alt->dx1, grid_alt->r});
    multi_array<double, 2> prod_KLD({grid_alt->dx1, grid_alt->r});
    multi_array<double, 2> kl_dot({grid_alt->dx1, grid_alt->r});
    set_zero(kl_dot);

    Index alpha_c;
    std::vector<Index> vec_index_c_dep;

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        // Shift X1,2 for calculation of the coefficients
        ShiftMultiArrayRows(id_c, xx_shift, tmp_xx, -sigma_c[mu], reaction_system.reactions[mu]->minus_nu, grid);

        vec_index_c_dep.resize(partition.n_dep[mu].size());
        set_zero(prod_KLC);
        set_zero(prod_KLD);

        // Loop through all species in the partition on which the propensity for reaction mu depends
        for (Index i = 0; i < partition.dx_dep(mu); i++)
        {
            CombIndexToVecIndex(vec_index_c_dep, i, partition.n_dep[mu]);

            CalculateCoefficientsX<id>(c_coeff, d_coeff, lr_sol, blas, xx_shift, i, reaction_system, grid, partition2, mu, w_x_dep);
            // Loop through the remaining species in the partition
            for (Index k = 0; k < partition.dx_rem(mu); k++)
            {
                // TODO: Check if this function is really needed
                alpha_c = DepVecIndexRemCombIndexToCombIndex(vec_index_c_dep, k, partition.n_rem[mu], grid_alt->n1, partition.dep_vec[mu]);

                // Calculate matrix-vector multiplication of C2*K and D2*K
                for (Index j = 0; j < grid.r; j++)
                {
                    for (Index l = 0; l < grid.r; l++)
                    {
                        prod_KLC(alpha_c, j) += tau * tmp_xx_c(alpha_c, l) * c_coeff(j, l);
                        prod_KLD(alpha_c, j) += tau * tmp_xx_c(alpha_c, l) * d_coeff(j, l);
                    }
                }
            }
        }
        // Shift prod_KC
        ShiftMultiArrayRows(id, prod_KLC_shift, prod_KLC, sigma[mu], reaction_system.reactions[mu]->nu, grid);

        // Calculate k_dot = shift(C1,2 * K) - D1,2 * K
        kl_dot += prod_KLC_shift;
        kl_dot -= prod_KLD;
    }
    (id == 1) ? (lr_sol.X += kl_dot) : (lr_sol.V += kl_dot);
    delete grid_alt;
}

#endif