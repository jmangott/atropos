#include "kl_step_functions.hpp"
#include <chrono>

using std::vector;

void CalculateCoefficientsX(int id, multi_array<double, 2> &c_coeff, multi_array<double, 2> &d_coeff, const lr2<double> &lr_sol, blas_ops blas, const multi_array<double, 2> &xx_shift, Index alpha2, mysys reaction_system, grid_info grid, Index mu, multi_array<double, 3> &w_x)
{
    std::array<Index, 2> tmp_xx_dim;
    Index weight_dim;
    (id == 1) ? (tmp_xx_dim = lr_sol.V.shape()) : (tmp_xx_dim = lr_sol.X.shape());
    (id == 1) ? (weight_dim = grid.dx2) : (weight_dim = grid.dx1);
    multi_array<double, 2> tmp_xx(tmp_xx_dim);
    multi_array<double, 1> weight({weight_dim});

    if (id == 1)
    {
        tmp_xx = lr_sol.V;
        for (Index alpha1 = 0; alpha1 < grid.dx2; alpha1++)
        {
            weight(alpha1) = w_x(alpha2, alpha1, mu) * grid.h2_mult;
        }
    }
    else if (id == 2)
    {
        tmp_xx = lr_sol.X;
        for (Index alpha1 = 0; alpha1 < grid.dx1; alpha1++)
        {
            weight(alpha1) = w_x(alpha1, alpha2, mu) * grid.h1_mult;
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


// TODO: Write tests for this object
partition_info::partition_info(grid_info grid, mysys reaction_system) : dx_dep1({reaction_system.mu()}), dx_dep2({reaction_system.mu()}), dx_rem1({reaction_system.mu()}), dx_rem2({reaction_system.mu()}), n_dep1(reaction_system.mu()), n_dep2(reaction_system.mu()), n_rem1(reaction_system.mu()), n_rem2(reaction_system.mu()), dep_vec1(reaction_system.mu()), dep_vec2(reaction_system.mu())
{
    vector<Index> dep_vec_tot;
    vector<Index> vec_index_dep1, vec_index_dep2;
    multi_array<Index, 1> vec_index_rem1({grid.m1}), vec_index_rem2({grid.m2});
    multi_array<Index, 1> vec_index_zero1({grid.m1}), vec_index_zero2({grid.m2});

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        dep_vec_tot = reaction_system.reactions[mu]->depends_on;
        dx_dep1(mu) = 1;
        dx_dep2(mu) = 1;
        dx_rem1(mu) = 1;
        dx_rem2(mu) = 1;
        std::copy(grid.n1.begin(), grid.n1.end(), std::back_inserter(n_rem1[mu]));
        std::copy(grid.n2.begin(), grid.n2.end(), std::back_inserter(n_rem2[mu]));

        for (Index i = 0; i < grid.m1; i++)
            vec_index_zero1(i) = 0;
        for (Index i = 0; i < grid.m2; i++)
            vec_index_zero2(i) = 0;

        for (auto const &ele : dep_vec_tot)
        {
            // convert indices (m1, m1 + 1, ..., m1 + m2) to (0, 1, ..., m2 - 1)
            if ((0 <= ele) && (ele < grid.m1))
            {
                dep_vec1[mu].push_back(ele);
                dx_dep1(mu) *= grid.n1(ele);
                n_dep1[mu].push_back(grid.n1(ele));
                n_rem1[mu][ele] = 1;
            }
            else if ((grid.m1 <= ele) && (ele < grid.m1 + grid.m2))
            {
                dep_vec2[mu].push_back(ele - grid.m1);
                dx_dep2(mu) *= grid.n2(ele - grid.m1);
                n_dep2[mu].push_back(grid.n2(ele - grid.m1));
                n_rem2[mu][ele - grid.m1] = 1;
            }
        }

        for (auto const &ele : n_rem1[mu])
        {
            dx_rem1(mu) *= ele;
        }
        for (auto const &ele : n_rem2[mu])
        {
            dx_rem2(mu) *= ele;
        }
    }
}