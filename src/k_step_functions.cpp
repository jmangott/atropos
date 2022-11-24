#include "k_step_functions.hpp"

using std::vector;

multi_array<double, 1> CalculateWeightX(int id, multi_array<Index, 1> vec_index_c, mysys reaction_system, grid_info grid, Index mu)
{
    Index d = grid.m1 + grid.m2;
    grid_info *grid_alt;
    Index inc, inc_c;
    if (id == 1)
    {
        grid_alt = new grid_info(grid.m1, grid.m2, grid.r, grid.n1, grid.n2, grid.k1, grid.k2);
        inc = 0;
        inc_c = grid.m1;
    }
    else if (id == 2)
    {
        grid_alt = new grid_info(grid.m2, grid.m1, grid.r, grid.n2, grid.n1, grid.k2, grid.k1);
        inc = grid.m1;
        inc_c = 0;
    }
    else
    {
        std::cerr << "ERROR: id must be 0 or 1!" << endl;
        std::abort();
    }
    
    multi_array<double, 1> a_vec({grid_alt->dx1});
    multi_array<double, 1> state_c({grid_alt->m2});
    multi_array<double, 1> state({grid_alt->m1});
    vector<double> state_tot(d, 0.0);
    multi_array<Index, 1> vec_index({grid_alt->m1});
    multi_array<double, 1> w_x({grid_alt->dx1});

    state_c = VecIndexToState(vec_index_c, grid_alt->n2, grid_alt->lim2);

    for (Index i = 0; i < grid_alt->m2; i++)
    {
        state_tot[i + inc_c] = state_c(i);
    }

    // Calculate propensity vector
    for (Index alpha = 0; alpha < grid_alt->dx1; alpha++)
    {
        vec_index = CombIndexToVecIndex(alpha, grid_alt->n1);
        state = VecIndexToState(vec_index, grid_alt->n1, grid_alt->lim1);

        // `state_tot` is the concatenation of `state1` and `state2`
        for (Index i = 0; i < grid_alt->m1; i++)
            state_tot[i + inc] = state(i);
        a_vec(alpha) = reaction_system.reactions[mu]->propensity(state_tot);
    }

    // Calculate integration weight
    double h_xx_mult = 1;
    for (Index i = 0; i < grid_alt->m1; i++)
        h_xx_mult *= grid_alt->h1(i);
    for (Index i = 0; i < grid_alt->dx1; i++)
        w_x(i) = a_vec(i) * h_xx_mult;

    delete grid_alt;
    return w_x;
}


void CalculateCoefficientsX(int id, multi_array<double, 2> &c2, multi_array<double, 2> &d2, const lr2<double> &lr_sol, blas_ops blas, const multi_array<double, 2> &xx2_shift, multi_array<Index, 1> vec_index1, mysys reaction_system, grid_info grid, Index mu)
{
    multi_array<double, 1> w_x2({grid.dx2});
    w_x2 = CalculateWeightX(id, vec_index1, reaction_system, grid, mu);

    coeff(xx2_shift, lr_sol.V, w_x2, c2, blas);
    coeff(lr_sol.V, lr_sol.V, w_x2, d2, blas);
}


void PerformKLStep(int id, vector<Index> sigma1, vector<Index> sigma2, lr2<double> &lr_sol, blas_ops blas, mysys reaction_system, grid_info grid, double tau)
{
    grid_info *grid_alt;
    multi_array<double, 2> xx_shift;
    multi_array<double, 2> tmp_xx;
    vector<Index> sigma, sigma_c;
    if (id == 1)
    {
        grid_alt = new grid_info(grid.m1, grid.m2, grid.r, grid.n1, grid.n2, grid.k1, grid.k2);
        xx_shift.resize(lr_sol.X.shape());
        tmp_xx.resize(lr_sol.X.shape());
        tmp_xx = lr_sol.X;
        sigma = sigma1;
        sigma_c = sigma2;
    }
    else if (id == 2)
    {
        grid_alt = new grid_info(grid.m2, grid.m1, grid.r, grid.n2, grid.n1, grid.k2, grid.k1);
        xx_shift.resize(lr_sol.V.shape());
        tmp_xx.resize(lr_sol.V.shape());
        // cout << xx_shift.shape()[0] << xx_shift.shape()[1];
        tmp_xx = lr_sol.V;
        sigma = sigma2;
        sigma_c = sigma1;
    }
    else
    {
        std::cerr << "ERROR: id must be 0 or 1!" << endl;
        std::abort();
    }

    multi_array<double, 2> c_coeff({grid.r, grid.r});
    multi_array<double, 2> d_coeff({grid.r, grid.r});

    vector<Index> dep_vec_tot, dep_vec_c, dep_vec;
    vector<Index> n_xx_c_coeff;
    multi_array<Index, 1> n_xx_c_reduced({grid_alt->m2});
    vector<Index> vec_index_c_coeff;
    multi_array<Index, 1> vec_index_c_reduced({grid_alt->m2});
    multi_array<Index, 1> vec_index_c_zero({grid_alt->m2});
    vector<double> lim_xx_c_coeff;
    multi_array<double, 1> lim_xx_c_reduced({grid_alt->m2});
    Index dxx_c_coeff_mult, dxx_c_reduced_mult;
    Index alpha_c;

    multi_array<double, 2> prod_KC({grid_alt->dx2, grid_alt->r});
    multi_array<double, 2> prod_KC_shift({grid_alt->dx2, grid_alt->r});
    multi_array<double, 2> prod_KD({grid_alt->dx2, grid_alt->r});
    multi_array<double, 2> k_dot({grid_alt->dx2, grid_alt->r});
    set_zero(k_dot);

    for (std::vector<myreact *>::size_type mu = 0; mu < reaction_system.mu(); mu++)
    {
        // Shift X2 for calculation of the coefficients
        ShiftMultiArrayRows(xx_shift, tmp_xx, -sigma[mu]);

        dep_vec_tot = reaction_system.reactions[mu]->depends_on;
        dep_vec_c.clear();
        dep_vec.clear();
        dxx_c_coeff_mult = 1;
        dxx_c_reduced_mult = 1;
        n_xx_c_coeff.clear();
        n_xx_c_reduced = grid_alt->n2;
        lim_xx_c_coeff.clear();
        lim_xx_c_reduced = grid_alt->lim2;
        set_zero(prod_KC);
        set_zero(prod_KD);

        for (Index i = 0; i < grid.m1; i++)
            vec_index_c_zero(i) = 0;

        for (auto const &ele : dep_vec_tot)
        {
            if (ele < grid.m1)
            {
                dep_vec_c.push_back(ele);
                dxx_c_coeff_mult *= grid.n1(ele);
                n_xx_c_coeff.push_back(grid.n1(ele));
                n_xx_c_reduced(ele) = 1;
                lim_xx_c_coeff.push_back(grid.lim1(ele));
                lim_xx_c_reduced(ele) = 0.0;
            }
            else
            {
                dep_vec.push_back(ele);
            }
        }

        for (auto const &ele : n_xx_c_reduced)
        {
            dxx_c_reduced_mult *= ele;
        }

        // Loop through all species in the complement partition on which the propensity for reaction mu depends
        for (Index i = 0; i < dxx_c_coeff_mult; i++)
        {
            vec_index_c_coeff = CombIndexToVecIndex(i, n_xx_c_coeff);

            // Convert vec_index_c_coeff to a vector with size m1
            for (vector<Index>::size_type j = 0; j < dep_vec_c.size(); j++)
                vec_index_c_zero(dep_vec_c[j]) = vec_index_c_coeff[j];

            CalculateCoefficientsX(2, c_coeff, d_coeff, lr_sol, blas, xx_shift, vec_index_c_zero, reaction_system, grid, mu);

            // Loop through the remaining species in the complement partition
            for (Index k = 0; k < dxx_c_reduced_mult; k++)
            {
                vec_index_c_reduced = CombIndexToVecIndex(k, n_xx_c_reduced);

                // vec_index_c_reduced contains now the real population number
                for (vector<Index>::size_type l = 0; l < dep_vec_c.size(); l++)
                    vec_index_c_reduced(dep_vec_c[l]) = vec_index_c_coeff[l];
                alpha_c = VecIndexToCombIndex(vec_index_c_reduced, grid_alt->n2);

                // Calculate matrix-vector multiplication of C2*K and D2*K
                for (Index j = 0; j < grid.r; j++)
                {
                    for (Index l = 0; l < grid.r; l++)
                    {
                        prod_KC(alpha_c, j) += tau * c_coeff(j, l) * lr_sol.X(alpha_c, l);
                        prod_KD(alpha_c, j) += tau * d_coeff(j, l) * lr_sol.X(alpha_c, l);
                    }
                }
            }
        }
        // Shift prod_KC
        ShiftMultiArrayRows(prod_KC_shift, prod_KC, sigma_c[mu]);

        // Calculate k_dot = shift(C2 * K) - D2 * K
        k_dot += prod_KC_shift;
        k_dot -= prod_KD;
    }
    lr_sol.X += k_dot;
    delete grid_alt;
}