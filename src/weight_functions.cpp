#include "weight_functions.hpp"

using std::vector;

// TODO: Only use two weight vectors w_x1 and w_x2 in the program
void CalculateWeight(multi_array<double, 3> &w_x, mysys reaction_system, grid_info grid)
{
    Index comb_index_tot;
    vector<double> state(grid.m1 + grid.m2, 0.0);

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        for (Index alpha1 = 0; alpha1 < grid.dx1; alpha1++)
        {
            // Calculate integration weight
            for (Index alpha2 = 0; alpha2 < grid.dx2; alpha2++)
            {
                comb_index_tot = alpha1 + grid.dx1 * alpha2;
                CombIndexToState(state, comb_index_tot, grid.n, grid.lim);

                w_x(alpha1, alpha2, mu) = reaction_system.reactions[mu]->propensity(state);
            }
        }
    }
}


void CalculateWeightDep(vector<multi_array<double, 2>> &w_x_dep, mysys reaction_system, grid_info grid, partition_info partition)
{
    Index comb_index_tot, comb_index1, comb_index2;
    vector<double> state(grid.m1 + grid.m2, 0.0);

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        w_x_dep[mu].resize({partition.dx_dep1(mu), partition.dx_dep2(mu)});
        for (Index alpha1 = 0; alpha1 < partition.dx_dep1(mu); alpha1++)
        {
            comb_index1 = DepCombIndexToCombIndex(alpha1, partition.n_dep1[mu], grid.n1, partition.dep_vec1[mu]);
            // Calculate integration weight
            for (Index alpha2 = 0; alpha2 < partition.dx_dep2(mu); alpha2++)
            {
                comb_index2 = DepCombIndexToCombIndex(alpha2, partition.n_dep2[mu], grid.n2, partition.dep_vec2[mu]);
                comb_index_tot = comb_index1 + grid.dx1 * comb_index2;
                CombIndexToState(state, comb_index_tot, grid.n, grid.lim);
                w_x_dep[mu](alpha1, alpha2) = reaction_system.reactions[mu]->propensity(state);
            }
        }
    }
}