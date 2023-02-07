#include "weight_functions.hpp"

using std::vector;

// TODO: Only use two weight vectors w_x1 and w_x2 in the program
void CalculateWeightDep(vector<multi_array<double, 2>> &w_x_dep, mysys reaction_system, grid_info grid, partition_info<1> partition1, partition_info<2> partition2)
{
    Index comb_index_tot, comb_index1, comb_index2;
    vector<double> state(grid.m1 + grid.m2, 0.0);

    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        w_x_dep[mu].resize({partition1.dx_dep(mu), partition2.dx_dep(mu)});
        for (Index alpha1 = 0; alpha1 < partition1.dx_dep(mu); alpha1++)
        {
            comb_index1 = DepCombIndexToCombIndex(alpha1, partition1.n_dep[mu], grid.n1, partition1.dep_vec[mu]);
            // Calculate integration weight
            for (Index alpha2 = 0; alpha2 < partition2.dx_dep(mu); alpha2++)
            {
                comb_index2 = DepCombIndexToCombIndex(alpha2, partition2.n_dep[mu], grid.n2, partition2.dep_vec[mu]);
                comb_index_tot = comb_index1 + grid.dx1 * comb_index2;
                CombIndexToState(state, comb_index_tot, grid.n, grid.liml, grid.limr);
                w_x_dep[mu](alpha1, alpha2) = reaction_system.reactions[mu]->propensity(state);
            }
        }
    }
}