#include "weight_functions.hpp"

using std::vector;

// TODO: Only use two weight vectors w_x1 and w_x2 in the program
void CalculateWeightDep(vector<multi_array<double, 2>> &w_x_dep, double &min_prop, double &max_prop, mysys reaction_system, grid_info grid, partition_info<1> partition1, partition_info<2> partition2)
{
    Index comb_index_tot, comb_index1, comb_index2;
    vector<double> state(grid.d, 0.0);
    double propensity;

    // Initialize `min_prop` and `max_prop`
    comb_index_tot = 0;
    CombIndexToState(state, comb_index_tot, grid.n, grid.liml, grid.binsize);
    for (Index i = 0; i < grid.d; i++)
        state[i] += (grid.binsize(i) - 1.0) * 0.5;
    propensity = reaction_system.reactions[0]->propensity(state);
    min_prop = propensity;
    max_prop = propensity;

    // Set up `w_x_dep`
    for (Index mu = 0; mu < reaction_system.mu(); mu++)
    {
        for (Index alpha1 = 0; alpha1 < partition1.dx_dep(mu); alpha1++)
        {
            comb_index1 = DepCombIndexToCombIndex(alpha1, partition1.n_dep[mu], grid.n1, partition1.dep_vec[mu]);
            for (Index alpha2 = 0; alpha2 < partition2.dx_dep(mu); alpha2++)
            {
                comb_index2 = DepCombIndexToCombIndex(alpha2, partition2.n_dep[mu], grid.n2, partition2.dep_vec[mu]);
                comb_index_tot = comb_index1 + grid.dx1 * comb_index2;
                CombIndexToState(state, comb_index_tot, grid.n, grid.liml, grid.binsize);

                // Calculate propensity for mid of bin
                for (Index i = 0; i < grid.d; i++)
                    state[i] += (grid.binsize(i) - 1.0) * 0.5;
                propensity = reaction_system.reactions[mu]->propensity(state);

                w_x_dep[mu](alpha1, alpha2) = propensity;
                if (propensity < min_prop) min_prop = propensity;
                if (propensity > max_prop) max_prop = propensity;
            }
        }
    }
}