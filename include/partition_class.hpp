#ifndef PARTITION_CLASS_HPP
#define PARTITION_CLASS_HPP

#include<iostream>
#include<vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "grid_class.hpp"
#include "reaction_class.hpp"

// TODO: Write tests for this object

struct partition_base
{
    multi_array<Index, 1> dx_dep;
    multi_array<Index, 1> dx_rem;
    std::vector<std::vector<Index>> n_dep;
    std::vector<std::vector<Index>> n_rem;
    std::vector<std::vector<Index>> dep_vec;

    partition_base(grid_info _grid, mysys _reaction_system) : dx_dep({_reaction_system.mu()}), dx_rem({_reaction_system.mu()}), n_dep(_reaction_system.mu()), n_rem(_reaction_system.mu()), dep_vec(_reaction_system.mu()){};

    void partition_common_init(grid_info grid, mysys reaction_system, Index inc)
    {
        std::vector<Index> dep_vec_tot;
        std::vector<Index> vec_index_dep;
        multi_array<Index, 1> vec_index_rem({grid.m1});
        multi_array<Index, 1> vec_index_zero({grid.m1});

        for (Index mu = 0; mu < reaction_system.mu(); mu++)
        {
            dep_vec_tot = reaction_system.reactions[mu]->depends_on;
            dx_dep(mu) = 1;
            dx_rem(mu) = 1;
            std::copy(grid.n1.begin(), grid.n1.end(), std::back_inserter(n_rem[mu]));

            for (Index i = 0; i < grid.m1; i++)
                vec_index_zero(i) = 0;

            for (auto const &ele : dep_vec_tot)
            {
                // convert indices (m1, m1 + 1, ..., m1 + m2) to (0, 1, ..., m2 - 1)
                auto ele_inc = ele - inc;
                if ((0 <= ele_inc) && (ele_inc < grid.m1))
                {
                    dep_vec[mu].push_back(ele_inc);
                    dx_dep(mu) *= grid.n(ele_inc);
                    n_dep[mu].push_back(grid.n(ele_inc));
                    n_rem[mu][ele_inc] = 1;
                }
            }

            for (auto const &ele : n_rem[mu])
            {
                dx_rem(mu) *= ele;
            }
        }
    }
};

template <Index id>
struct partition_info {};

template <>
struct partition_info<1> : partition_base
{
    partition_info(grid_info _grid, mysys _reaction_system) : partition_base(_grid, _reaction_system)
    {
        partition_base::partition_common_init(_grid, _reaction_system, 0);
    }
};

template <>
struct partition_info<2> : partition_base
{
    partition_info(grid_info _grid, mysys _reaction_system) : partition_base(_grid, _reaction_system)
    {
        grid_info grid_alt(_grid.m2, _grid.m1, _grid.r, _grid.n2, _grid.n1, _grid.binsize2, _grid.binsize1, _grid.liml2, _grid.liml1);
        partition_base::partition_common_init(grid_alt, _reaction_system, _grid.m1);
    }
};

#endif