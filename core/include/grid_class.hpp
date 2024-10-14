#ifndef GRID_CLASS_HPP
#define GRID_CLASS_HPP

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

// struct for storing the grid parameters of one node of the hierarchical problem in a
// compact form
struct grid_parms {
    Index d;
    Index n_reactions;
    std::vector<Index> n;
    std::vector<Index> binsize;
    std::vector<double> liml;
    Ensign::multi_array<bool, 2> dep;
    Ensign::multi_array<Index, 2> nu;
    std::vector<int> species;

    Index dx = 1;
    double h_mult = 1.0;
    std::vector<Index> dx_dep;
    std::vector<Index> shift;
    std::vector<std::vector<Index>> idx_dep;
    std::vector<std::vector<Index>> n_dep;

    grid_parms(Index _d, Index _n_reactions)
        : d(_d), n_reactions(_n_reactions), n(_d), binsize(_d), liml(_d),
          dep({_n_reactions, _d}), nu({_n_reactions, _d}), species(_d),
          dx_dep(_n_reactions), shift(_n_reactions), idx_dep(_n_reactions),
          n_dep(_n_reactions) {};

    grid_parms(std::vector<Index> _n, std::vector<Index> _binsize,
               std::vector<double> _liml, Ensign::multi_array<bool, 2> _dep,
               Ensign::multi_array<Index, 2> _nu, std::vector<int> _species)
        : d(_n.size()), n_reactions(_dep.shape()[0]), n(_n), binsize(_binsize),
          liml(_liml), dep(_dep), nu(_nu), species(_species), dx_dep(n_reactions),
          shift(n_reactions), idx_dep(n_reactions), n_dep(n_reactions) {};

    void Initialize();
};

#endif
