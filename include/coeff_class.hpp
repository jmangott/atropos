#ifndef COEFF_CLASS_HPP
#define COEFF_CLASS_HPP

#include<memory>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "grid_class.hpp"


struct coeff
{
    std::vector<std::vector<double>> propensity;

    explicit coeff(const Index _n_reactions) : propensity(_n_reactions) {}
};

struct cme_internal_coeff : coeff
{
    multi_array<double, 4> E, F;
    multi_array<double, 6> G, H;

    cme_internal_coeff(const Index _n_reactions, const Index r_in, const std::array<Index, 2> r_out) 
    : coeff(_n_reactions)
    , E({r_out[0], r_out[1], r_out[0], r_out[1]})
    , F({r_out[0], r_out[1], r_out[0], r_out[1]})
    , G({r_in, r_out[0], r_out[1], r_in, r_out[0], r_out[1]})
    , H({r_in, r_out[0], r_out[1], r_in, r_out[0], r_out[1]})
    {}
};

struct cme_external_coeff : coeff
{
    std::vector<multi_array<double, 3>> C, D;
    
    explicit cme_external_coeff(const Index n_reactions)
    : coeff(n_reactions)
    , C(n_reactions)
    , D(n_reactions)
    {}
};

#endif