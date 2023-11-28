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

    cme_internal_coeff(const Index _n_reactions, const Index r_in, const Index r_out) 
    : coeff(_n_reactions)
    , E({r_out, r_out, r_out, r_out})
    , F({r_out, r_out, r_out, r_out})
    , G({r_in, r_out, r_out, r_in, r_out, r_out})
    , H({r_in, r_out, r_out, r_in, r_out, r_out})
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