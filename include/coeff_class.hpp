#ifndef COEFF_CLASS_HPP
#define COEFF_CLASS_HPP

#include<memory>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "grid_class.hpp"

struct cme_coeff
{
    std::vector<std::vector<double>> propensity;
    multi_array<double, 3> A, B;
    multi_array<double, 4> E, F;

    cme_coeff(const Index _n_reactions, const Index _r_in)
    : propensity(_n_reactions) 
    , A({_n_reactions, _r_in, _r_in})
    , B({_n_reactions, _r_in, _r_in})
    , E({_r_in, _r_in, _r_in, _r_in})
    , F({_r_in, _r_in, _r_in, _r_in})
    {}
};

struct cme_internal_coeff
{
    multi_array<double, 6> G, H;

    cme_internal_coeff(const Index _r_in, const std::array<Index, 2> _r_out) 
    : G({_r_in, _r_out[0], _r_out[1], _r_in, _r_out[0], _r_out[1]})
    , H({_r_in, _r_out[0], _r_out[1], _r_in, _r_out[0], _r_out[1]})
    {}
};

struct cme_external_coeff
{
    std::vector<multi_array<double, 3>> C, D;

    cme_external_coeff(const Index _n_reactions)
    : C(_n_reactions)
    , D(_n_reactions)
    {}
};

#endif