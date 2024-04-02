#ifndef COEFF_CLASS_HPP
#define COEFF_CLASS_HPP

#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "grid_class.hpp"

// TODO: rename A->C, B->D, A_bar->A and B_bar->B
// TODO: delete F and H coefficients
struct cme_coeff
{
    std::vector<multi_array<double, 2>> A, B, A_bar, B_bar;
    multi_array<double, 4> E, F;

    cme_coeff(const Index _n_reactions, const Index _r_in)
    : A(_n_reactions)
    , B(_n_reactions)
    , A_bar(_n_reactions)
    , B_bar(_n_reactions)
    , E({_r_in, _r_in, _r_in, _r_in})
    , F({_r_in, _r_in, _r_in, _r_in})
    {
        for (Index mu = 0; mu < _n_reactions; ++mu)
        {
            A[mu].resize({_r_in, _r_in});
            B[mu].resize({_r_in, _r_in});
            A_bar[mu].resize({_r_in, _r_in});
            B_bar[mu].resize({_r_in, _r_in});
        }
    }
};

struct cme_internal_coeff
{
    multi_array<double, 4> G, H;

    cme_internal_coeff(const Index _r_in, const std::array<Index, 2> _r_out) 
    : G({prod(_r_out), prod(_r_out), _r_in, _r_in})
    , H({prod(_r_out), prod(_r_out), _r_in, _r_in})
    {}
};

struct cme_external_coeff
{
    std::vector<std::vector<double>> propensity;

    cme_external_coeff(const Index _n_reactions)
    : propensity(_n_reactions)
    {}
};

#endif