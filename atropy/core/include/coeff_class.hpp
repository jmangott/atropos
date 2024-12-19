#ifndef COEFF_CLASS_HPP
#define COEFF_CLASS_HPP

#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "grid_class.hpp"

// TODO: rename A->a, B->b, A_bar->A and B_bar->B
struct cme_coeff {
    std::vector<Ensign::multi_array<double, 2>> A, B, A_bar, B_bar;

    cme_coeff(const Index _n_reactions, const Index _r_in)
        : A(_n_reactions), B(_n_reactions), A_bar(_n_reactions), B_bar(_n_reactions)
    {
        for (Index mu = 0; mu < _n_reactions; ++mu) {
            A[mu].resize({_r_in, _r_in});
            B[mu].resize({_r_in, _r_in});
            A_bar[mu].resize({_r_in, _r_in});
            B_bar[mu].resize({_r_in, _r_in});
        }
    }
};

#endif
