#include "grid_class.hpp"

void grid_parms::Initialize()
{
    // dx, h_mult
    dx = 1;
    h_mult = 1.0;
    for (int i = 0; i < d; ++i) {
        dx *= n[i];
        h_mult *= (double)binsize[i];
    }

    // dx_dep, shift, idx_dep, n_dep
    for (Index mu = 0; mu < n_reactions; ++mu) {
        dx_dep[mu] = 1;
        shift[mu] = 0;
        int stride = 1;
        for (Index i = 0; i < d; ++i) {
            if (dep(mu, i) == true) {
                dx_dep[mu] *= n[i];
                idx_dep[mu].push_back(i);
                n_dep[mu].push_back(n[i]);
            }
            shift[mu] += nu(mu, i) * stride;
            stride *= n[i];
        }
    }
}
