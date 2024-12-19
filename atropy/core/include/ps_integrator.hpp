#ifndef PS_INTEGRATOR_HPP
#define PS_INTEGRATOR_HPP

#include "integrators.hpp"

struct ps_integrator : integrator_base {
    ps_integrator(
        const Ensign::Matrix::blas_ops& _blas,
        const std::map<std::string, integration_method*>& _integration_methods)
        : integrator_base(_blas, _integration_methods)
    {
    }

    void operator()(cme_internal_node* const node, const double tau) const;

    template <Index id>
    void SubflowPhi(cme_internal_node* const node, const double tau) const;

    void SubflowPsi(cme_internal_node* const node, const double tau) const;
};

#endif
