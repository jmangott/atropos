#ifndef BUG_INTEGRATOR_HPP
#define BUG_INTEGRATOR_HPP

#include "integrators.hpp"

struct bug_integrator : integrator_base {
    bug_integrator(
        const Ensign::Matrix::blas_ops& _blas,
        const std::map<std::string, integration_method*>& _integration_methods,
        const double _theta)
        : integrator_base(_blas, _integration_methods), theta(_theta)
    {
    }

    void operator()(cme_internal_node* const node, const double tau) const;

    template <Index id>
    void SubflowPhi(cme_internal_node* const node, const double tau) const;

    void SubflowPsi(cme_internal_node* const node, const double tau) const;

    void SubflowTheta(cme_internal_node* const node) const;

  private:
    const double theta;
};

#endif
