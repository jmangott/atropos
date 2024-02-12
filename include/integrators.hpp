#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP

#include <generic/matrix.hpp>

#include "integration_methods.hpp"
#include "timer_class.hpp"
#include "tree_class.hpp"

struct integrator_base
{
    integrator_base(const blas_ops &_blas, const std::map<std::string, integration_method*> &_integration_methods)
    : blas(_blas)
    , integration_methods(_integration_methods)
    {}

    const blas_ops blas;
    const std::map<std::string, integration_method*> integration_methods;
};

struct ttn_integrator : integrator_base
{
    ttn_integrator(const blas_ops &_blas, const std::map<std::string, integration_method*> &_integration_methods)
    : integrator_base(_blas, _integration_methods)
    {}

    void operator()(cme_internal_node * const node, const double tau) const;

    template <Index id>
    void SubflowPhi(cme_internal_node *const node, const double tau) const;

    void SubflowPsi(cme_internal_node *const node, const double tau) const;
};

#endif