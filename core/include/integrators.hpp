#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP

#include <generic/matrix.hpp>

#include "integration_methods.hpp"
#include "timer_class.hpp"
#include "tree_class.hpp"

struct integrator_base {
    integrator_base(
        const blas_ops& _blas,
        const std::map<std::string, integration_method*>& _integration_methods)
        : blas(_blas), integration_methods(_integration_methods)
    {
    }

    const blas_ops blas;
    const std::map<std::string, integration_method*> integration_methods;
};

#endif
