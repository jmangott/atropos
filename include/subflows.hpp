#ifndef SUBFLOWS_HPP
#define SUBFLOWS_HPP

#include <generic/matrix.hpp>

#include "integration_methods.hpp"
#include "timer_class.hpp"
#include "tree_class.hpp"

void TTNIntegrator(cme_internal_node *node, const blas_ops &blas, const double tau, const integration_method &method);

template <Index id>
void SubflowPhi(cme_internal_node * const node, const blas_ops &blas, const double tau, const integration_method &method);

void SubflowPsi(cme_internal_node * const node, const blas_ops &blas, const double tau);

#endif