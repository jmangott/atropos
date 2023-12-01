#ifndef SUBFLOWS_HPP
#define SUBFLOWS_HPP

// #include <chrono>
// #include <fstream>
// #include <functional>
// #include <iostream>
// #include <sstream>
// #include <stdexcept>
// #include <string>
// #include <vector>

#include <generic/matrix.hpp>
// #include <generic/storage.hpp>
// #include <lr/coefficients.hpp>
// #include <lr/lr.hpp>

#include "integrators.hpp"
#include "tree_class.hpp"

// TODO: Verify transposed quantities -> should be ok
// TODO: check if input from netcdf file is transposed in (r_out, r_out)-dimensions or transpose it in Python interface
template <Index id>
void SubflowPhi(cme_internal_node * const node, const blas_ops &blas, const double tau);

void SubflowPsi(cme_internal_node * const node, const blas_ops &blas, const double tau);

void SubflowK(cme_external_node * const node, const blas_ops &blas, const double tau);

void SubflowS(cme_node* const node, const blas_ops &blas, const double tau);

#endif