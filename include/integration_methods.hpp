#ifndef INTEGRATION_METHODS_HPP
#define INTEGRATION_METHODS_HPP

#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/IterativeSolvers>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "timer_class.hpp"

struct integration_method
{
    integration_method() = default;
    virtual ~integration_method() = default;
    virtual void integrate(multi_array<double, 2> &arr, const std::function<multi_array<double, 2>(const multi_array<double, 2> &)> &rhs, const double tau) const = 0;
};

struct explicit_euler : integration_method
{
    explicit_euler() = default;

    void integrate(multi_array<double, 2> &arr, const std::function<multi_array<double, 2>(const multi_array<double, 2> &)> &rhs, const double tau) const override
    {
        arr += rhs(arr) * tau;
    }
};

class matrix_free;
using Eigen::SparseMatrix;

namespace Eigen
{
    namespace internal
    {
        template <>
        struct traits<matrix_free> : public Eigen::internal::traits<Eigen::SparseMatrix<double>>
        {
        };
    }
}

class matrix_free : public Eigen::EigenBase<matrix_free>
{
public:
    typedef double Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;
    enum
    {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    Eigen::Index rows() const { return prod(arr_shape); }
    Eigen::Index cols() const { return prod(arr_shape); }

    template <typename Rhs>
    Eigen::Product<matrix_free, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs> &x) const
    {
        return Eigen::Product<matrix_free, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
    }

    // Custom API
    matrix_free(const std::array<Index, 2> _arr_shape, const std::function<multi_array<double, 2>(const multi_array<double, 2>&)> &_eqn, const double _tau)
        : arr_shape(_arr_shape)
        , eqn(_eqn)
        , tau(_tau)
        {}

    std::array<Index, 2> arr_shape;
    const std::function<multi_array<double, 2>(multi_array<double, 2>)> eqn;
    double tau;
};

namespace Eigen
{
    namespace internal
    {
        template <typename Rhs>
        struct generic_product_impl<matrix_free, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
            : generic_product_impl_base<matrix_free, Rhs, generic_product_impl<matrix_free, Rhs>>
        {
            typedef typename Product<matrix_free, Rhs>::Scalar Scalar;
            template <typename Dest>
            static void scaleAndAddTo(Dest &dst, const matrix_free &lhs, const Rhs &rhs, const Scalar &scale)
            {
                // From the Eigen documentation:
                // This method should implement "dst += alpha * lhs * rhs" inplace,
                // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
                assert(scale == Scalar(1) && "scaling is not implemented");

                double tau = lhs.tau;

                multi_array<double, 2> arr(lhs.arr_shape);
                Eigen::Map<Rhs>(arr.data(), rhs.size()) = rhs;

                arr -= lhs.eqn(arr) * tau;
                dst = Eigen::Map<Dest>(arr.data(), prod(arr.shape()));
            }
        };
    }
}

struct implicit_euler : integration_method
{
    implicit_euler() = default;

    void integrate(multi_array<double, 2> &arr, const std::function<multi_array<double, 2>(const multi_array<double, 2> &)> &rhs, const double tau) const override
    {
        matrix_free A(arr.shape(), rhs, tau);
        Eigen::GMRES<matrix_free, Eigen::IdentityPreconditioner> gmres;
        Eigen::VectorXd x, b;
        b = Eigen::Map<Eigen::VectorXd>(arr.data(), prod(arr.shape()));

        gmres.compute(A);
        x = gmres.solve(b);

        Eigen::Map<Eigen::VectorXd>(arr.data(), x.size()) = x;
    }
};

#endif