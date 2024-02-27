#ifndef INTEGRATION_METHODS_HPP
#define INTEGRATION_METHODS_HPP

#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "matrix_free.hpp"
#include "timer_class.hpp"
//TODO: use a function pointer instead of std::function for performance
struct integration_method
{
    integration_method() = default;
    virtual ~integration_method() = default;
    virtual void integrate(multi_array<double, 2> &arr, const std::function<multi_array<double, 2>(const multi_array<double, 2> &)> &rhs, const double tau) const = 0;
    virtual std::string get_name() const = 0;
};

struct explicit_euler : integration_method
{
    explicit_euler(const unsigned int _substeps)
    : substeps(_substeps)
    {};

    void integrate(multi_array<double, 2> &arr, const std::function<multi_array<double, 2>(const multi_array<double, 2> &)> &rhs, const double tau) const override
    {
        double tau_substep = tau / substeps;
        for (auto i = 0U; i < substeps; ++i)
        {
            arr += rhs(arr) * tau_substep;
        }
    }

    std::string get_name() const override
    {
        return "explicit_euler (" + std::to_string(substeps) + " substeps)";
    }

    const unsigned int substeps;
};

struct implicit_euler : integration_method
{
    implicit_euler(const unsigned int _substeps)
    : substeps(_substeps)
    {};

    void integrate(multi_array<double, 2> &arr, const std::function<multi_array<double, 2>(const multi_array<double, 2> &)> &rhs, const double tau) const override
    {
        Eigen::GMRES<matrix_free, Eigen::IdentityPreconditioner> gmres;
        Eigen::VectorXd x, b;
        b = Eigen::Map<Eigen::VectorXd>(arr.data(), prod(arr.shape()));

        double tau_substep = tau / substeps;
        matrix_free A(arr.shape(), rhs, tau_substep);
        gmres.compute(A);

        for (auto i = 0U; i < substeps; ++i)
        {
            x = gmres.solve(b);
            b = x;
        }

        Eigen::Map<Eigen::VectorXd>(arr.data(), x.size()) = x;
    }

    std::string get_name() const override
    {
        return "implicit_euler (" + std::to_string(substeps) + " substeps)";
    }

    const unsigned int substeps;
};

struct crank_nicolson : integration_method
{
    crank_nicolson(const unsigned int _substeps)
    : substeps(_substeps)
    {};

    void integrate(multi_array<double, 2> &arr, const std::function<multi_array<double, 2>(const multi_array<double, 2> &)> &rhs, const double tau) const override
    {
        double tau_half = 0.5 * tau;

        Eigen::GMRES<matrix_free, Eigen::IdentityPreconditioner> gmres;
        Eigen::VectorXd x, b;
        multi_array<double, 2> b_arr(arr.shape());

        double tau_substep = tau_half / substeps;
        matrix_free A(arr.shape(), rhs, tau_substep);
        gmres.compute(A);

        for (auto i = 0U; i < substeps; ++i)
        {
            b_arr = arr;
            b_arr += rhs(arr) * tau_half;

            b = Eigen::Map<Eigen::VectorXd>(b_arr.data(), prod(b_arr.shape()));

            x = gmres.solve(b);
            
            Eigen::Map<Eigen::VectorXd>(arr.data(), x.size()) = x;
        }

    }

    std::string get_name() const override
    {
        return "crank_nicolson (" + std::to_string(substeps) + " substeps)";
    }

    const unsigned int substeps;
};

#endif