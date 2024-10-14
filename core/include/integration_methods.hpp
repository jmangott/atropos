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
#include <generic/timer.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "matrix_free.hpp"
// TODO: use a function pointer instead of std::function for performance
struct integration_method {
    integration_method() = default;
    virtual ~integration_method() = default;
    virtual void integrate(Ensign::multi_array<double, 2>& arr,
                           const std::function<Ensign::multi_array<double, 2>(
                               const Ensign::multi_array<double, 2>&)>& rhs,
                           const double tau) const = 0;
    virtual std::string get_name() const = 0;
};

struct explicit_euler : integration_method {
    explicit_euler(const unsigned int _substeps) : substeps(_substeps) {};

    void integrate(Ensign::multi_array<double, 2>& arr,
                   const std::function<Ensign::multi_array<double, 2>(
                       const Ensign::multi_array<double, 2>&)>& rhs,
                   const double tau) const override
    {
        double tau_substep = tau / substeps;
        for (auto i = 0U; i < substeps; ++i) {
            arr += rhs(arr) * tau_substep;
        }
    }

    std::string get_name() const override
    {
        return "explicit_euler (" + std::to_string(substeps) + " substeps)";
    }

    const unsigned int substeps;
};

struct rk4 : integration_method {
    rk4(const unsigned int _substeps) : substeps(_substeps) {};

    void integrate(Ensign::multi_array<double, 2>& arr,
                   const std::function<Ensign::multi_array<double, 2>(
                       const Ensign::multi_array<double, 2>&)>& rhs,
                   const double tau) const override
    {
        Ensign::multi_array<double, 2> k1(arr.shape());
        Ensign::multi_array<double, 2> k2(arr.shape());
        Ensign::multi_array<double, 2> k3(arr.shape());
        Ensign::multi_array<double, 2> k4(arr.shape());

        double tau_substep = tau / substeps;
        for (auto i = 0U; i < substeps; ++i) {
            k1 = rhs(arr);
            k2 = rhs(arr + k1 * (0.5 * tau_substep));
            k3 = rhs(arr + k2 * (0.5 * tau_substep));
            k4 = rhs(arr + k3 * tau_substep);

            arr += (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (tau_substep / 6.0);
        }
    }

    std::string get_name() const override
    {
        return "rk4 (" + std::to_string(substeps) + " substeps)";
    }

    const unsigned int substeps;
};

struct implicit_euler : integration_method {
    implicit_euler(const unsigned int _substeps) : substeps(_substeps) {};

    void integrate(Ensign::multi_array<double, 2>& arr,
                   const std::function<Ensign::multi_array<double, 2>(
                       const Ensign::multi_array<double, 2>&)>& rhs,
                   const double tau) const override
    {
        Eigen::GMRES<matrix_free, Eigen::IdentityPreconditioner> gmres;
        Eigen::VectorXd x, b;
        b = Eigen::Map<Eigen::VectorXd>(arr.data(), Ensign::prod(arr.shape()));

        double tau_substep = tau / substeps;
        matrix_free A(arr.shape(), rhs, tau_substep);
        gmres.compute(A);

        for (auto i = 0U; i < substeps; ++i) {
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

struct crank_nicolson : integration_method {
    crank_nicolson(const unsigned int _substeps) : substeps(_substeps) {};

    void integrate(Ensign::multi_array<double, 2>& arr,
                   const std::function<Ensign::multi_array<double, 2>(
                       const Ensign::multi_array<double, 2>&)>& rhs,
                   const double tau) const override
    {
        double tau_half = 0.5 * tau;

        Eigen::GMRES<matrix_free, Eigen::IdentityPreconditioner> gmres;
        Eigen::VectorXd x, b;
        Ensign::multi_array<double, 2> b_arr(arr.shape());

        double tau_substep = tau_half / substeps;
        matrix_free A(arr.shape(), rhs, tau_substep);
        gmres.compute(A);

        for (auto i = 0U; i < substeps; ++i) {
            b_arr = arr;
            b_arr += rhs(arr) * tau_half;

            b = Eigen::Map<Eigen::VectorXd>(b_arr.data(), Ensign::prod(b_arr.shape()));

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
