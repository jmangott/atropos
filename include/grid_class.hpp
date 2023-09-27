#ifndef GRID_CLASS_HPP
#define GRID_CLASS_HPP

#include <stdexcept>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "reaction_class.hpp"

// struct for storing the grid parameters in a compact form
struct grid_info
{
    Index d;
    Index m1;
    Index m2;
    Index r;
    multi_array<Index, 1> n1;
    multi_array<Index, 1> n2;
    multi_array<Index, 1> n;
    multi_array<Index, 1> binsize1;
    multi_array<Index, 1> binsize2;
    multi_array<Index, 1> binsize;
    multi_array<double, 1> liml1;
    multi_array<double, 1> liml2;
    multi_array<double, 1> liml;
    Index dx1, dx2;
    Index h1_mult, h2_mult, h_mult;

    // Initialize grid limits, grid spacing and the total number of grid points
    void grid_common_init();

    grid_info(Index _m1, Index _m2, Index _r, Index _n, Index _binsize, std::vector<double> _liml1, std::vector<double> _liml2);

    grid_info(Index _m1, Index _m2, Index _r, multi_array<Index, 1> _n1, multi_array<Index, 1> _n2, multi_array<Index, 1> _binsize1, multi_array<Index, 1> _binsize2, multi_array<double, 1> _liml1, multi_array<double, 1> _liml2);

    grid_info(Index _m1, Index _m2, Index _r, const Index _n1[], const Index _n2[], const Index _binsize1[], const Index _binsize2[], const double _liml1[], const double _liml2[]);
};


// struct for storing the grid parameters of one node of the hierarchical problem in a compact form

struct grid_parms
{
    Index d;
    Index mu;
    std::vector<Index> n;
    std::vector<Index> binsize;
    std::vector<double> liml;
    multi_array<bool, 2> dep;

    grid_parms(Index _d, Index _mu) : d(_d), mu(_mu), n(_d), binsize(_d), liml(_d), dep({_mu, _d}) {};

    grid_parms(std::vector<Index> _n, std::vector<Index> _binsize, std::vector<double> _liml, multi_array<bool, 2> _dep) : d(_n.size()), mu(_dep.shape()[0]), n(_n), binsize(_binsize), liml(_liml), dep(_dep) {};

    Index dx() const
    {
        Index result = 1;
        for (const auto &ele : n)
        {
            result *= ele;
        }
        return result;
    }

    Index h_mult() const
    {
        Index result = 1;
        for (const auto &ele : binsize)
        {
            result *= ele;
        }
        return result;
    }

    Index dx_dep(Index mu) const
    {
        Index result = 1;
        for (Index i = 0; i < d; i++)
        {
            if (dep(mu, i) == true)
            {
                result *= n[i];
            }
        }
        return result;
    }
};

#endif