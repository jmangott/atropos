#ifndef GRID_CLASS_HPP
#define GRID_CLASS_HPP

#include <stdexcept>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>


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

    grid_info(Index _m1, Index _m2, Index _r, Index _n, Index _k, std::vector<double> _liml1, std::vector<double> _liml2);

    grid_info(Index _m1, Index _m2, Index _r, multi_array<Index, 1> _n1, multi_array<Index, 1> _n2, multi_array<Index, 1> _k1, multi_array<Index, 1> _k2, multi_array<double, 1> _liml1, multi_array<double, 1> _liml2);

    grid_info(Index _m1, Index _m2, Index _r, std::vector<Index> _n1, std::vector<Index> _n2, std::vector<Index> _k1, std::vector<Index> _k2, std::vector<double> _liml1, std::vector<double> _liml2);
};


#endif