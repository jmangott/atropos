#ifndef GRID_CLASS_HPP
#define GRID_CLASS_HPP

#include <stdexcept>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>


// struct for storing the grid parameters in a compact form
struct grid_info
{
    Index m1;
    Index m2;
    Index r;
    multi_array<Index, 1> n1;
    multi_array<Index, 1> n2;
    multi_array<Index, 1> n;
    multi_array<Index, 1> k1;
    multi_array<Index, 1> k2;
    multi_array<Index, 1> k;
    multi_array<double, 1> lim1;
    multi_array<double, 1> lim2;
    multi_array<double, 1> lim;
    multi_array<double, 1> h1;
    multi_array<double, 1> h2;
    multi_array<double, 1> h;
    Index dx1, dx2;
    double h1_mult, h2_mult;

    // Initialize grid limits, grid spacing and the total number of grid points
    void grid_common_init();

    grid_info(Index _m1, Index _m2, Index _r, Index _n, Index _k);

    grid_info(Index _m1, Index _m2, Index _r, multi_array<Index, 1> _n1, multi_array<Index, 1> _n2, multi_array<Index, 1> _k1, multi_array<Index, 1> _k2);

    grid_info(Index _m1, Index _m2, Index _r, std::vector<Index> _n1, std::vector<Index> _n2, std::vector<Index> _k1, std::vector<Index> _k2);
};


#endif