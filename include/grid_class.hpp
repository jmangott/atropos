#ifndef GRID_CLASS_HPP
#define GRID_CLASS_HPP

#include <stdexcept>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>


// struct for storing the grid parameters in a compact form
template<Index m1, Index m2>
struct grid_info
{
    Index r;
    multi_array<Index, 1> n1;
    multi_array<Index, 1> n2;
    multi_array<Index, 1> k1;
    multi_array<Index, 1> k2;
    multi_array<double, 1> lim1;
    multi_array<double, 1> lim2;
    multi_array<double, 1> h1;
    multi_array<double, 1> h2;
    Index dx1, dx2;


    // Initialize grid limits, grid spacing and the total number of grid points
    void common_init()
    {
        // Calculate the grid limits and spacing from the number of grid points and the grid density
        for (Index i = 0; i < m1; i++)
        {
            if ((n1(i) - 1) % k1(i) != 0)
            {
                std::cerr << "ERROR: (n1 - 1) must be a multiple of k1!" << std::endl;
                std::abort();
            }
            else
            {
                lim1(i) = (n1(i) - 1) / k1(i);
                h1(i) = 1.0 / k1(i);
            }
        }

        for (Index i = 0; i < m2; i++)
        {
            if ((n2(i) - 1) % k2(i) != 0)
            {
                std::cerr << "ERROR: (n2 - 1) must be a multiple of k2!" << std::endl;
                std::abort();
            }
            else
            {
                lim2(i) = (n2(i) - 1) / k2(i);
                h2(i) = 1.0 / k2(i);
            }
        }

        // Calculate the total number of grid points dxx1_mult and dxx2_mult
        dx1 = 1;
        dx2 = 1;
        for (Index i = 0; i < m1; i++)
            dx1 *= n1(i);
        for (Index i = 0; i < m2; i++)
            dx2 *= n2(i);
    }


    grid_info(Index _r, Index _n, Index _k) : r(_r), n1({m1}), n2({m2}), k1({m1}), k2({m2}), lim1({m1}), lim2({m2}), h1({m1}), h2({m2})
    {
        // Use n to initialize the number of grid points n_xx1 and n_xx2 uniformly for all species
        for (Index i = 0; i < m1; i++)
            n1(i) = _n;
        for (Index i = 0; i < m2; i++)
            n2(i) = _n;

        // Use k to initialize the grid density k_xx1 and k_xx2 uniformly for all species
        for (Index i = 0; i < m1; i++)
            k1(i) = _k;
        for (Index i = 0; i < m2; i++)
            k2(i) = _k;

        common_init();
    }


    grid_info(Index _r, multi_array<Index, 1> _n1, multi_array<Index, 1> _n2, multi_array<Index, 1> _k1, multi_array<Index, 1> _k2) : r(_r), n1(_n1), n2(_n2), k1(_k1), k2(_k2), lim1({m1}), lim2({m2}), h1({m1}), h2({m2})
    {
        common_init();
    }
};


#endif