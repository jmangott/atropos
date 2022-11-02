#ifndef GRID_CLASS_HPP
#define GRID_CLASS_HPP

#include <stdexcept>

#include <generic/storage.hpp>


// struct for storing the grid parameters in a compact form
template <Index m1, Index m2>
struct grid_info
{
    Index r;
    array<Index, m1> n1;
    array<Index, m2> n2;
    array<Index, m1> k1;
    array<Index, m2> k2;
    array<double, m1> lim1;
    array<double, m2> lim2;
    array<double, m1> h1;
    array<double, m2> h2;
    Index dx1, dx2;

    grid_info(Index _r, Index _n, Index _k) : r(_r)
    {
        // TODO: implement a function that enables reading in of custom values for n_xx1, n_xx2 and k_xx1, k_xx2
        // Use n to initialize the number of grid points n_xx1 and n_xx2 uniformly for all species
        for (auto &ele : n1)
            ele = _n;
        for (auto &ele : n2)
            ele = _n;

        // Use k to initialize the grid density k_xx1 and k_xx2 uniformly for all species
        for (auto &ele : k1)
            ele = _k;
        for (auto &ele : k2)
            ele = _k;

        // Calculate the grid limits and spacing from the number of grid points and the grid density
        for (Index i = 0; i < m1; i++)
        {
            if ((n1[i] - 1) % k1[i] != 0)
            {
                std::cerr << "ERROR: (n_xx1 - 1) must be a multiple of k_xx1!";
                std::abort();
            }
            else
            {
                lim1[i] = (n1[i] - 1) / k1[i];
                h1[i] = 1.0 / k1[i];
            }
        }

        for (Index i = 0; i < m2; i++)
        {
            if ((n2[i] - 1) % k2[i] != 0)
            {
                std::cerr << "ERROR: (n_xx2 - 1) must be a multiple of k_xx2!";
                std::abort();
            }
            else
            {
                lim2[i] = (n2[i] - 1) / k2[i];
                h2[i] = 1.0 / k2[i];
            }
        }

        // Calculate the total number of grid points dxx1_mult and dxx2_mult
        dx1 = 1;
        dx2 = 1;
        for (auto const ele : n1)
            dx1 *= ele;
        for (auto const ele : n2)
            dx2 *= ele;
    }
};


// Initialize all necessary objects for datum generation and the LR algorithm
void InitializeAuxiliaryObjects(multi_array<Index, 1> &n_xx1, multi_array<Index, 1> &n_xx2, multi_array<Index, 1> &k_xx1, multi_array<Index, 1> &k_xx2, multi_array<double, 1> &lim_xx1, multi_array<double, 1> &lim_xx2, multi_array<double, 1> &h_xx1, multi_array<double, 1> &h_xx2, Index &dxx1_mult, Index &dxx2_mult, Index m1, Index m2, Index n, Index k);

#endif