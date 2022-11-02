#include "grid_class.hpp"

void InitializeAuxiliaryObjects(multi_array<Index, 1> &n_xx1, multi_array<Index, 1> &n_xx2, multi_array<Index, 1> &k_xx1, multi_array<Index, 1> &k_xx2, multi_array<double, 1> &lim_xx1, multi_array<double, 1> &lim_xx2, multi_array<double, 1> &h_xx1, multi_array<double, 1> &h_xx2, Index &dxx1_mult, Index &dxx2_mult, Index m1, Index m2, Index n, Index k)
{
    // TODO: implement a function that enables reading in of custom values for n_xx1, n_xx2 and k_xx1, k_xx2
    // Use n to initialize the number of grid points n_xx1 and n_xx2 uniformly for all species
    for (auto &ele : n_xx1)
        ele = n;
    for (auto &ele : n_xx2)
        ele = n;

    // Use k to initialize the grid density k_xx1 and k_xx2 uniformly for all species
    for (auto &ele : k_xx1)
        ele = k;
    for (auto &ele : k_xx2)
        ele = k;

    // Calculate the grid limits and spacing from the number of grid points and the grid density
    for (Index ii = 0; ii < m1; ii++)
    {
        if ((n_xx1(ii) - 1) % k_xx1(ii) != 0)
        {
            std::cerr << "ERROR: (n_xx1 - 1) must be a multiple of k_xx1!";
            std::abort();
        }
        else
        {
            lim_xx1(ii) = (n_xx1(ii) - 1) / k_xx1(ii);
            h_xx1(ii) = 1.0 / k_xx1(ii);
        }
    }

    for (Index ii = 0; ii < m2; ii++)
    {
        if ((n_xx2(ii) - 1) % k_xx2(ii) != 0)
        {
            std::cerr << "ERROR: (n_xx2 - 1) must be a multiple of k_xx2!";
            std::abort();
        }
        else
        {
            lim_xx2(ii) = (n_xx2(ii) - 1) / k_xx2(ii);
            h_xx2(ii) = 1.0 / k_xx2(ii);
        }
    }

    // Calculate the total number of grid points dxx1_mult and dxx2_mult
    dxx1_mult = 1;
    dxx2_mult = 1;
    for (auto const ele : n_xx1)
        dxx1_mult *= ele;
    for (auto const ele : n_xx2)
        dxx2_mult *= ele;
}
