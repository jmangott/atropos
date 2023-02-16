#include "grid_class.hpp"

using std::endl;

// Initialize grid limits, grid spacing and the total number of grid points
void grid_info::grid_common_init()
{
    dx1 = 1;
    dx2 = 1;
    h1_mult = 1;
    h2_mult = 1;

    for (Index i = 0; i < m1; i++)
    {
        if ((n1(i) - 1) % k1(i) != 0)
        {
            std::cerr << "ERROR: (n1 - 1) must be a multiple of k1!" << endl;
            std::abort();
        }
        else
        {
            limr1(i) = (n1(i) - 1) / k1(i) + liml1(i);
            h1(i) = 1.0 / k1(i);

            n(i) = n1(i);
            k(i) = k1(i);
            liml(i) = liml1(i);
            limr(i) = limr1(i);
            h(i) = h1(i);

            dx1 *= n1(i);
            h1_mult *= h1(i);
        }
    }

    for (Index i = 0; i < m2; i++)
    {
        if ((n2(i) - 1) % k2(i) != 0)
        {
            std::cerr << "ERROR: (n2 - 1) must be a multiple of k2!" << endl;
            std::abort();
        }
        else
        {
            limr2(i) = (n2(i) - 1) / k2(i) + liml2(i);
            h2(i) = 1.0 / k2(i);

            n(i + m1) = n2(i);
            k(i + m1) = k2(i);
            liml(i + m1) = liml2(i);
            limr(i + m1) = limr2(i);
            h(i + m1) = h2(i);

            dx2 *= n2(i);
            h2_mult *= h2(i);
        }
    }
}

grid_info::grid_info(Index _m1, Index _m2, Index _r, Index _n, Index _k, vector<double> _liml1, vector<double> _liml2) : m1(_m1), m2(_m2), r(_r), n1({_m1}), n2({_m2}), n({_m1 + _m2}), k1({_m1}), k2({_m2}), k({_m1 + _m2}), liml1({_m1}), liml2({_m2}), liml({_m1 + _m2}), limr1({_m1}), limr2({_m2}), limr({_m1 + _m2}), h1({_m1}), h2({_m2}), h({_m1 + _m2})
{
    // Use n and k to initialize the number of grid points and grid density uniformly for all species
    std::fill(n1.begin(), n1.end(), _n);
    std::fill(n2.begin(), n2.end(), _n);
    std::fill(k1.begin(), k1.end(), _k);
    std::fill(k2.begin(), k2.end(), _k);
    std::copy(_liml1.begin(), _liml1.end(), liml1.begin());
    std::copy(_liml2.begin(), _liml2.end(), liml2.begin());

    grid_common_init();
}

grid_info::grid_info(Index _m1, Index _m2, Index _r, multi_array<Index, 1> _n1, multi_array<Index, 1> _n2, multi_array<Index, 1> _k1, multi_array<Index, 1> _k2, multi_array<double, 1> _liml1, multi_array<double, 1> _liml2) : m1(_m1), m2(_m2), r(_r), n1(_n1), n2(_n2), n({_m1 + _m2}), k1(_k1), k2(_k2), k({_m1 + _m2}), liml1(_liml1), liml2(_liml2), liml({_m1 + _m2}), limr1({_m1}), limr2({_m2}), limr({_m1 + _m2}), h1({_m1}), h2({_m2}), h({_m1 + _m2})
{
    grid_common_init();
}

grid_info::grid_info(Index _m1, Index _m2, Index _r, vector<Index> _n1, vector<Index> _n2, vector<Index> _k1, vector<Index> _k2, vector<double> _liml1, vector<double> _liml2) : m1(_m1), m2(_m2), r(_r), n1({_m1}), n2({_m2}), n({_m1 + _m2}), k1({_m1}), k2({_m2}), k({_m1 + _m2}), liml1({_m1}), liml2({_m2}), liml({_m1 + _m2}), limr1({_m1}), limr2({_m2}), limr({_m1 + _m2}), h1({_m1}), h2({_m2}), h({_m1 + m2})
{
    std::copy(_n1.begin(), _n1.end(), n1.begin());
    std::copy(_n2.begin(), _n2.end(), n2.begin());
    std::copy(_k1.begin(), _k1.end(), k1.begin());
    std::copy(_k2.begin(), _k2.end(), k2.begin());
    std::copy(_liml1.begin(), _liml1.end(), liml1.begin());
    std::copy(_liml2.begin(), _liml2.end(), liml2.begin());

    grid_common_init();
}