#include "grid_class.hpp"

using std::endl;

// Initialize grid limits, grid spacing and the total number of grid points
void grid_info::grid_common_init()
{
    d = m1 + m2;
    dx1 = 1;
    dx2 = 1;
    h1_mult = 1;
    h2_mult = 1;
    h_mult = 1;

    for (Index i = 0; i < m1; i++)
    {
        n(i) = n1(i);
        binsize(i) = binsize1(i);
        liml(i) = liml1(i);

        dx1 *= n1(i);
        h1_mult *= binsize1(i);
        h_mult *= binsize1(i);
    }

    for (Index i = 0; i < m2; i++)
    {
        n(i + m1) = n2(i);
        binsize(i + m1) = binsize2(i);
        liml(i + m1) = liml2(i);

        dx2 *= n2(i);
        h2_mult *= binsize2(i);
        h_mult *= binsize2(i);
    }
}

grid_info::grid_info(Index _m1, Index _m2, Index _r, Index _n, Index _binsize, vector<double> _liml1, vector<double> _liml2) : m1(_m1), m2(_m2), r(_r), n1({_m1}), n2({_m2}), n({_m1 + _m2}), binsize1({_m1}), binsize2({_m2}), binsize({_m1 + _m2}), liml1({_m1}), liml2({_m2}), liml({_m1 + _m2})
{
    // Use `n` and `binsize` to initialize the number of grid points and grid density uniformly for all species
    std::fill(n1.begin(), n1.end(), _n);
    std::fill(n2.begin(), n2.end(), _n);
    std::fill(binsize1.begin(), binsize1.end(), _binsize);
    std::fill(binsize2.begin(), binsize2.end(), _binsize);
    std::copy(_liml1.begin(), _liml1.end(), liml1.begin());
    std::copy(_liml2.begin(), _liml2.end(), liml2.begin());

    grid_common_init();
}

grid_info::grid_info(Index _m1, Index _m2, Index _r, multi_array<Index, 1> _n1, multi_array<Index, 1> _n2, multi_array<Index, 1> _binsize1, multi_array<Index, 1> _binsize2, multi_array<double, 1> _liml1, multi_array<double, 1> _liml2) : m1(_m1), m2(_m2), r(_r), n1(_n1), n2(_n2), n({_m1 + _m2}), binsize1(_binsize1), binsize2(_binsize2), binsize({_m1 + _m2}), liml1(_liml1), liml2(_liml2), liml({_m1 + _m2})
{
    grid_common_init();
}

grid_info::grid_info(Index _m1, Index _m2, Index _r, const Index _n1[], const Index _n2[], const Index _binsize1[], const Index _binsize2[], const double _liml1[], const double _liml2[]) : m1(_m1), m2(_m2), r(_r), n1({_m1}), n2({_m2}), n({_m1 + _m2}), binsize1({_m1}), binsize2({_m2}), binsize({_m1 + _m2}), liml1({_m1}), liml2({_m2}), liml({_m1 + _m2})
{
    std::copy(_n1, _n1 + _m1, n1.begin());
    std::copy(_n2, _n2 + _m2, n2.begin());
    std::copy(_binsize1, _binsize1 + _m1, binsize1.begin());
    std::copy(_binsize2, _binsize2 + _m2, binsize2.begin());
    std::copy(_liml1, _liml1 + _m1, liml1.begin());
    std::copy(_liml2, _liml2 + _m2, liml2.begin());

    grid_common_init();
}