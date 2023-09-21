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

grid_parms::grid_parms(vector<Index> _n, vector<Index> _binsize, vector<double> _liml, vector<vector<Index>> _dep) : n(_n), binsize(_binsize), liml(_liml), d(_n.size()), dep(_dep), n_dep(_dep.size()), n_rem(_dep.size()), dx_dep(_dep.size()), dx_rem(_dep.size())
{
    dx = 1;
    h_mult = 1;

    for (Index i = 0; i < d; i++)
    {
        dx *= n[i];
        h_mult *= binsize[i];
    }
    for (size_t mu = 0; mu < _dep.size(); mu++)
    {
        dx_dep[mu] = 1;
        dx_rem[mu] = 1;
        n_rem[mu] = n;
        for (auto const &dep_ele : dep[mu])
        {
            n_dep[mu].push_back(n[dep_ele]);
            n_rem[mu][dep_ele] = 1;
        }
    }
}

// grid_parms::grid_parms(grid_parms _grid1, grid_parms _grid2) : d(_grid1.d + _grid2.d) 
// {
//     n = _grid1.n;
//     binsize = _grid1.binsize;
//     liml = _grid1.liml;

//     n.insert(n.end(), _grid2.n.begin(), _grid2.n.end());
//     binsize.insert(binsize.end(), _grid2.binsize.begin(), _grid2.binsize.end());
//     liml.insert(liml.end(), _grid2.liml.begin(), _grid2.liml.end());

//     dep = _grid1.dep;
//     n_dep = _grid1.n_dep;
//     n_rem = _grid1.n_rem;
//     dx_dep = _grid1.dx_dep;
//     dx_rem = _grid1.dx_rem;

//     for (size_t mu = 0; mu < _grid1.dep.size(); mu++)
//     {
//         dep[mu].insert(dep[mu].end(), _grid2.dep[mu].begin(), _grid2.dep[mu].end());
//         n_dep[mu].insert(n_dep[mu].end(), _grid2.n_dep[mu].begin(), _grid2.n_dep[mu].end());
//         n_rem[mu].insert(n_rem[mu].end(), _grid2.n_rem[mu].begin(), _grid2.n_rem[mu].end());
//         dx_dep[mu] *= _grid2.dx_dep[mu];
//         dx_rem[mu] *= _grid2.dx_rem[mu];
//     }
// }