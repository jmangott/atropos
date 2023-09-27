#ifndef IO_FUNCTIONS_HPP
#define IO_FUNCTIONS_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <netcdf.h>

#include <generic/netcdf.hpp>
#include <generic/storage.hpp>

#include <lr/lr.hpp>

#include "coeff_class.hpp"
#include "grid_class.hpp"
#include "reaction_class.hpp"
#include "tree_class.hpp"

#define ERROR_NETCDF(e)                                                         \
    {                                                                           \
        printf("NetCDF Error %s:%i: %s\n", __FILE__, __LINE__, nc_strerror(e)); \
        exit(1);                                                                \
    }

// Write `lr_sol` to a .netCDF file
void WriteNC(std::string fn, const lr2<double> &lr_sol, vector<string> names, grid_info grid, double *t, double *dt);


// Read a .netCDF file to multi_arrays `xx1`, `xx2` and `ss`
void ReadNC(std::string fn, multi_array<double, 2> &xx1, multi_array<double, 2> &xx2, multi_array<double, 2> &ss, Index &n_basisfunctions);


grid_parms ReadGridParms(int ncid)
{
    int retval;

    // read dimensions
    int id_d, id_mu;
    if ((retval = nc_inq_dimid(ncid, "d", &id_d)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_dimid(ncid, "mu", &id_mu)))
        ERROR_NETCDF(retval);

    size_t d_t, mu_t;
    char tmp[NC_MAX_NAME + 1];
    if ((retval = nc_inq_dim(ncid, id_d, tmp, &d_t)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_dim(ncid, id_mu, tmp, &mu_t)))
        ERROR_NETCDF(retval);

    Index d = (Index)d_t;
    Index mu = (Index)mu_t;

    // read variables
    int id_n, id_binsize, id_liml, id_dep, id_S;

    if ((retval = nc_inq_varid(ncid, "n", &id_n)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_varid(ncid, "binsize", &id_binsize)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_varid(ncid, "liml", &id_liml)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_varid(ncid, "dep", &id_dep)))
        ERROR_NETCDF(retval);

    grid_parms grid(d, mu);
    multi_array<signed char, 2> dep_int({mu, d});

    if ((retval = nc_get_var_long(ncid, id_n, grid.n.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_long(ncid, id_binsize, grid.binsize.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_double(ncid, id_liml, grid.liml.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_schar(ncid, id_dep, dep_int.data())))
        ERROR_NETCDF(retval);

    std::copy(dep_int.begin(), dep_int.end(), grid.dep.begin());

    return grid;
}

template<class T>
void ReadInternalNode(int ncid, internal_node<T>* node)
{
    int retval;

    // read rank
    int id_r;
    if ((retval = nc_inq_dimid(ncid, "r", &id_r)))
        ERROR_NETCDF(retval);

    size_t r_t;
    char tmp[NC_MAX_NAME + 1];
    if ((retval = nc_inq_dim(ncid, id_r, tmp, &r_t)))
        ERROR_NETCDF(retval);

    Index r = (Index) r_t;

    // read grid
    grid_info grid = ReadGridParms(ncid);

    // create node
    tree.root_node = new cme_internal_node<T>(grid, r);
    
    // read Q and S
    int id_Q, id_S;
    if ((retval = nc_inq_varid(ncid, "Q", &id_Q)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_varid(ncid, "S", &id_S)))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_double(ncid, id_Q, tree.root_node->Q.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_double(ncid, id_S, tree.root_node->S.data())))
        ERROR_NETCDF(retval);
}


template<class T>
void ReadExternalNode(int ncid, external_node<T>* node)
{
    int retval;

    // read grid
    grid_info grid = ReadGridParms(ncid);

    // create node
    tree.root_node = new cme_external_node<T>(grid);

    // read X
    int id_X;
    if ((retval = nc_inq_varid(ncid, "X", &id_X)))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_double(ncid, id_X, tree.root_node->X.data())))
        ERROR_NETCDF(retval);
}


template<class T>
void __ReadHierarchicalNC(int ncid, std::string id, node<T>* node)
{
    int retval, grp_ncid;
    if (retval = nc_inq_grp_ncid(ncid, (id + "0").c_str(), &grp_ncid))
    {
        ReadInternalNode(node);
        __ReadHierarchicalNC(grp_ncid, id + "0", node->left);
    }
    else
    {
        ReadExternalNode(node);
    }

    if (retval = nc_inq_grp_ncid(ncid, (id + "1").c_str(), &grp_ncid))
    {
        ReadInternalNode(node);
        __ReadHierarchicalNC(grp_ncid, id + "1", node->right);
    }
    else
    {
        ReadExternalNode(node);
    }
}


// Read a .netCDF file with a hierarchical tree structure
template<class T>
void ReadHierarchicalNC(std::string fn, lr_tree<T> &tree)
{
    int retval, ncid;

    if ((retval = nc_open(fn.c_str(), NC_NOWRITE, &ncid)))
        ERROR_NETCDF(retval);

    // read rank
    int id_r;
    if ((retval = nc_inq_dimid(ncid, "r", &id_r)))
        ERROR_NETCDF(retval);

    size_t r_t;
    char tmp[NC_MAX_NAME + 1];
    if ((retval = nc_inq_dim(ncid, id_r, tmp, &r_t)))
        ERROR_NETCDF(retval);

    Index r = (Index)r_t;

    // read grid
    grid_info grid = ReadGridParms(ncid);

    // create node
    tree.root_node = new cme_root<T>(grid, r);

    // read S
    int id_S;
    if ((retval = nc_inq_varid(ncid, "S", &id_S)))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_double(ncid, id_S, tree.root_node->S.data())))
        ERROR_NETCDF(retval);

    int grp_ncid;

    if (retval = nc_inq_grp_ncid(ncid, (id + "0").c_str(), &grp_ncid))
        ERROR_NETCDF(retval);
    __ReadHierarchicalNC(tree.root_node, grp_ncid);

    if (retval = nc_inq_grp_ncid(ncid, (id + "1").c_str(), &grp_ncid))
        ERROR_NETCDF(retval);
    __ReadHierarchicalNC(tree.root_node, grp_ncid);

    if ((retval = nc_close(ncid)))
        ERROR_NETCDF(retval);
    return;
}

#endif