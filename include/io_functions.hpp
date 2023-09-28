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


grid_parms ReadGridParms(int ncid);


template<class T>
void ReadInternalNode(int ncid, std::string id, node<T>* parent_node, node<T>* &node)
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
    grid_parms grid = ReadGridParms(ncid);

    // create node
    if (id == "0" || id == "1")
    {
        node = new cme_internal_node<T>(dynamic_cast<cme_root<T> *>(parent_node), grid, r);
    }
    else
    {
        node = new cme_internal_node<T>(dynamic_cast<cme_internal_node<T> *>(parent_node), grid, r);
    }
    
    // read Q and S
    int id_Q, id_S;
    if ((retval = nc_inq_varid(ncid, "Q", &id_Q)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_varid(ncid, "S", &id_S)))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_double(ncid, id_Q, dynamic_cast<cme_internal_node<T>*>(node)->Q.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_double(ncid, id_S, dynamic_cast<cme_internal_node<T> *>(node)->S.data())))
        ERROR_NETCDF(retval);
}

template <class T>
void ReadExternalNode(int ncid, std::string id, node<T> *parent_node, node<T>* &node)
{
    int retval;

    // read grid
    grid_parms grid = ReadGridParms(ncid);

    // create node
    if (id == "0" || id == "1")
    {
        node = new cme_external_node<T>(dynamic_cast<cme_root<T>*>(parent_node), grid);
    }
    else
    {
        node = new cme_external_node<T>(dynamic_cast<cme_internal_node<T>*>(parent_node), grid);
    }

    // read X
    int id_X;
    if ((retval = nc_inq_varid(ncid, "X", &id_X)))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_double(ncid, id_X, dynamic_cast<cme_external_node<T>*>(node)->X.data())))
        ERROR_NETCDF(retval);
}


template<class T>
void __ReadHierarchicalNC(int ncid, std::string id, node<T>* parent_node, node<T>* &node)
{
    int retval, grp_ncid;
    if ((retval = nc_inq_ncid(ncid, (id + "0").c_str(), &grp_ncid)))
    {
        ReadExternalNode(ncid, id, parent_node, node);
    }
    else
    {
        ReadInternalNode(ncid, id, parent_node, node);
        __ReadHierarchicalNC(grp_ncid, id + "0", node, dynamic_cast<internal_node<T> *>(node)->left);
    }

    if ((retval = nc_inq_ncid(ncid, (id + "1").c_str(), &grp_ncid)))
    {
        ReadExternalNode(ncid, id, parent_node, node);
    }
    else
    {
        ReadInternalNode(ncid, id, parent_node, node);
        __ReadHierarchicalNC(grp_ncid, id + "1", node, node->right);
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
    grid_parms grid = ReadGridParms(ncid);

    // create node
    tree.root_node = new cme_root<T>(grid, r);

    // read S
    int id_S;
    if ((retval = nc_inq_varid(ncid, "S", &id_S)))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_double(ncid, id_S, tree.root_node->S.data())))
        ERROR_NETCDF(retval);

    int grp_ncid;

    if ((retval = nc_inq_ncid(ncid, "0", &grp_ncid)))
        ERROR_NETCDF(retval);
    __ReadHierarchicalNC(grp_ncid, "0", tree.root_node, tree.root_node->left);

    if ((retval = nc_inq_ncid(ncid, "1", &grp_ncid)))
        ERROR_NETCDF(retval);
    __ReadHierarchicalNC(grp_ncid, "1", tree.root_node, tree.root_node->right);

    if ((retval = nc_close(ncid)))
        ERROR_NETCDF(retval);
    return;
}

#endif