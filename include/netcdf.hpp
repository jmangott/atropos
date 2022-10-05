#pragma once
#include <netcdf.h>

#include <Eigen/Dense>
using Eigen::MatrixXd;

#include "common.hpp"

#define ERROR_NETCDF(e)                                                         \
    {                                                                           \
        printf("NetCDF Error %s:%i: %s\n", __FILE__, __LINE__, nc_strerror(e)); \
        exit(1);                                                                \
    }

void write_nc(std::string fn, const MatrixXd &U, const MatrixXd &S, const MatrixXd &V, vector<std::string> names)
{
    int retval, ncid;

    if ((retval = nc_create(fn.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid)))
        ERROR_NETCDF(retval);

    int id_d, id_r, id_n_part1, id_n_part2;
    if ((retval = nc_def_dim(ncid, "d", names.size(), &id_d)))
        ERROR_NETCDF(retval);
    if ((retval = nc_def_dim(ncid, "r", U.cols(), &id_r)))
        ERROR_NETCDF(retval);
    if ((retval = nc_def_dim(ncid, "n1", U.rows(), &id_n_part1)))
        ERROR_NETCDF(retval);
    if ((retval = nc_def_dim(ncid, "n2", V.rows(), &id_n_part2)))
        ERROR_NETCDF(retval);

    int varid_U, varid_V, varid_S, varid_names;
    // U
    int dimids_U[2] = {id_r, id_n_part1};
    if ((retval = nc_def_var(ncid, "U", NC_DOUBLE, 2, dimids_U, &varid_U)))
        ERROR_NETCDF(retval);
    // V
    int dimids_V[2] = {id_r, id_n_part2};
    if ((retval = nc_def_var(ncid, "V", NC_DOUBLE, 2, dimids_V, &varid_V)))
        ERROR_NETCDF(retval);
    // S
    int dimids_S[2] = {id_r, id_r};
    if ((retval = nc_def_var(ncid, "S", NC_DOUBLE, 2, dimids_S, &varid_S)))
        ERROR_NETCDF(retval);
    // names
    if ((retval = nc_def_var(ncid, "names", NC_STRING, 1, &id_d, &varid_names)))
        ERROR_NETCDF(retval);

    if ((retval = nc_enddef(ncid)))
        ERROR_NETCDF(retval);

    if ((retval = nc_put_var_double(ncid, varid_U, U.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_put_var_double(ncid, varid_V, V.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_put_var_double(ncid, varid_S, S.data())))
        ERROR_NETCDF(retval);
    size_t count = 1;
    for (size_t i = 0; i < names.size(); i++)
    {
        const char *name = names[i].c_str();
        if ((retval = nc_put_vara_string(ncid, varid_names, &i, &count, &name)))
            ERROR_NETCDF(retval);
    }

    if ((retval = nc_close(ncid)))
        ERROR_NETCDF(retval);
}

void read_nc(std::string fn, MatrixXd &U, MatrixXd &S, MatrixXd &V, vector<std::string> &names)
{
    int retval, ncid;

    if ((retval = nc_open(fn.c_str(), NC_NOWRITE, &ncid)))
        ERROR_NETCDF(retval);

    // read dimensions
    int id_d, id_r, id_n_part1, id_n_part2;
    if ((retval = nc_inq_dimid(ncid, "d", &id_d)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_dimid(ncid, "r", &id_r)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_dimid(ncid, "n1", &id_n_part1)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_dimid(ncid, "n2", &id_n_part2)))
        ERROR_NETCDF(retval);

    size_t d, r, n1, n2;
    char tmp[NC_MAX_NAME + 1];
    if ((retval = nc_inq_dim(ncid, id_d, tmp, &d)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_dim(ncid, id_r, tmp, &r)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_dim(ncid, id_n_part1, tmp, &n1)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_dim(ncid, id_n_part2, tmp, &n2)))
        ERROR_NETCDF(retval);

    U.resize(n1, r);
    V.resize(n2, r);
    S.resize(r, r);

    // read variables
    int id_U, id_V, id_S, id_names;
    if ((retval = nc_inq_varid(ncid, "U", &id_U)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_varid(ncid, "V", &id_V)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_varid(ncid, "S", &id_S)))
        ERROR_NETCDF(retval);
    if ((retval = nc_inq_varid(ncid, "names", &id_names)))
        ERROR_NETCDF(retval);

    if ((retval = nc_get_var_double(ncid, id_U, U.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_double(ncid, id_V, V.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_get_var_double(ncid, id_S, S.data())))
        ERROR_NETCDF(retval);

    // read names
    vector<char *> ptr_names;
    for (ind i = 0; i < (ind)d; i++)
        ptr_names.push_back((char *)malloc(sizeof(char) * (NC_MAX_NAME + 1)));
    nc_get_var_string(ncid, id_names, &ptr_names[0]);

    for (ind i = 0; i < (ind)d; i++)
    {
        names.push_back(ptr_names[i]);
        free(ptr_names[i]);
    }

    if ((retval = nc_close(ncid)))
        ERROR_NETCDF(retval);
}