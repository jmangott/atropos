#include "io_functions.hpp"

using std::ifstream;
using std::string;
using std::stringstream;


void WriteNC(string fn, const lr2<double> &lr_sol, vector<string> names, grid_info grid, double *t, double *dt)
{
    int retval, ncid;

    multi_array<long long, 1> n1(grid.n1.shape()), n2(grid.n2.shape());
    multi_array<long long, 1> binsize(grid.binsize.shape());

    for (Index i = 0; i < grid.m1; i++)
        n1(i) = (long long) grid.n1(i);

    for (Index i = 0; i < grid.m2; i++)
        n2(i) = (long long) grid.n2(i);

    for (Index i = 0; i < grid.d; i++)
        binsize(i) = (long long) grid.binsize(i);

    if ((retval = nc_create(fn.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid)))
        NETCDF_ERROR(retval);

    // Dimensions
    int id_r, id_d, id_m1, id_m2, id_dx1, id_dx2;
    if ((retval = nc_def_dim(ncid, "r", grid.r, &id_r)))
        NETCDF_ERROR(retval);
    if ((retval = nc_def_dim(ncid, "d", grid.d, &id_d)))
        NETCDF_ERROR(retval);
    if ((retval = nc_def_dim(ncid, "m1", grid.m1, &id_m1)))
        NETCDF_ERROR(retval);
    if ((retval = nc_def_dim(ncid, "m2", grid.m2, &id_m2)))
        NETCDF_ERROR(retval);
    if ((retval = nc_def_dim(ncid, "dx1", grid.dx1, &id_dx1)))
        NETCDF_ERROR(retval);
    if ((retval = nc_def_dim(ncid, "dx2", grid.dx2, &id_dx2)))
        NETCDF_ERROR(retval);

    // Variables
    int varid_X, varid_V, varid_S, varid_names, varid_n1, varid_n2, varid_binsize, varid_liml, varid_t, varid_dt;
    // X
    int dimids_X[2] = {id_r, id_dx1};
    if ((retval = nc_def_var(ncid, "X", NC_DOUBLE, 2, dimids_X, &varid_X)))
        NETCDF_ERROR(retval);
    // V
    int dimids_V[2] = {id_r, id_dx2};
    if ((retval = nc_def_var(ncid, "V", NC_DOUBLE, 2, dimids_V, &varid_V)))
        NETCDF_ERROR(retval);
    // S
    int dimids_S[2] = {id_r, id_r};
    if ((retval = nc_def_var(ncid, "S", NC_DOUBLE, 2, dimids_S, &varid_S)))
        NETCDF_ERROR(retval);
    // names
    if ((retval = nc_def_var(ncid, "names", NC_STRING, 1, &id_d, &varid_names)))
        NETCDF_ERROR(retval);
    // n1
    if ((retval = nc_def_var(ncid, "n1", NC_INT64, 1, &id_m1, &varid_n1)))
        NETCDF_ERROR(retval);
    // n2
    if ((retval = nc_def_var(ncid, "n2", NC_INT64, 1, &id_m2, &varid_n2)))
        NETCDF_ERROR(retval);
    // binsize
    if ((retval = nc_def_var(ncid, "binsize", NC_INT64, 1, &id_d, &varid_binsize)))
        NETCDF_ERROR(retval);
    // liml
    if ((retval = nc_def_var(ncid, "liml", NC_DOUBLE, 1, &id_d, &varid_liml)))
        NETCDF_ERROR(retval);
    // t
    if ((retval = nc_def_var(ncid, "t", NC_DOUBLE, 0, 0, &varid_t)))
        NETCDF_ERROR(retval);
    // dt
    if ((retval = nc_def_var(ncid, "dt", NC_DOUBLE, 0, 0, &varid_dt)))
        NETCDF_ERROR(retval);

    if ((retval = nc_enddef(ncid)))
        NETCDF_ERROR(retval);

    if ((retval = nc_put_var_double(ncid, varid_X, lr_sol.X.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_put_var_double(ncid, varid_V, lr_sol.V.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_put_var_double(ncid, varid_S, lr_sol.S.data())))
        NETCDF_ERROR(retval);
    size_t count = 1;
    for (size_t i = 0; i < names.size(); i++)
    {
        const char *name = names[i].c_str();
        if ((retval = nc_put_vara_string(ncid, varid_names, &i, &count, &name)))
            NETCDF_ERROR(retval);
    }
    if ((retval = nc_put_var_longlong(ncid, varid_n1, n1.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_put_var_longlong(ncid, varid_n2, n2.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_put_var_longlong(ncid, varid_binsize, binsize.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_put_var_double(ncid, varid_liml, grid.liml.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_put_var_double(ncid, varid_t, t)))
        NETCDF_ERROR(retval);
    if ((retval = nc_put_var_double(ncid, varid_dt, dt)))
        NETCDF_ERROR(retval);

    if ((retval = nc_close(ncid)))
        NETCDF_ERROR(retval);
}


void ReadNC(string fn, multi_array<double, 2> &xx1, multi_array<double, 2> &xx2, multi_array<double, 2> &ss, Index &n_basisfunctions)
{
    int retval, ncid;

    if ((retval = nc_open(fn.c_str(), NC_NOWRITE, &ncid)))
        NETCDF_ERROR(retval);

    // read dimensions
    int id_r, id_dx1, id_dx2, id_n_basisfunctions;
    if ((retval = nc_inq_dimid(ncid, "r", &id_r)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_dimid(ncid, "dx1", &id_dx1)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_dimid(ncid, "dx2", &id_dx2)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_dimid(ncid, "n_basisfunctions", &id_n_basisfunctions)))
        NETCDF_ERROR(retval);

    size_t r_t, dx1_t, dx2_t, n_basisfunctions_t;
    char tmp[NC_MAX_NAME + 1];
    if ((retval = nc_inq_dim(ncid, id_r, tmp, &r_t)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_dim(ncid, id_dx1, tmp, &dx1_t)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_dim(ncid, id_dx2, tmp, &dx2_t)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_dim(ncid, id_n_basisfunctions, tmp, &n_basisfunctions_t)))
        NETCDF_ERROR(retval);

    Index r = (Index) r_t;
    Index dx1 = (Index) dx1_t;
    Index dx2 = (Index) dx2_t;
    n_basisfunctions = (Index) n_basisfunctions_t;

    xx1.resize({dx1, n_basisfunctions});
    xx2.resize({dx2, n_basisfunctions});
    ss.resize({r, r});

    // read variables
    int id_xx1, id_xx2, id_ss;
    if ((retval = nc_inq_varid(ncid, "xx1", &id_xx1)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_varid(ncid, "xx2", &id_xx2)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_varid(ncid, "ss", &id_ss)))
        NETCDF_ERROR(retval);

    if ((retval = nc_get_var_double(ncid, id_xx1, xx1.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_double(ncid, id_xx2, xx2.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_double(ncid, id_ss, ss.data())))
        NETCDF_ERROR(retval);

    if ((retval = nc_close(ncid)))
        NETCDF_ERROR(retval);
}