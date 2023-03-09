#include "io_functions.hpp"

using std::ifstream;
using std::string;
using std::stringstream;


void ReadCSV(multi_array<double, 2> &output_array, string filename)
{
    Index i = 0;
    Index j = 0;
    string line;
    ifstream input_file(filename);
    stringstream ss_line;
    string element;

    while (getline(input_file, line))
    {
        ss_line.str(line);
        while (getline(ss_line, element, ','))
        {
            output_array(i, j) = stod(element);
            j++;
        }
        i++;
        j = 0;
        ss_line.clear();
    }
    input_file.close();
}


void WriteCSV(const multi_array<double, 2> &input_array, string filename)
{
    stringstream filename_type;
    filename_type << filename << ".csv";
    ofstream output_file(filename_type.str());
    Index n_rows, n_cols;
    n_rows = input_array.shape()[0];
    n_cols = input_array.shape()[1];
    for (Index i = 0; i < n_rows; i++)
    {
        for (Index j = 0; j < n_cols; j++)
        {
            output_file << input_array(i, j) << ',';
        }
        output_file << endl;
    }
    output_file.close();
}


// void WriteNetCDF(string fn, const lr2<double> &lr_sol, const mysys &reaction_system, grid_info grid, double *t, double *dt)
// {
//     nc_writer ncw(fn, {grid.r, grid.m1, grid.m2, grid.dx1, grid.dx2}, {"r", "m1", "m2", "dx1", "dx2"});
//     ncw.add_var("n1", {"m1"});
//     ncw.add_var("n2", {"m2"});
//     ncw.add_var("X", {"r", "dx1"});
//     ncw.add_var("S", {"r", "r"});
//     ncw.add_var("V", {"r", "dx2"});

//     ncw.start_write_mode();

//     vector<double> vec_r(1), vec_m1(1), vec_m2(1), vec_dx1(1), vec_dx2(1);
//     vec_r[0] = grid.r;
//     vec_m1[0] = grid.m1;
//     vec_m2[0] = grid.m2;
//     vec_dx1[0] = grid.dx1;
//     vec_dx2[0] = grid.dx2;

//     vector<double> vec_n1(grid.m1), vec_n2(grid.m2);
//     for (Index i = 0; i < grid.m1; i++)
//         vec_n1[i] = grid.n1(i);
//     for (Index i = 0; i < grid.m2; i++)
//         vec_n2[i] = grid.n2(i);

//     ncw.write("t", t);
//     ncw.write("t", dt);

//     // ncw.write("names", reaction_system.species_names.data());

//     // ncw.write("r", vec_r.data());
//     // ncw.write("m1", vec_m1.data());
//     // ncw.write("m2", vec_m2.data());
//     // ncw.write("dx1", vec_dx1.data());
//     // ncw.write("dx2", vec_dx2.data());
//     ncw.write("n1", vec_n1.data());
//     ncw.write("n2", vec_n2.data());
//     ncw.write("liml1", grid.liml1.data());
//     ncw.write("liml2", grid.liml2.data());
//     ncw.write("limr1", grid.limr1.data());
//     ncw.write("limr2", grid.limr2.data());

//     ncw.write("X", lr_sol.X.data());
//     ncw.write("S", lr_sol.S.data());
//     ncw.write("V", lr_sol.V.data());
// }

#define ERROR_NETCDF(e)                                                         \
    {                                                                           \
        printf("NetCDF Error %s:%i: %s\n", __FILE__, __LINE__, nc_strerror(e)); \
        exit(1);                                                                \
    }

void WriteNC(string fn, const lr2<double> &lr_sol, vector<string> names, grid_info grid, double *t, double *dt)
{
    int retval, ncid;

    multi_array<long long, 1> n1(grid.n1.shape()), n2(grid.n2.shape());

    for (Index i = 0; i < grid.m1; i++)
        n1(i) = (long long) grid.n1(i);

    for (Index i = 0; i < grid.m2; i++)
        n2(i) = (long long) grid.n2(i);

    if ((retval = nc_create(fn.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid)))
        ERROR_NETCDF(retval);

    // Dimensions
    int id_r, id_d, id_m1, id_m2, id_dx1, id_dx2, id_np;
    if ((retval = nc_def_dim(ncid, "r", grid.r, &id_r)))
        ERROR_NETCDF(retval);
    if ((retval = nc_def_dim(ncid, "d", grid.m1 + grid.m2, &id_d)))
        ERROR_NETCDF(retval);
    if ((retval = nc_def_dim(ncid, "m1", grid.m1, &id_m1)))
        ERROR_NETCDF(retval);
    if ((retval = nc_def_dim(ncid, "m2", grid.m2, &id_m2)))
        ERROR_NETCDF(retval);
    if ((retval = nc_def_dim(ncid, "dx1", grid.dx1, &id_dx1)))
        ERROR_NETCDF(retval);
    if ((retval = nc_def_dim(ncid, "dx2", grid.dx2, &id_dx2)))
        ERROR_NETCDF(retval);
    if ((retval = nc_def_dim(ncid, "two", 2, &id_np)))
        ERROR_NETCDF(retval);

    // Variables
    int varid_X, varid_V, varid_S, varid_names, varid_n1, varid_n2, varid_lim1, varid_lim2, varid_t, varid_dt;
    // X
    int dimids_X[2] = {id_r, id_dx1};
    if ((retval = nc_def_var(ncid, "X", NC_DOUBLE, 2, dimids_X, &varid_X)))
        ERROR_NETCDF(retval);
    // V
    int dimids_V[2] = {id_r, id_dx2};
    if ((retval = nc_def_var(ncid, "V", NC_DOUBLE, 2, dimids_V, &varid_V)))
        ERROR_NETCDF(retval);
    // S
    int dimids_S[2] = {id_r, id_r};
    if ((retval = nc_def_var(ncid, "S", NC_DOUBLE, 2, dimids_S, &varid_S)))
        ERROR_NETCDF(retval);
    // names
    if ((retval = nc_def_var(ncid, "names", NC_STRING, 1, &id_d, &varid_names)))
        ERROR_NETCDF(retval);
    // n1
    if ((retval = nc_def_var(ncid, "n1", NC_INT64, 1, &id_m1, &varid_n1)))
        ERROR_NETCDF(retval);
    // n2
    if ((retval = nc_def_var(ncid, "n2", NC_INT64, 1, &id_m2, &varid_n2)))
        ERROR_NETCDF(retval);
    // lim1
    int dimids_lim1[2] = {id_np, id_m1};
    if ((retval = nc_def_var(ncid, "lim1", NC_DOUBLE, 2, dimids_lim1, &varid_lim1)))
        ERROR_NETCDF(retval);
    // lim2
    int dimids_lim2[2] = {id_np, id_m2};
    if ((retval = nc_def_var(ncid, "lim2", NC_DOUBLE, 2, dimids_lim2, &varid_lim2)))
        ERROR_NETCDF(retval);
    // t
    if ((retval = nc_def_var(ncid, "t", NC_DOUBLE, 0, 0, &varid_t)))
        ERROR_NETCDF(retval);
    // dt
    if ((retval = nc_def_var(ncid, "dt", NC_DOUBLE, 0, 0, &varid_dt)))
        ERROR_NETCDF(retval);

    if ((retval = nc_enddef(ncid)))
        ERROR_NETCDF(retval);

    if ((retval = nc_put_var_double(ncid, varid_X, lr_sol.X.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_put_var_double(ncid, varid_V, lr_sol.V.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_put_var_double(ncid, varid_S, lr_sol.S.data())))
        ERROR_NETCDF(retval);
    size_t count = 1;
    for (size_t i = 0; i < names.size(); i++)
    {
        const char *name = names[i].c_str();
        if ((retval = nc_put_vara_string(ncid, varid_names, &i, &count, &name)))
            ERROR_NETCDF(retval);
    }
    if ((retval = nc_put_var_longlong(ncid, varid_n1, n1.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_put_var_longlong(ncid, varid_n2, n2.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_put_var_double(ncid, varid_lim1, grid.lim1.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_put_var_double(ncid, varid_lim2, grid.lim2.data())))
        ERROR_NETCDF(retval);
    if ((retval = nc_put_var_double(ncid, varid_t, t)))
        ERROR_NETCDF(retval);
    if ((retval = nc_put_var_double(ncid, varid_dt, dt)))
        ERROR_NETCDF(retval);
}

// void WriteNC(string fn, const lr2<double> &lr_sol, vector<string> names, grid_info grid, double *t, double *dt)
// {
//     multi_array<long long, 1> n1(grid.n1.shape()), n2(grid.n2.shape());

//     for (Index i = 0; i < grid.m1; i++)
//         n1(i) = (long long) grid.n1(i);

//     for (Index i = 0; i < grid.m2; i++)
//         n2(i) = (long long) grid.n2(i);

//     try
//     {
//         netCDF::NcFile file(fn, netCDF::NcFile::replace);
//         auto r_dim = file.addDim("r", grid.r);
//         auto d_dim = file.addDim("d", grid.d);
//         auto m1_dim = file.addDim("m1", grid.m1);
//         auto m2_dim = file.addDim("m2", grid.m2);
//         auto dx1_dim = file.addDim("dx1", grid.dx1);
//         auto dx2_dim = file.addDim("dx2", grid.dx2);

//         // NOTE: Ensign uses column-major order, but this interface uses row-major order,
//         // therefore the data has to be transposed!
//         auto x_data = file.addVar("X", netCDF::ncDouble, {r_dim, dx1_dim});
//         auto s_data = file.addVar("S", netCDF::ncDouble, {r_dim, r_dim});
//         auto v_data = file.addVar("V", netCDF::ncDouble, {r_dim, dx2_dim});
//         auto names_data = file.addVar("names", netCDF::ncString, {d_dim});
//         auto n1_data = file.addVar("n1", netCDF::ncInt64, {m1_dim});
//         auto n2_data = file.addVar("n2", netCDF::ncInt64, {m2_dim});
//         auto lim1_data = file.addVar("lim1", netCDF::ncDouble, {2, m1_dim});
//         auto lim2_data = file.addVar("lim2", netCDF::ncDouble, {2, m2_dim});
//         auto t_data = file.addVar("t", netCDF::ncDouble);
//         auto dt_data = file.addVar("dt", netCDF::ncDouble);

//         x_data.putVar(lr_sol.X.data());
//         s_data.putVar(lr_sol.S.data());
//         v_data.putVar(lr_sol.V.data());
//         names_data.putVar(names.data());
//         n1_data.putVar(n1.data());
//         n2_data.putVar(n2.data());
//         lim1_data.putVar(grid.lim1.data());
//         lim2_data.putVar(grid.lim2.data());
//         t_data.putVar(t);
//         dt_data.putVar(dt);
//     }
//     catch (netCDF::exceptions::NcException &e)
//     {
//         std::cout << e.what() << std::endl;
//     }
// }

// void read_nc(std::string fn, MatrixXd &U, MatrixXd &S, MatrixXd &V, vector<std::string> &names)
// {
//     int retval, ncid;

//     if ((retval = nc_open(fn.c_str(), NC_NOWRITE, &ncid)))
//         ERROR_NETCDF(retval);

//     // read dimensions
//     int id_d, id_r, id_n_part1, id_n_part2;
//     if ((retval = nc_inq_dimid(ncid, "d", &id_d)))
//         ERROR_NETCDF(retval);
//     if ((retval = nc_inq_dimid(ncid, "r", &id_r)))
//         ERROR_NETCDF(retval);
//     if ((retval = nc_inq_dimid(ncid, "n1", &id_n_part1)))
//         ERROR_NETCDF(retval);
//     if ((retval = nc_inq_dimid(ncid, "n2", &id_n_part2)))
//         ERROR_NETCDF(retval);

//     size_t d, r, n1, n2;
//     char tmp[NC_MAX_NAME + 1];
//     if ((retval = nc_inq_dim(ncid, id_d, tmp, &d)))
//         ERROR_NETCDF(retval);
//     if ((retval = nc_inq_dim(ncid, id_r, tmp, &r)))
//         ERROR_NETCDF(retval);
//     if ((retval = nc_inq_dim(ncid, id_n_part1, tmp, &n1)))
//         ERROR_NETCDF(retval);
//     if ((retval = nc_inq_dim(ncid, id_n_part2, tmp, &n2)))
//         ERROR_NETCDF(retval);

//     U.resize(n1, r);
//     V.resize(n2, r);
//     S.resize(r, r);

//     // read variables
//     int id_U, id_V, id_S, id_names;
//     if ((retval = nc_inq_varid(ncid, "U", &id_U)))
//         ERROR_NETCDF(retval);
//     if ((retval = nc_inq_varid(ncid, "V", &id_V)))
//         ERROR_NETCDF(retval);
//     if ((retval = nc_inq_varid(ncid, "S", &id_S)))
//         ERROR_NETCDF(retval);
//     if ((retval = nc_inq_varid(ncid, "names", &id_names)))
//         ERROR_NETCDF(retval);

//     if ((retval = nc_get_var_double(ncid, id_U, U.data())))
//         ERROR_NETCDF(retval);
//     if ((retval = nc_get_var_double(ncid, id_V, V.data())))
//         ERROR_NETCDF(retval);
//     if ((retval = nc_get_var_double(ncid, id_S, S.data())))
//         ERROR_NETCDF(retval);

//     // read names
//     vector<char *> ptr_names;
//     for (ind i = 0; i < (ind)d; i++)
//         ptr_names.push_back((char *)malloc(sizeof(char) * (NC_MAX_NAME + 1)));
//     nc_get_var_string(ncid, id_names, &ptr_names[0]);

//     for (ind i = 0; i < (ind)d; i++)
//     {
//         names.push_back(ptr_names[i]);
//         free(ptr_names[i]);
//     }

//     if ((retval = nc_close(ncid)))
//         ERROR_NETCDF(retval);
// }