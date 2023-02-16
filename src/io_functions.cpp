#include "io_functions.hpp"

using std::ifstream;
using std::stringstream;

// TODO: Use .netcdf instead of .csv
void ReadInMultiArray(multi_array<double, 2> &output_array, string filename)
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


void WriteOutMultiArray(const multi_array<double, 2> &input_array, string filename)
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


// void save_lr(string fn, const lr2<double> &lr_sol, const grid_info &gi)
// {
//     nc_writer ncw(fn, {gi.N_xx[0], gi.N_xx[1], gi.N_xx[2], gi.N_vv[0], gi.N_vv[1], gi.N_vv[2], gi.r}, {"x", "y", "z", "v", "w", "u", "r"});
//     ncw.add_var("r", {"r"});
//     ncw.add_var("x", {"x"});
//     ncw.add_var("y", {"y"});
//     ncw.add_var("z", {"z"});
//     ncw.add_var("u", {"u"});
//     ncw.add_var("v", {"v"});
//     ncw.add_var("w", {"w"});
//     ncw.add_var("X", {"r", "z", "y", "x"});
//     ncw.add_var("S", {"r", "r"});
//     ncw.add_var("V", {"r", "u", "w", "v"});

//     ncw.start_write_mode();

//     vector<double> vec_r(gi.r);
//     for (Index i = 0; i < gi.r; i++)
//         vec_r[i] = i;

//     vector<double> vec_x(gi.N_xx[0]), vec_y(gi.N_xx[1]), vec_z(gi.N_xx[2]);
//     for (Index i = 0; i < gi.N_xx[0]; i++)
//         vec_x[i] = gi.x(0, i);
//     for (Index i = 0; i < gi.N_xx[1]; i++)
//         vec_y[i] = gi.x(1, i);
//     for (Index i = 0; i < gi.N_xx[2]; i++)
//         vec_z[i] = gi.x(2, i);

//     vector<double> vec_v(gi.N_vv[0]), vec_w(gi.N_vv[1]), vec_u(gi.N_vv[2]);
//     for (Index i = 0; i < gi.N_vv[0]; i++)
//         vec_v[i] = gi.v(0, i);
//     for (Index i = 0; i < gi.N_vv[1]; i++)
//         vec_w[i] = gi.v(1, i);
//     for (Index i = 0; i < gi.N_vv[2]; i++)
//         vec_u[i] = gi.v(2, i);

//     ncw.write("r", vec_r.data());
//     ncw.write("x", vec_x.data());
//     ncw.write("y", vec_y.data());
//     ncw.write("z", vec_z.data());
//     ncw.write("v", vec_v.data());
//     ncw.write("w", vec_w.data());
//     ncw.write("u", vec_u.data());

//     ncw.write("X", lr_sol.X.data());
//     ncw.write("S", lr_sol.S.data());
//     ncw.write("V", lr_sol.V.data());
// }