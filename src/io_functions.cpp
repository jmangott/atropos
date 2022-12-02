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