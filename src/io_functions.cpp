#include "io_functions.hpp"

using std::ifstream;
using std::stringstream;

// TODO: Use .netcdf instead of .csv
void ReadInMultiArray(multi_array<double, 2> &output_array, string filename)
{
    Index ii = 0;
    Index jj = 0;
    string line;
    ifstream input_file(filename);
    stringstream ss_line;
    string element;

    while (getline(input_file, line))
    {
        ss_line.str(line);
        while (getline(ss_line, element, ','))
        {
            output_array(ii, jj) = stod(element);
            jj++;
        }
        ii++;
        jj = 0;
        ss_line.clear();
    }
    input_file.close();
}