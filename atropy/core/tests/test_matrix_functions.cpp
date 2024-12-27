#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include <algorithm>
#include <cmath>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <generic/tensor.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "matrix.hpp"
#include "tree_class.hpp"

class generator {
  private:
    Index i = 0;

  public:
    inline Index operator()() { return ++i; };
};

TEST_CASE("ShiftRows", "[ShiftRows]")
{
    Index d = 2;
    Index n_reactions = 4;
    std::vector<Index> n(d);
    std::vector<Index> binsize(d);
    std::vector<double> liml(d);
    Ensign::multi_array<bool, 2> dep({n_reactions, d});
    Ensign::multi_array<Index, 2> nu({n_reactions, d});
    std::vector<int> species(d);

    n = {4, 3};
    binsize = {1, 1};
    liml = {0.0, 0.0};
    species = {0, 1};

    std::fill(std::begin(dep), std::end(dep), false);
    dep(0, 0) = true;
    dep(0, 1) = true;
    dep(1, 1) = true;
    dep(2, 1) = true;
    dep(3, 0) = true;

    std::fill(std::begin(nu), std::end(nu), 0);
    nu(0, 0) = -1;
    nu(0, 1) = -1;
    nu(1, 1) = -1;
    nu(2, 0) = 1;
    nu(3, 1) = 1;

    grid_parms grid(n, binsize, liml, dep, nu, species);
    grid.Initialize();

    // TEST 1
    Index n_rows = grid.dx;
    Ensign::multi_array<double, 2> input_array({n_rows, 1}), output_array({n_rows, 1});
    Ensign::multi_array<double, 2> comparison_array({n_rows, 1});

    Ensign::Matrix::set_zero(input_array);
    for (Index i = 0; i < n_rows; ++i)
        input_array(i, 0) = (double)i + 1;

    Ensign::Matrix::set_zero(comparison_array);
    comparison_array(0, 0) = 6.0;
    comparison_array(1, 0) = 7.0;
    comparison_array(2, 0) = 8.0;
    comparison_array(3, 0) = 0.0; // !
    comparison_array(4, 0) = 10.0;
    comparison_array(5, 0) = 11.0;
    comparison_array(6, 0) = 12.0;

    Matrix::ShiftRows<1>(output_array, input_array, grid, 0);
    REQUIRE(bool(output_array == comparison_array));

    // TEST 2
    input_array.resize({n_rows, 2}), output_array.resize({n_rows, 2});
    comparison_array.resize({n_rows, 2});

    Ensign::Matrix::set_zero(input_array);
    for (Index i = 0; i < n_rows; ++i) {
        input_array(i, 0) = (double)i + 1;
        input_array(i, 1) = (double)i + 2;
    }

    Ensign::Matrix::set_zero(comparison_array);
    comparison_array(5, 0) = 1.0;
    comparison_array(6, 0) = 2.0;
    comparison_array(7, 0) = 3.0;
    comparison_array(8, 0) = 0.0; // !
    comparison_array(9, 0) = 5.0;
    comparison_array(10, 0) = 6.0;
    comparison_array(11, 0) = 7.0;

    comparison_array(5, 1) = 2.0;
    comparison_array(6, 1) = 3.0;
    comparison_array(7, 1) = 4.0;
    comparison_array(8, 1) = 0.0; // !
    comparison_array(9, 1) = 6.0;
    comparison_array(10, 1) = 7.0;
    comparison_array(11, 1) = 8.0;

    Matrix::ShiftRows<-1>(output_array, input_array, grid, 0);

    REQUIRE(bool(output_array == comparison_array));
}
