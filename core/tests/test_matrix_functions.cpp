#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <algorithm>
#include <cmath>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "matrix.hpp"
#include "tree_class.hpp"

class generator
{
    private:
        Index i = 0;

    public:
        inline Index operator()()
        {
            return ++i;
        };
};

TEST_CASE("Matricize_Tensorize", "[Matricize_Tensorize]")
{
    Index d0 = 5, d1 = 7, d2 = 9;
    multi_array<Index, 3> ten({d0, d1, d2}), ten2({d0, d1, d2}), ten1({d0, d1, d2}), ten0({d0, d1, d2});
    multi_array<Index, 2> mat2({d0 * d1, d2}), mat2_ref({d0 * d1, d2});
    multi_array<Index, 2> mat1({d2 * d0, d1}), mat1_ref({d2 * d0, d1});
    multi_array<Index, 2> mat0({d1 * d2, d0}), mat0_ref({d1 * d2, d0});

    std::generate(std::begin(ten), std::end(ten), generator{});

    Matrix::Matricize<2>(ten, mat2);
    Matrix::Matricize<1>(ten, mat1);
    Matrix::Matricize<0>(ten, mat0);

    generator generator2{};
    for (Index k = 0; k < d2; ++k)
    {
        for (Index j = 0; j < d1; ++j)
        {
            for (Index i = 0; i < d0; ++i)
            {
                mat2_ref(i + j * d0, k) = generator2();
            }
        }
    }
    
    generator generator1{};
    for (Index k = 0; k < d2; ++k)
    {
        for (Index j = 0; j < d1; ++j)
        {
            for (Index i = 0; i < d0; ++i)
            {
                mat1_ref(i + k * d0, j) = generator1();
            }
        }
    }

    generator generator0{};
    for (Index k = 0; k < d2; ++k)
    {
        for (Index j = 0; j < d1; ++j)
        {
            for (Index i = 0; i < d0; ++i)
            {
                mat0_ref(j + k * d1, i) = generator0();
            }
        }
    }

    REQUIRE(bool(mat2 == mat2_ref));
    REQUIRE(bool(mat1 == mat1_ref));
    REQUIRE(bool(mat0 == mat0_ref));

    Matrix::Tensorize<2>(mat2, ten2);
    Matrix::Tensorize<1>(mat1, ten1);
    Matrix::Tensorize<0>(mat0, ten0);

    REQUIRE(bool(ten2 == ten));
    REQUIRE(bool(ten1 == ten));
    REQUIRE(bool(ten0 == ten));
}

TEST_CASE("RemoveElement", "[RemoveElement]")
{
    std::array<Index, 10> start_vec;
    std::array<Index, 9> vec1, vec2, vec3, vec1_ref, vec2_ref, vec3_ref;

    std::generate(std::begin(start_vec), std::end(start_vec), generator{});

    Matrix::RemoveElement(std::begin(start_vec), std::end(start_vec), std::begin(vec1), 0);
    Matrix::RemoveElement(std::begin(start_vec), std::end(start_vec), std::begin(vec2), 5);
    Matrix::RemoveElement(std::begin(start_vec), std::end(start_vec), std::begin(vec3), 9);

    vec1_ref = {2, 3, 4, 5, 6, 7, 8, 9, 10};
    vec2_ref = {1, 2, 3, 4, 5, 7, 8, 9, 10};
    vec3_ref = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    REQUIRE(bool(vec1 == vec1_ref));
    REQUIRE(bool(vec2 == vec2_ref));
    REQUIRE(bool(vec3 == vec3_ref));

}

TEST_CASE("Orthogonalize", "[Orthogonalize]")
{
    Index r = 5, dx = 6;
    Index n_basisfunctions = 1;
    multi_array<double, 2> mat({dx, r}), mat_ref({dx, r});
    set_zero(mat);
    set_zero(mat_ref);
    mat(0, 0) = 1.0;
    std::function<double (double*, double*)> ip;
    blas_ops blas;

    multi_array<double, 2> Q(mat), R({r, r});
    R = Matrix::Orthogonalize(Q, n_basisfunctions, 1.0, blas);

    multi_array<double, 2> Q2({r, r}), id_r({r, r});
    set_identity(id_r);
    blas.matmul_transa(Q, Q, Q2);
    blas.matmul(Q, R, mat_ref);

    REQUIRE(bool(Q2 == id_r));
    REQUIRE(bool(mat == mat_ref));
}

TEST_CASE("ShiftRows", "[ShiftRows]")
{
    Index d = 2;
    Index n_reactions = 4;
    std::vector<Index> n(d);
    std::vector<Index> binsize(d);
    std::vector<double> liml(d);
    multi_array<bool, 2> dep({n_reactions, d});
    multi_array<Index, 2> nu({n_reactions, d});
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
    multi_array<double, 2> input_array({n_rows, 1}), output_array({n_rows, 1});
    multi_array<double, 2> comparison_array({n_rows, 1});

    set_zero(input_array);
    for (Index i = 0; i < n_rows; ++i)
        input_array(i, 0) = (double)i + 1;

    set_zero(comparison_array);
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

    set_zero(input_array);
    for (Index i = 0; i < n_rows; ++i)
    {
        input_array(i, 0) = (double) i + 1;
        input_array(i, 1) = (double) i + 2;
    }

    set_zero(comparison_array);
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