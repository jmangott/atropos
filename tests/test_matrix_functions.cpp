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

using std::cout;
using std::endl;

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

    Matrix::Matricize(ten, mat2, 2);
    Matrix::Matricize(ten, mat1, 1);
    Matrix::Matricize(ten, mat0, 0);

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

    Matrix::Tensorize(mat2, ten2, 2);
    Matrix::Tensorize(mat1, ten1, 1);
    Matrix::Tensorize(mat0, ten0, 0);

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

    for(auto &el : vec1) cout << el << " ";
    cout << endl;

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
    ip = inner_product_from_const_weight(1.0, dx);
    blas_ops blas;

    multi_array<double, 2> Q(mat), R({r, r});
    R = Matrix::Orthogonalize(Q, n_basisfunctions, ip, blas);

    multi_array<double, 2> Q2({r, r}), id_r({r, r});
    set_identity(id_r);
    blas.matmul_transa(Q, Q, Q2);
    blas.matmul(Q, R, mat_ref);

    REQUIRE(bool(Q2 == id_r));
    REQUIRE(bool(mat == mat_ref));
}