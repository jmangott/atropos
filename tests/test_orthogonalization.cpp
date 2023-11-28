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

class initial_distribution
{
    private:
        double liml;
        double binsize;
        Index i = 0;
    public:
        initial_distribution(const double _liml, const Index _binsize) : liml(_liml), binsize(_binsize) {};
        inline double operator()()
        {
            double output = std::exp(-std::pow(i + liml, 2));
            i = i + binsize;
            return output;
        };
};

TEST_CASE("orthogonalization", "[orthogonalization]")
{
    Index r = 5, r0 = 4;
    Index n_basisfunctions = 1, n_basisfunctions0 = 1;
    Index n_reactions = 1;

    Index val_n00 = 11, val_n01 = 11, val_n1 = 21;
    double val_liml00 = 0.0, val_liml01 = 0.0, val_liml1 = 0.0;
    Index val_binsize00 = 1, val_binsize01 = 1, val_binsize1 = 1;

    vector<Index> n {val_n00, val_n01, val_n1};
    vector<Index> n0 {val_n00, val_n01};
    vector<Index> n1 {val_n1};
    vector<Index> n00 {val_n00};
    vector<Index> n01 {val_n01};

    vector<double> liml {val_liml00, val_liml01, val_liml1};
    vector<double> liml0 {val_liml00, val_liml01};
    vector<double> liml1 {val_liml1};
    vector<double> liml00 {val_liml00};
    vector<double> liml01 {val_liml01};

    vector<Index> binsize {val_binsize00, val_binsize01, val_binsize1};
    vector<Index> binsize0 {val_binsize00, val_binsize01};
    vector<Index> binsize1 {val_binsize1};
    vector<Index> binsize00 {val_binsize00};
    vector<Index> binsize01 {val_binsize01};

    multi_array<bool, 2> dep({n_reactions, n.size()});
    multi_array<bool, 2> dep0({n_reactions, n0.size()});
    multi_array<bool, 2> dep1({n_reactions, n1.size()});
    multi_array<bool, 2> dep00({n_reactions, n00.size()});
    multi_array<bool, 2> dep01({n_reactions, n01.size()});

    multi_array<Index, 2> nu({n_reactions, n.size()});
    multi_array<Index, 2> nu0({n_reactions, n0.size()});
    multi_array<Index, 2> nu1({n_reactions, n1.size()});
    multi_array<Index, 2> nu00({n_reactions, n00.size()});
    multi_array<Index, 2> nu01({n_reactions, n01.size()});

    grid_parms grid(n, binsize, liml, dep, nu);
    grid_parms grid0(n0, binsize0, liml0, dep0, nu0);
    grid_parms grid1(n1, binsize1, liml1, dep1, nu1);
    grid_parms grid00(n00, binsize00, liml00, dep00, nu00);
    grid_parms grid01(n01, binsize01, liml01, dep01, nu01);
    grid.Initialize();
    grid0.Initialize();
    grid1.Initialize();
    grid00.Initialize();
    grid01.Initialize();

    multi_array<double, 3> p({val_n00, val_n01, val_n1}), p_ortho({val_n00, val_n01, val_n1});
    multi_array<double, 3> Q({r, r, 1}), Q0({r0, r0, r});
    multi_array<double, 2> X00({val_n00, r0}), X01({val_n01, r0}), X1({val_n1, r});
    multi_array<double, 2> Q0_mat({r0 * r0, r});

    // Initialize Q and Q0
    std::fill(std::begin(Q), std::end(Q), 0.0);
    std::fill(std::begin(Q0), std::end(Q0), 0.0);
    Q(0, 0, 0) = 1.0;
    Q0(0, 0, 1) = 1.0;

    // Initialize and normalize X00, X01, X1
    set_zero(X00);
    set_zero(X01);
    set_zero(X1);
    std::function<double(double*, double*)> ip00 = inner_product_from_const_weight(grid00.h_mult, grid00.dx);
    std::function<double(double*, double*)> ip01 = inner_product_from_const_weight(grid01.h_mult, grid01.dx);
    std::function<double(double*, double*)> ip1 = inner_product_from_const_weight(grid1.h_mult, grid1.dx);
    std::function<double(double *, double *)> ip0 = inner_product_from_const_weight(1.0, r0 * r0);

    std::generate(std::begin(X00), std::begin(X00) + val_n00, initial_distribution(val_liml00, val_binsize00));
    double gamma00 = ip00(X00.extract({0}), X00.extract({0}));
    X00 /= gamma00;

    std::generate(std::begin(X01), std::begin(X01) + val_n01, initial_distribution(val_liml01, val_binsize01));
    double gamma01 = ip01(X01.extract({0}), X01.extract({0}));
    X01 /= gamma01;

    std::generate(std::begin(X1), std::begin(X1) + val_n1, initial_distribution(val_liml1, val_binsize1));
    double gamma1 = ip1(X1.extract({0}), X1.extract({0}));
    X1 /= gamma1;

    Matrix::Matricize(Q0, Q0_mat, 2);
    double gamma0 = ip0(Q0_mat.extract({0}), Q0_mat.extract({0}));
    Q0 /= gamma0;

    // Construct cme_lr_tree
    cme_internal_node* root = new cme_internal_node("", nullptr, grid, 1, {r, r}, 1);
    cme_internal_node* node0 = new cme_internal_node("0", root, grid0, r, {r0, r0}, n_basisfunctions);
    cme_external_node* node1 = new cme_external_node("1", root, grid1, r, n_basisfunctions);
    cme_external_node* node00 = new cme_external_node("00", node0, grid00, r0, n_basisfunctions0);
    cme_external_node* node01 = new cme_external_node("01", node0, grid01, r0, n_basisfunctions0);

    root->Q = Q;
    node0->Q = Q0;
    node1->X = X1;
    node00->X = X00;
    node01->X = X01;

    root->child[0] = node0;
    root->child[1] = node1;
    root->child[0]->child[0] = node00;
    root->child[0]->child[1] = node01;
    cme_lr_tree tree(root);

    // Calculate probability distribution
    std::fill(std::begin(p), std::end(p), 0.0);
    
    for (Index i = 0; i < r; ++i)
    {
        for (Index j = 0; j < r; ++j)
        {
            for (Index i0 = 0; i0 < r0; ++i0)
            {
                for (Index j0 = 0; j0 < r0; ++j0)
                {
                    for (Index x00 = 0; x00 < val_n00; ++x00)
                    {
                        for (Index x01 = 0; x01 < val_n01; ++x01)
                        {
                            for (Index x1 = 0; x1 < val_n1; ++x1)
                            {
                                p(x00, x01, x1) += Q(i, j, 0) * Q0(i0, j0, i) * X00(x00, i0) * X01(x01, j0) * X1(x1, j);
                            }
                        }
                    }
                }
            }
        }
    }

    double norm = std::accumulate(std::begin(p), std::end(p), 0.0);
    p /= norm;

    // Normalize Q
    tree.root->Q /= norm;

    // Check if the probability distribution remains the same under orthogonalization
    std::fill(std::begin(p_ortho), std::end(p_ortho), 0.0);
    tree.Orthogonalize();

    for (Index i = 0; i < r; ++i)
    {
        for (Index j = 0; j < r; ++j)
        {
            for (Index i0 = 0; i0 < r0; ++i0)
            {
                for (Index j0 = 0; j0 < r0; ++j0)
                {
                    for (Index x00 = 0; x00 < val_n00; ++x00)
                    {
                        for (Index x01 = 0; x01 < val_n01; ++x01)
                        {
                            for (Index x1 = 0; x1 < val_n1; ++x1)
                            {
                                p_ortho(x00, x01, x1) 
                                += tree.root->Q(i, j, 0) 
                                * ((cme_internal_node*) tree.root->child[0])->Q(i0, j0, i) 
                                * ((cme_external_node*) tree.root->child[0]->child[0])->X(x00, i0) 
                                * ((cme_external_node*) tree.root->child[0]->child[1])->X(x01, j0) 
                                * ((cme_external_node*) tree.root->child[1])->X(x1, j);
                            }
                        }
                    }
                }
            }
        }
    }

    REQUIRE(bool(p == p_ortho));

    // Check if the Xs and Qs are orthonormal
    blas_ops blas;

    multi_array<double, 2> X00_ortho({r0, r0}), X01_ortho({r0, r0}), X1_ortho({r, r}), Q_ortho({1, 1}), Q0_ortho({r, r});
    set_zero(X00_ortho);
    set_zero(X01_ortho);
    set_zero(X1_ortho);
    set_zero(Q0_ortho);

    blas.matmul_transa(
        ((cme_external_node*) tree.root->child[0]->child[0])->X,
        ((cme_external_node*) tree.root->child[0]->child[0])->X,
        X00_ortho
        );

    blas.matmul_transa(
        ((cme_external_node*) tree.root->child[0]->child[1])->X,
        ((cme_external_node*) tree.root->child[0]->child[1])->X,
        X01_ortho
        );

    blas.matmul_transa(
        ((cme_external_node*) tree.root->child[1])->X,
        ((cme_external_node*) tree.root->child[1])->X,
        X1_ortho
        );

    Matrix::Matricize(((cme_internal_node*) tree.root->child[0])->Q, Q0_mat, 2);
    blas.matmul_transa(Q0_mat, Q0_mat, Q0_ortho);

    multi_array<double, 2> id_r({r, r}), id_r0({r0, r0}), id_1({1, 1});
    set_identity(id_r);
    set_identity(id_r0);

    REQUIRE(bool(X00_ortho == id_r0));
    REQUIRE(bool(X01_ortho == id_r0));
    REQUIRE(bool(X1_ortho == id_r));
    REQUIRE(bool(Q0_ortho == id_r));
}