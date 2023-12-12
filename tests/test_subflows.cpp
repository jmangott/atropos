#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <algorithm>
#include <cmath>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "matrix.hpp"
#include "subflows.hpp"
#include "tree_class.hpp"

TEST_CASE("subflows", "[subflows]")
{

    Index r = 3, r0 = 2;
    Index n_basisfunctions = 1, n_basisfunctions0 = 1;
    Index n_reactions = 6;

    Index val_n00 = 2, val_n01 = 2, val_n1 = 3;
    double val_liml00 = 0.0, val_liml01 = 0.0, val_liml1 = 0.0;
    Index val_binsize00 = 1, val_binsize01 = 1, val_binsize1 = 1;

    vector<Index> n{val_n00, val_n01, val_n1};
    vector<Index> n0{val_n00, val_n01};
    vector<Index> n1{val_n1};
    vector<Index> n00{val_n00};
    vector<Index> n01{val_n01};

    vector<double> liml{val_liml00, val_liml01, val_liml1};
    vector<double> liml0{val_liml00, val_liml01};
    vector<double> liml1{val_liml1};
    vector<double> liml00{val_liml00};
    vector<double> liml01{val_liml01};

    vector<Index> binsize{val_binsize00, val_binsize01, val_binsize1};
    vector<Index> binsize0{val_binsize00, val_binsize01};
    vector<Index> binsize1{val_binsize1};
    vector<Index> binsize00{val_binsize00};
    vector<Index> binsize01{val_binsize01};

    multi_array<bool, 2> dep({n_reactions, (Index)n.size()});
    multi_array<bool, 2> dep0({n_reactions, (Index)n0.size()});
    multi_array<bool, 2> dep1({n_reactions, (Index)n1.size()});
    multi_array<bool, 2> dep00({n_reactions, (Index)n00.size()});
    multi_array<bool, 2> dep01({n_reactions, (Index)n01.size()});

    std::fill(std::begin(dep), std::end(dep), false);
    std::fill(std::begin(dep0), std::end(dep0), false);
    std::fill(std::begin(dep1), std::end(dep1), false);
    std::fill(std::begin(dep00), std::end(dep00), false);
    std::fill(std::begin(dep01), std::end(dep01), false);

    dep(0, 0) = true;
    dep(1, 1) = true;
    dep(2, 2) = true;
    dep(3, 0) = true;
    dep(4, 1) = true;
    dep(5, 2) = true;

    dep0(0, 0) = true;
    dep0(1, 1) = true;
    dep0(3, 0) = true;
    dep0(4, 1) = true;

    dep1(2, 0) = true;
    dep1(5, 0) = true;

    dep00(0, 0) = true;
    dep00(3, 0) = true;

    dep01(1, 0) = true;
    dep01(4, 0) = true;

    multi_array<Index, 2> nu({n_reactions, (Index)n.size()});
    multi_array<Index, 2> nu0({n_reactions, (Index)n0.size()});
    multi_array<Index, 2> nu1({n_reactions, (Index)n1.size()});
    multi_array<Index, 2> nu00({n_reactions, (Index)n00.size()});
    multi_array<Index, 2> nu01({n_reactions, (Index)n01.size()});

    std::fill(std::begin(nu), std::end(nu), 0);
    std::fill(std::begin(nu0), std::end(nu0), 0);
    std::fill(std::begin(nu1), std::end(nu1), 0);
    std::fill(std::begin(nu00), std::end(nu00), 0);
    std::fill(std::begin(nu01), std::end(nu01), 0);

    nu(0, 0) = -1;
    nu(1, 1) = -1;
    nu(2, 2) = -1;
    nu(3, 0) = 1;
    nu(4, 1) = 1;
    nu(5, 2) = 1;

    nu0(0, 0) = -1;
    nu0(1, 1) = -1;
    nu0(3, 0) = 1;
    nu0(4, 1) = 1;

    nu1(2, 0) = -1;
    nu1(5, 0) = 1;

    nu00(0, 0) = -1;
    nu00(3, 0) = 1;

    nu01(1, 0) = -1;
    nu01(4, 0) = 1;

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

    std::vector<std::vector<double>> propensity(grid.n_reactions);
    std::vector<std::vector<double>> propensity0(grid.n_reactions);
    std::vector<std::vector<double>> propensity1(grid.n_reactions);
    std::vector<std::vector<double>> propensity00(grid.n_reactions);
    std::vector<std::vector<double>> propensity01(grid.n_reactions);

    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        propensity[mu].resize(grid.dx_dep[mu]);
        propensity0[mu].resize(grid.dx_dep[mu]);
        propensity1[mu].resize(grid.dx_dep[mu]);
        propensity00[mu].resize(grid.dx_dep[mu]);
        propensity01[mu].resize(grid.dx_dep[mu]);
    }

    propensity[0] = {0.0, 1.0};
    propensity[1] = {0.0, 1.0};
    propensity[2] = {0.0, 1.0, 2.0};
    propensity[3] = {1.0, 0.5};
    propensity[4] = {1.0, 0.5};
    propensity[5] = {1.0, 0.5, 1.0/3};

    propensity0[0] = {0.0, 1.0};
    propensity0[1] = {0.0, 1.0};
    propensity0[2] = {1.0};
    propensity0[3] = {1.0, 0.5};
    propensity0[4] = {1.0, 0.5};
    propensity0[5] = {1.0};

    propensity1[0] = {1.0};
    propensity1[1] = {1.0};
    propensity1[2] = {0.0, 1.0, 2.0};
    propensity1[3] = {1.0};
    propensity1[4] = {1.0};
    propensity1[5] = {1.0, 0.5, 1.0 / 3};

    propensity00[0] = {0.0, 1.0};
    propensity00[1] = {1.0};
    propensity00[2] = {1.0};
    propensity00[3] = {1.0, 0.5};
    propensity00[4] = {1.0};
    propensity00[5] = {1.0};

    propensity01[0] = {1.0};
    propensity01[1] = {0.0, 1.0};
    propensity01[2] = {1.0};
    propensity01[3] = {1.0};
    propensity01[4] = {1.0, 0.5};
    propensity01[5] = {1.0};

    multi_array<double, 3> Q({r, r, 1}), Q0({r0, r0, r});
    multi_array<double, 2> X00({val_n00, r0}), X01({val_n01, r0}), X1({val_n1, r});
    multi_array<double, 2> Q0_mat({r0 * r0, r});

    // Initialize Q and Q0
    std::fill(std::begin(Q), std::end(Q), 0.0);
    std::fill(std::begin(Q0), std::end(Q0), 0.0);
    Q(0, 0, 0) = 1.0;
    Q0(0, 0, 0) = 1.0;

    // Initialize and normalize X00, X01, X1
    set_zero(X00);
    set_zero(X01);
    set_zero(X1);

    X00(0, 0) = 1.0;
    X00(1, 0) = 1.0;
    X00 *= std::exp(-0.25);

    X01(0, 0) = 1.0;
    X01(1, 0) = 1.0;
    X01 *= std::exp(-0.25);

    X1(0, 0) = 1.0;
    X1(1, 0) = 1.0;
    X1(2, 0) = std::exp(-2.0);
    X1 *= std::exp(-0.25);

    // Construct cme_lr_tree
    cme_internal_node *root = new cme_internal_node("", nullptr, grid, 1, {r, r}, 1);
    cme_internal_node *node0 = new cme_internal_node("0", root, grid0, r, {r0, r0}, n_basisfunctions);
    cme_external_node *node1 = new cme_external_node("1", root, grid1, r, n_basisfunctions);
    cme_external_node *node00 = new cme_external_node("00", node0, grid00, r0, n_basisfunctions0);
    cme_external_node *node01 = new cme_external_node("01", node0, grid01, r0, n_basisfunctions0);

    root->Q = Q;
    node0->Q = Q0;
    node1->X = X1;
    node00->X = X00;
    node01->X = X01;

    root->coefficients.propensity = propensity;
    node0->coefficients.propensity = propensity0;
    node1->coefficients.propensity = propensity1;
    node00->coefficients.propensity = propensity00;
    node01->coefficients.propensity = propensity01;

    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        root->coefficients.A(mu, 0, 0) = 1.0;
        root->coefficients.B(mu, 0, 0) = 1.0;
    }
    root->coefficients.E(0, 0, 0, 0) = 1.0;
    root->coefficients.F(0, 0, 0, 0) = 1.0;

    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        node00->external_coefficients.C[mu].resize({node00->grid.dx_dep[mu], node00->RankIn(), node00->RankIn()});
        node00->external_coefficients.D[mu].resize({node00->grid.dx_dep[mu], node00->RankIn(), node00->RankIn()});

        node01->external_coefficients.C[mu].resize({node01->grid.dx_dep[mu], node01->RankIn(), node01->RankIn()});
        node01->external_coefficients.D[mu].resize({node01->grid.dx_dep[mu], node01->RankIn(), node01->RankIn()});

        node1->external_coefficients.C[mu].resize({node1->grid.dx_dep[mu], node1->RankIn(), node1->RankIn()});
        node1->external_coefficients.D[mu].resize({node1->grid.dx_dep[mu], node1->RankIn(), node1->RankIn()});
    }

    root->child[0] = node0;
    root->child[1] = node1;
    root->child[0]->child[0] = node00;
    root->child[0]->child[1] = node01;

    cme_lr_tree tree(root);

    // Calculate probability distribution
    multi_array<double, 3> p({val_n00, val_n01, val_n1}), p_ortho({val_n00, val_n01, val_n1});
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
                                p_ortho(x00, x01, x1) += tree.root->Q(i, j, 0) * ((cme_internal_node *)tree.root->child[0])->Q(i0, j0, i) * ((cme_external_node *)tree.root->child[0]->child[0])->X(x00, i0) * ((cme_external_node *)tree.root->child[0]->child[1])->X(x01, j0) * ((cme_external_node *)tree.root->child[1])->X(x1, j);
                            }
                        }
                    }
                }
            }
        }
    }

    REQUIRE(bool(p == p_ortho));

    // Orthogonalize Q and X manually
    std::fill(std::begin(Q), std::end(Q), 0.0);
    std::fill(std::begin(Q0), std::end(Q0), 0.0);
    set_zero(Q0_mat);

    // TODO: correct error in notes
    Q(0, 0, 0) = (2.0 * std::exp(-0.5) * std::exp(-0.25) * sqrt(2.0 + std::exp(-4.0)));
    Q0_mat(0, 0) = 1.0;
    Q0_mat(1, 1) = 1.0;
    Q0_mat(2, 2) = 1.0;
    Matrix::Tensorize(Q0_mat, Q0, 2);

    set_zero(X00);
    set_zero(X01);
    set_zero(X1);
    
    X00(0, 0) = 1.0;
    X00(1, 0) = 1.0;
    X00(0, 1) = 1.0;
    X00(1, 1) =-1.0;
    X00 /= sqrt(2.0);

    X01(0, 0) = 1.0;
    X01(1, 0) = 1.0;
    X01(0, 1) = 1.0;
    X01(1, 1) =-1.0;
    X01 /= sqrt(2.0);

    X1(0, 0) = 1.0;
    X1(1, 0) = 1.0;
    X1(2, 0) = std::exp(-2.0);
    X1(0, 1) = sqrt(1.0 + 0.5 * std::exp(-4.0));
    X1(1, 1) =-sqrt(1.0 + 0.5 * std::exp(-4.0));
    X1(2, 1) = 0.0;
    X1(0, 2) = std::exp(-2.0) / sqrt(2.0);
    X1(1, 2) = std::exp(-2.0) / sqrt(2.0);
    X1(2, 2) = -sqrt(2.0);
    X1 /= sqrt(2.0 + std::exp(-4.0));

    root->Q = Q;
    node0->Q = Q0;
    node1->X = X1;
    node00->X = X00;
    node01->X = X01;

    // Check if this yields the same probability distribution
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

    REQUIRE(bool(p == p_ortho));

    blas_ops blas;
    SubflowPhi<0>(root, blas, 1.0);

    // for (Index i = 0; i < root->RankIn(); ++i)
    // {
    //     for (Index i0 = 0; i0 < root->RankOut()[0]; ++i0)
    //     {
    //         for (Index i1 = 0; i1 < root->RankOut()[1]; ++i1)
    //         {
    //             cout << root->Q(i0, i1, i) << " ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }

    // cout << endl;

    // for (Index i = 0; i < root->child[0]->RankIn(); ++i)
    // {
    //     for (Index i0 = 0; i0 < ((cme_internal_node*) root->child[0])->RankOut()[0]; ++i0)
    //     {
    //         for (Index i1 = 0; i1 < ((cme_internal_node*) root->child[0])->RankOut()[1]; ++i1)
    //         {
    //             cout << ((cme_internal_node*) root->child[0])->Q(i0, i1, i) << " ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }
    // cout << endl;



    // blas_ops blas;
    // TTNIntegrator(root, blas, 1.0);
    // SubflowPhi<0>(root, blas, 1.0);

    // for (Index mu = 0; mu < node1->grid.n_reactions; ++mu)
    // {
    //     cout << "mu: " << mu << endl;
    //     for (Index x = 0; x < node1->grid.dx_dep[mu]; ++x)
    //     {
    //         cout << "x: " << x << endl;
    //         for (Index i = 0; i < node1->RankIn(); ++i)
    //         {
    //             for (Index j = 0; j < node1->RankIn(); ++j)
    //             {
    //                 cout << node1->external_coefficients.C[mu](x, i, j) << " ";
    //             }
    //             cout << endl;
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }
}