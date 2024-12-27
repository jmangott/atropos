#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <generic/tensor.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "matrix.hpp"
#include "tree_class.hpp"

TEST_CASE("tree_h2", "[tree_h2]")
{
    Ensign::Matrix::blas_ops blas;

    Index r = 3, r0 = 2;
    Index n_basisfunctions = 1, n_basisfunctions0 = 1;
    Index n_reactions = 6;

    Index val_n00 = 2, val_n01 = 2, val_n1 = 3;
    double val_liml00 = 0.0, val_liml01 = 0.0, val_liml1 = 0.0;
    Index val_binsize00 = 1, val_binsize01 = 1, val_binsize1 = 1;
    int val_species00 = 0, val_species01 = 1, val_species1 = 2;

    std::vector<Index> n{val_n00, val_n01, val_n1};
    std::vector<Index> n0{val_n00, val_n01};
    std::vector<Index> n1{val_n1};
    std::vector<Index> n00{val_n00};
    std::vector<Index> n01{val_n01};

    std::vector<double> liml{val_liml00, val_liml01, val_liml1};
    std::vector<double> liml0{val_liml00, val_liml01};
    std::vector<double> liml1{val_liml1};
    std::vector<double> liml00{val_liml00};
    std::vector<double> liml01{val_liml01};

    std::vector<Index> binsize{val_binsize00, val_binsize01, val_binsize1};
    std::vector<Index> binsize0{val_binsize00, val_binsize01};
    std::vector<Index> binsize1{val_binsize1};
    std::vector<Index> binsize00{val_binsize00};
    std::vector<Index> binsize01{val_binsize01};

    std::vector<int> species{val_species00, val_species01, val_species1};
    std::vector<int> species0{val_species00, val_species01};
    std::vector<int> species1{val_species1};
    std::vector<int> species00{val_species00};
    std::vector<int> species01{val_species01};

    Ensign::multi_array<bool, 2> dep({n_reactions, (Index)n.size()});
    Ensign::multi_array<bool, 2> dep0({n_reactions, (Index)n0.size()});
    Ensign::multi_array<bool, 2> dep1({n_reactions, (Index)n1.size()});
    Ensign::multi_array<bool, 2> dep00({n_reactions, (Index)n00.size()});
    Ensign::multi_array<bool, 2> dep01({n_reactions, (Index)n01.size()});

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

    Ensign::multi_array<Index, 2> nu({n_reactions, (Index)n.size()});
    Ensign::multi_array<Index, 2> nu0({n_reactions, (Index)n0.size()});
    Ensign::multi_array<Index, 2> nu1({n_reactions, (Index)n1.size()});
    Ensign::multi_array<Index, 2> nu00({n_reactions, (Index)n00.size()});
    Ensign::multi_array<Index, 2> nu01({n_reactions, (Index)n01.size()});

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

    grid_parms grid(n, binsize, liml, dep, nu, species);
    grid_parms grid0(n0, binsize0, liml0, dep0, nu0, species0);
    grid_parms grid1(n1, binsize1, liml1, dep1, nu1, species1);
    grid_parms grid00(n00, binsize00, liml00, dep00, nu00, species00);
    grid_parms grid01(n01, binsize01, liml01, dep01, nu01, species01);
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

    for (Index mu = 0; mu < grid.n_reactions; ++mu) {
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
    propensity[5] = {1.0, 0.5, 1.0 / 3};

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

    Ensign::multi_array<double, 3> Q({r, r, 1}), Q0({r0, r0, r});
    Ensign::multi_array<double, 2> X00({val_n00, r0}), X01({val_n01, r0}),
        X1({val_n1, r});
    Ensign::multi_array<double, 2> Q0_mat({r0 * r0, r});

    // Initialize Q and Q0
    std::fill(std::begin(Q), std::end(Q), 0.0);
    std::fill(std::begin(Q0), std::end(Q0), 0.0);
    Q(0, 0, 0) = 1.0;
    Q0(0, 0, 0) = 1.0;

    // Initialize and normalize X00, X01, X1
    Ensign::Matrix::set_zero(X00);
    Ensign::Matrix::set_zero(X01);
    Ensign::Matrix::set_zero(X1);

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
    cme_internal_node* root = new cme_internal_node("", nullptr, grid, 1, {r, r}, 1);
    cme_internal_node* node0 =
        new cme_internal_node("0", root, grid0, r, {r0, r0}, n_basisfunctions);
    cme_external_node* node1 =
        new cme_external_node("1", root, grid1, r, n_basisfunctions);
    cme_external_node* node00 =
        new cme_external_node("00", node0, grid00, r0, n_basisfunctions0);
    cme_external_node* node01 =
        new cme_external_node("01", node0, grid01, r0, n_basisfunctions0);

    root->Q = Q;
    node0->Q = Q0;
    node1->X = X1;
    node00->X = X00;
    node01->X = X01;

    node1->propensity = propensity1;
    node00->propensity = propensity00;
    node01->propensity = propensity01;

    for (Index mu = 0; mu < grid.n_reactions; ++mu) {
        root->coefficients.A[mu](0, 0) = 1.0;
        root->coefficients.B[mu](0, 0) = 1.0;
    }

    root->child[0] = node0;
    root->child[1] = node1;
    root->child[0]->child[0] = node00;
    root->child[0]->child[1] = node01;

    cme_lr_tree tree(root);

    // Calculate probability distribution
    Ensign::multi_array<double, 3> p({val_n00, val_n01, val_n1}),
        p_ortho({val_n00, val_n01, val_n1});
    std::fill(std::begin(p), std::end(p), 0.0);

    for (Index i = 0; i < r; ++i) {
        for (Index j = 0; j < r; ++j) {
            for (Index i0 = 0; i0 < r0; ++i0) {
                for (Index j0 = 0; j0 < r0; ++j0) {
                    for (Index x00 = 0; x00 < val_n00; ++x00) {
                        for (Index x01 = 0; x01 < val_n01; ++x01) {
                            for (Index x1 = 0; x1 < val_n1; ++x1) {
                                p(x00, x01, x1) += Q(i, j, 0) * Q0(i0, j0, i) *
                                                   X00(x00, i0) * X01(x01, j0) *
                                                   X1(x1, j);
                            }
                        }
                    }
                }
            }
        }
    }

    // Check if the probability distribution remains the same under orthogonalization
    std::fill(std::begin(p_ortho), std::end(p_ortho), 0.0);
    tree.Orthogonalize(blas);

    for (Index i = 0; i < r; ++i) {
        for (Index j = 0; j < r; ++j) {
            for (Index i0 = 0; i0 < r0; ++i0) {
                for (Index j0 = 0; j0 < r0; ++j0) {
                    for (Index x00 = 0; x00 < val_n00; ++x00) {
                        for (Index x01 = 0; x01 < val_n01; ++x01) {
                            for (Index x1 = 0; x1 < val_n1; ++x1) {
                                p_ortho(x00, x01, x1) +=
                                    tree.root->Q(i, j, 0) *
                                    ((cme_internal_node*)tree.root->child[0])
                                        ->Q(i0, j0, i) *
                                    ((cme_external_node*)tree.root->child[0]->child[0])
                                        ->X(x00, i0) *
                                    ((cme_external_node*)tree.root->child[0]->child[1])
                                        ->X(x01, j0) *
                                    ((cme_external_node*)tree.root->child[1])->X(x1, j);
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
    Ensign::Matrix::set_zero(Q0_mat);

    Q(0, 0, 0) = (2.0 * std::exp(-0.5) * std::exp(-0.25) * sqrt(2.0 + std::exp(-4.0)));
    Q0_mat(0, 0) = 1.0;
    Q0_mat(1, 1) = 1.0;
    Q0_mat(2, 2) = 1.0;
    Ensign::Tensor::tensorize<2>(Q0_mat, Q0);

    Ensign::Matrix::set_zero(X00);
    Ensign::Matrix::set_zero(X01);
    Ensign::Matrix::set_zero(X1);

    X00(0, 0) = 1.0;
    X00(1, 0) = 1.0;
    X00(0, 1) = 1.0;
    X00(1, 1) = -1.0;
    X00 /= sqrt(2.0);

    X01(0, 0) = 1.0;
    X01(1, 0) = 1.0;
    X01(0, 1) = 1.0;
    X01(1, 1) = -1.0;
    X01 /= sqrt(2.0);

    X1(0, 0) = 1.0;
    X1(1, 0) = 1.0;
    X1(2, 0) = std::exp(-2.0);
    X1(0, 1) = sqrt(1.0 + 0.5 * std::exp(-4.0));
    X1(1, 1) = -sqrt(1.0 + 0.5 * std::exp(-4.0));
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
    tree.InitializeAB_bar(blas);

    // Check if this yields the same probability distribution
    std::fill(std::begin(p), std::end(p), 0.0);

    for (Index i = 0; i < r; ++i) {
        for (Index j = 0; j < r; ++j) {
            for (Index i0 = 0; i0 < r0; ++i0) {
                for (Index j0 = 0; j0 < r0; ++j0) {
                    for (Index x00 = 0; x00 < val_n00; ++x00) {
                        for (Index x01 = 0; x01 < val_n01; ++x01) {
                            for (Index x1 = 0; x1 < val_n1; ++x1) {
                                p(x00, x01, x1) += Q(i, j, 0) * Q0(i0, j0, i) *
                                                   X00(x00, i0) * X01(x01, j0) *
                                                   X1(x1, j);
                            }
                        }
                    }
                }
            }
        }
    }

    REQUIRE(bool(p == p_ortho));

    Ensign::multi_array<double, 3> G(
        {root->RankOut()[0], root->RankOut()[1], root->RankIn()});
    Ensign::multi_array<double, 2> S0({root->RankOut()[0], root->RankOut()[0]});

    std::fill(std::begin(Q), std::end(Q), 0.0);
    std::fill(std::begin(G), std::end(G), 0.0);
    Ensign::Matrix::set_zero(S0);

    Q(0, 0, 0) = 2.0 * std::exp(-0.75) * sqrt(2.0 + std::exp(-4.0));
    G(0, 0, 0) = 1.0;
    G(1, 1, 0) = 1.0;
    G(2, 2, 0) = 1.0;
    S0(0, 0) = 2.0 * std::exp(-0.75) * sqrt(2.0 + std::exp(-4.0));

    Ensign::multi_array<double, 3> G0(
        {node0->RankOut()[0], node0->RankOut()[1], node0->RankIn()});
    Ensign::multi_array<double, 2> S00({node0->RankOut()[0], node0->RankOut()[0]});

    std::fill(std::begin(Q0), std::end(Q0), 0.0);
    std::fill(std::begin(G0), std::end(G0), 0.0);
    Ensign::Matrix::set_zero(S00);

    Q0(0, 0, 0) = 2.0 * std::exp(-0.75) * sqrt(2.0 + std::exp(-4.0));
    G0(0, 0, 0) = 1.0;
    G0(1, 0, 1) = 1.0;
    S00(0, 0) = 2.0 * std::exp(-0.75) * sqrt(2.0 + std::exp(-4.0));

    root->Q = Q;
    root->G = G;

    node0->Q = Q0;
    node0->G = G0;
    node0->S = S0;

    node00->S = S00;
    node00->X = X00;
    node01->X = X01;
    node1->X = X1;

    root->CalculateAB<0>(blas);
    node0->CalculateAB<0>(blas);

    std::vector<Ensign::multi_array<double, 2>> A0_comparison(node0->grid.n_reactions);
    std::vector<Ensign::multi_array<double, 2>> B0_comparison(node0->grid.n_reactions);
    std::vector<Ensign::multi_array<double, 2>> A00_comparison(
        node00->grid.n_reactions);
    std::vector<Ensign::multi_array<double, 2>> B00_comparison(
        node00->grid.n_reactions);
    std::vector<Ensign::multi_array<double, 2>> A1_bar_comparison(
        node1->grid.n_reactions);
    std::vector<Ensign::multi_array<double, 2>> B1_bar_comparison(
        node1->grid.n_reactions);
    std::vector<Ensign::multi_array<double, 2>> A01_bar_comparison(
        node01->grid.n_reactions);
    std::vector<Ensign::multi_array<double, 2>> B01_bar_comparison(
        node01->grid.n_reactions);

    for (Index mu = 0; mu < root->grid.n_reactions; ++mu) {
        A0_comparison[mu].resize({node0->RankIn(), node0->RankIn()});
        B0_comparison[mu].resize({node0->RankIn(), node0->RankIn()});
        A00_comparison[mu].resize({node00->RankIn(), node00->RankIn()});
        B00_comparison[mu].resize({node00->RankIn(), node00->RankIn()});
        A1_bar_comparison[mu].resize({node1->RankIn(), node1->RankIn()});
        B1_bar_comparison[mu].resize({node1->RankIn(), node1->RankIn()});
        A01_bar_comparison[mu].resize({node01->RankIn(), node01->RankIn()});
        B01_bar_comparison[mu].resize({node01->RankIn(), node01->RankIn()});

        Ensign::Matrix::set_zero(A0_comparison[mu]);
        Ensign::Matrix::set_zero(B0_comparison[mu]);
        Ensign::Matrix::set_zero(A00_comparison[mu]);
        Ensign::Matrix::set_zero(B00_comparison[mu]);
        Ensign::Matrix::set_zero(A1_bar_comparison[mu]);
        Ensign::Matrix::set_zero(B1_bar_comparison[mu]);
        Ensign::Matrix::set_zero(A01_bar_comparison[mu]);
        Ensign::Matrix::set_zero(B01_bar_comparison[mu]);
    }

    // Calculate A1_bar_comparison and B1_bar_comparison
    A1_bar_comparison[0](0, 0) = 1.0;
    A1_bar_comparison[0](1, 1) = 1.0;
    A1_bar_comparison[0](2, 2) = 1.0;
    A1_bar_comparison[1](0, 0) = 1.0;
    A1_bar_comparison[1](1, 1) = 1.0;
    A1_bar_comparison[1](2, 2) = 1.0;
    A1_bar_comparison[3](0, 0) = 1.0;
    A1_bar_comparison[3](1, 1) = 1.0;
    A1_bar_comparison[3](2, 2) = 1.0;
    A1_bar_comparison[4](0, 0) = 1.0;
    A1_bar_comparison[4](1, 1) = 1.0;
    A1_bar_comparison[4](2, 2) = 1.0;

    A1_bar_comparison[2](0, 0) =
        1.0 / (2.0 + std::exp(-4.0)) * (1.0 + 2.0 * std::exp(-2.0));
    A1_bar_comparison[2](0, 1) =
        -1.0 / (2.0 + std::exp(-4.0)) * sqrt(1.0 + 0.5 * std::exp(-4.0));
    A1_bar_comparison[2](0, 2) =
        1.0 / (2.0 + std::exp(-4.0)) * (std::exp(-2.0) - 4.0) / sqrt(2.0);
    A1_bar_comparison[2](1, 0) = 1.0 / (2.0 + std::exp(-4.0)) *
                                 sqrt(1.0 + 0.5 * std::exp(-4.0)) *
                                 (1.0 - 2.0 * std::exp(-2.0));
    A1_bar_comparison[2](1, 1) = -0.5;
    A1_bar_comparison[2](1, 2) = 1.0 / (2.0 + std::exp(-4.0)) *
                                 sqrt(1.0 + 0.5 * std::exp(-4.0)) *
                                 (std::exp(-2.0) + 4.0) / sqrt(2.0);
    A1_bar_comparison[2](2, 0) = 1.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) *
                                 (1.0 + 2.0 * std::exp(-2.0)) / sqrt(2.0);
    A1_bar_comparison[2](2, 1) = -1.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) *
                                 sqrt(1.0 + 0.5 * std::exp(-4.0)) / sqrt(2.0);
    A1_bar_comparison[2](2, 2) =
        1.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) * 0.5 * (std::exp(-2.0) - 4.0);

    A1_bar_comparison[5](0, 0) =
        1.0 / (2.0 + std::exp(-4.0)) * (1.0 + 0.5 * std::exp(-2.0));
    A1_bar_comparison[5](0, 1) = 1.0 / (2.0 + std::exp(-4.0)) *
                                 sqrt(1.0 + 0.5 * std::exp(-4.0)) *
                                 (1.0 - 0.5 * std::exp(-2.0));
    A1_bar_comparison[5](0, 2) = 1.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) *
                                 (1.0 + 0.5 * std::exp(-2.0)) / sqrt(2.0);
    A1_bar_comparison[5](1, 0) =
        -1.0 / (2.0 + std::exp(-4.0)) * sqrt(1.0 + 0.5 * std::exp(-4.0));
    A1_bar_comparison[5](1, 1) = -0.5;
    A1_bar_comparison[5](1, 2) = -1.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) *
                                 sqrt(1.0 + 0.5 * std::exp(-4.0)) / sqrt(2.0);
    A1_bar_comparison[5](2, 0) =
        1.0 / (2.0 + std::exp(-4.0)) * (std::exp(-2.0) - 1.0) / sqrt(2.0);
    A1_bar_comparison[5](2, 1) = 1.0 / (2.0 + std::exp(-4.0)) *
                                 sqrt(1.0 + 0.5 * std::exp(-4.0)) *
                                 (1.0 + std::exp(-2.0)) / sqrt(2.0);
    A1_bar_comparison[5](2, 2) =
        1.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) * 0.5 * (std::exp(-2.0) - 1.0);

    B1_bar_comparison[0](0, 0) = 1.0;
    B1_bar_comparison[0](1, 1) = 1.0;
    B1_bar_comparison[0](2, 2) = 1.0;
    B1_bar_comparison[1](0, 0) = 1.0;
    B1_bar_comparison[1](1, 1) = 1.0;
    B1_bar_comparison[1](2, 2) = 1.0;
    B1_bar_comparison[3](0, 0) = 1.0;
    B1_bar_comparison[3](1, 1) = 1.0;
    B1_bar_comparison[3](2, 2) = 1.0;
    B1_bar_comparison[4](0, 0) = 1.0;
    B1_bar_comparison[4](1, 1) = 1.0;
    B1_bar_comparison[4](2, 2) = 1.0;

    B1_bar_comparison[2](0, 0) =
        1.0 / (2.0 + std::exp(-4.0)) * (1.0 + 2.0 * std::exp(-4.0));
    B1_bar_comparison[2](0, 1) =
        -1.0 / (2.0 + std::exp(-4.0)) * sqrt(1.0 + 0.5 * std::exp(-4.0));
    B1_bar_comparison[2](1, 0) = B1_bar_comparison[2](0, 1);
    B1_bar_comparison[2](0, 2) =
        -3.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) / sqrt(2);
    B1_bar_comparison[2](2, 0) = B1_bar_comparison[2](0, 2);
    B1_bar_comparison[2](1, 1) = 0.5;
    B1_bar_comparison[2](1, 2) = -1.0 / (2.0 + std::exp(-4.0)) *
                                 sqrt(1.0 + 0.5 * std::exp(-4.0)) * std::exp(-2.0) /
                                 sqrt(2.0);
    B1_bar_comparison[2](2, 1) = B1_bar_comparison[2](1, 2);
    B1_bar_comparison[2](2, 2) =
        1.0 / (2.0 + std::exp(-4.0)) * (0.5 * std::exp(-4.0) + 4.0);

    B1_bar_comparison[5](0, 0) =
        1.0 / (2.0 + std::exp(-4.0)) * (1.5 + std::exp(-4.0) / 3.0);
    B1_bar_comparison[5](0, 1) =
        1.0 / (2.0 + std::exp(-4.0)) * 0.5 * sqrt(1.0 + 0.5 * std::exp(-4.0));
    B1_bar_comparison[5](1, 0) = B1_bar_comparison[5](0, 1);
    B1_bar_comparison[5](0, 2) =
        1.0 / (2.0 + std::exp(-4.0)) * 5.0 * std::exp(-2.0) / (6.0 * sqrt(2));
    B1_bar_comparison[5](2, 0) = B1_bar_comparison[5](0, 2);
    B1_bar_comparison[5](1, 1) = 0.75;
    B1_bar_comparison[5](1, 2) = 1.0 / (2.0 + std::exp(-4.0)) * 0.5 *
                                 sqrt(1.0 + 0.5 * std::exp(-4.0)) * std::exp(-2.0) /
                                 sqrt(2.0);
    B1_bar_comparison[5](2, 1) = B1_bar_comparison[5](1, 2);
    B1_bar_comparison[5](2, 2) =
        1.0 / (2.0 + std::exp(-4.0)) * (0.75 * std::exp(-4.0) + 2.0 / 3.0);

    for (Index mu = 0; mu < root->grid.n_reactions; ++mu) {
        REQUIRE(
            bool(tree.root->child[1]->coefficients.A_bar[mu] == A1_bar_comparison[mu]));
        REQUIRE(
            bool(tree.root->child[1]->coefficients.B_bar[mu] == B1_bar_comparison[mu]));
    }

    // Calculate A01_bar_comparison and B01_bar_comparison
    A01_bar_comparison[0](0, 0) = 1.0;
    A01_bar_comparison[0](1, 1) = 1.0;
    A01_bar_comparison[2](0, 0) = 1.0;
    A01_bar_comparison[2](1, 1) = 1.0;
    A01_bar_comparison[3](0, 0) = 1.0;
    A01_bar_comparison[3](1, 1) = 1.0;
    A01_bar_comparison[5](0, 0) = 1.0;
    A01_bar_comparison[5](1, 1) = 1.0;

    A01_bar_comparison[1](0, 0) = 0.5;
    A01_bar_comparison[1](0, 1) = -0.5;
    A01_bar_comparison[1](1, 0) = 0.5;
    A01_bar_comparison[1](1, 1) = -0.5;

    A01_bar_comparison[4](0, 0) = 0.5;
    A01_bar_comparison[4](0, 1) = 0.5;
    A01_bar_comparison[4](1, 0) = -0.5;
    A01_bar_comparison[4](1, 1) = -0.5;

    B01_bar_comparison[0](0, 0) = 1.0;
    B01_bar_comparison[0](1, 1) = 1.0;
    B01_bar_comparison[2](0, 0) = 1.0;
    B01_bar_comparison[2](1, 1) = 1.0;
    B01_bar_comparison[3](0, 0) = 1.0;
    B01_bar_comparison[3](1, 1) = 1.0;
    B01_bar_comparison[5](0, 0) = 1.0;
    B01_bar_comparison[5](1, 1) = 1.0;

    B01_bar_comparison[1](0, 0) = 0.5;
    B01_bar_comparison[1](0, 1) = -0.5;
    B01_bar_comparison[1](1, 0) = B01_bar_comparison[1](0, 1);
    B01_bar_comparison[1](1, 1) = 0.5;

    B01_bar_comparison[4](0, 0) = 0.75;
    B01_bar_comparison[4](0, 1) = 0.25;
    B01_bar_comparison[4](1, 0) = B01_bar_comparison[4](0, 1);
    B01_bar_comparison[4](1, 1) = 0.75;

    for (Index mu = 0; mu < root->grid.n_reactions; ++mu) {
        REQUIRE(bool(tree.root->child[0]->child[1]->coefficients.A_bar[mu] ==
                     A01_bar_comparison[mu]));
        REQUIRE(bool(tree.root->child[0]->child[1]->coefficients.B_bar[mu] ==
                     B01_bar_comparison[mu]));
    }

    // Calculate A0_comparison and B0_comparison
    A0_comparison[0] = A1_bar_comparison[0];
    A0_comparison[1] = A1_bar_comparison[1];
    A0_comparison[2] = A1_bar_comparison[2];
    A0_comparison[3] = A1_bar_comparison[3];
    A0_comparison[4] = A1_bar_comparison[4];
    A0_comparison[5] = A1_bar_comparison[5];

    B0_comparison[0] = B1_bar_comparison[0];
    B0_comparison[1] = B1_bar_comparison[1];
    B0_comparison[2] = B1_bar_comparison[2];
    B0_comparison[3] = B1_bar_comparison[3];
    B0_comparison[4] = B1_bar_comparison[4];
    B0_comparison[5] = B1_bar_comparison[5];

    // Calculate A00_comparison and B00_comparison
    A00_comparison[0](0, 0) = 1.0;
    A00_comparison[0](1, 1) = 1.0;

    A00_comparison[1](0, 0) = 0.5;
    A00_comparison[1](1, 1) = 0.5;

    A00_comparison[2](0, 0) = (1.0 + 2.0 * std::exp(-2.0)) / (2.0 + std::exp(-4.0));
    A00_comparison[2](0, 1) =
        -sqrt(1.0 + 0.5 * std::exp(-4.0)) / (2.0 + std::exp(-4.0));
    A00_comparison[2](1, 0) = sqrt(1.0 + 0.5 * std::exp(-4.0)) *
                              (1.0 - 2.0 * std::exp(-2.0)) / (2.0 + std::exp(-4.0));
    A00_comparison[2](1, 1) = -0.5;

    A00_comparison[3](0, 0) = 1.0;
    A00_comparison[3](1, 1) = 1.0;

    A00_comparison[4](0, 0) = 0.5;
    A00_comparison[4](1, 1) = 0.5;

    A00_comparison[5](0, 0) = (1.0 + 0.5 * std::exp(-2.0)) / (2.0 + std::exp(-4.0));
    A00_comparison[5](0, 1) = sqrt(1.0 + 0.5 * std::exp(-4.0)) *
                              (1.0 - 0.5 * std::exp(-2.0)) / (2.0 + std::exp(-4.0));
    A00_comparison[5](1, 0) =
        -sqrt(1.0 + 0.5 * std::exp(-4.0)) / (2.0 + std::exp(-4.0));
    A00_comparison[5](1, 1) = -0.5;

    B00_comparison[0](0, 0) = 1.0;
    B00_comparison[0](1, 1) = 1.0;

    B00_comparison[1](0, 0) = 0.5;
    B00_comparison[1](1, 1) = 0.5;

    B00_comparison[2](0, 0) = (1.0 + 2.0 * std::exp(-4.0)) / (2.0 + std::exp(-4.0));
    B00_comparison[2](0, 1) =
        -sqrt(1.0 + 0.5 * std::exp(-4.0)) / (2.0 + std::exp(-4.0));
    B00_comparison[2](1, 0) = B00_comparison[2](0, 1);
    B00_comparison[2](1, 1) = 0.5;

    B00_comparison[3](0, 0) = 1.0;
    B00_comparison[3](1, 1) = 1.0;

    B00_comparison[4](0, 0) = 0.75;
    B00_comparison[4](1, 1) = 0.75;

    B00_comparison[5](0, 0) = (1.5 + std::exp(-4.0) / 3.0) / (2.0 + std::exp(-4.0));
    B00_comparison[5](0, 1) =
        0.5 * sqrt(1.0 + 0.5 * std::exp(-4.0)) / (2.0 + std::exp(-4.0));
    B00_comparison[5](1, 0) = B00_comparison[5](0, 1);
    B00_comparison[5](1, 1) = 0.75;

    for (Index mu = 0; mu < root->grid.n_reactions; ++mu) {
        REQUIRE(bool(tree.root->child[0]->coefficients.A[mu] == A0_comparison[mu]));
        REQUIRE(bool(tree.root->child[0]->coefficients.B[mu] == B0_comparison[mu]));
    }

    for (Index mu = 0; mu < root->grid.n_reactions; ++mu) {
        REQUIRE(bool(tree.root->child[0]->child[0]->coefficients.A[mu] ==
                     A00_comparison[mu]));
        REQUIRE(bool(tree.root->child[0]->child[0]->coefficients.B[mu] ==
                     B00_comparison[mu]));
    }

    Ensign::multi_array<double, 2> K00(X00.shape());
    Ensign::multi_array<double, 2> K00_dot(X00.shape());
    Ensign::multi_array<double, 2> K00_dot_comparison(X00.shape());

    blas.matmul(X00, S00, K00);
    K00_dot = CalculateKDot(K00, node00, blas);

    // Calculate K00_comparison
    double alpha = sqrt(2.0) * std::exp(-0.75) * sqrt(2.0 + std::exp(-4.0));
    double beta = (2.0 * (std::exp(-2.0) - std::exp(-4.0)) + 0.5 * std::exp(-2.0) -
                   std::exp(-4.0) / 3.0 - 0.5) /
                      (2.0 + std::exp(-4.0)) -
                  0.25;
    double gamma = (0.5 - 2.0 * std::exp(-2.0)) * sqrt(1.0 + 0.5 * std::exp(-4.0)) /
                   (2.0 + std::exp(-4.0));

    K00_dot_comparison(0, 0) = alpha * beta;
    K00_dot_comparison(0, 1) = alpha * gamma;
    K00_dot_comparison(1, 0) = alpha * (beta - 0.5);
    K00_dot_comparison(1, 1) = alpha * gamma;

    REQUIRE(bool(K00_dot == K00_dot_comparison));
}
