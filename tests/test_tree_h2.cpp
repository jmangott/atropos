#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <algorithm>
#include <cmath>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "integrators.hpp"
#include "matrix.hpp"
#include "tree_class.hpp"

TEST_CASE("tree_h2", "[tree_h2]")
{
    unsigned int substeps = 1;
    std::map<std::string, integration_method *> integrations_methods;
    integrations_methods["K"] = new explicit_euler(substeps);
    integrations_methods["S"] = new explicit_euler(substeps);
    integrations_methods["Q"] = new explicit_euler(substeps);

    blas_ops blas;
    ttn_integrator integrator(blas, integrations_methods);

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

    node1->external_coefficients.propensity = propensity1;
    node00->external_coefficients.propensity = propensity00;
    node01->external_coefficients.propensity = propensity01;

    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        root->coefficients.A[mu](0, 0) = 1.0;
        root->coefficients.B[mu](0, 0) = 1.0;
    }
    root->coefficients.E(0, 0, 0, 0) = 1.0;
    root->coefficients.F(0, 0, 0, 0) = 1.0;

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
    tree.Orthogonalize(blas);

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
    tree.InitializeAB_bar(blas);

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

    gram_schmidt gs(&blas);

    // Calculate G0_comparison (needed for A01_comparison)
    multi_array<double, 2> Qmat({node0->RankIn() * node0->RankOut()[1], node0->RankOut()[0]});
    multi_array<double, 3> Q0_comparison({node0->RankOut()[0], node0->RankOut()[1], node0->RankIn()});
    multi_array<double, 3> G0_comparison({node0->RankOut()[0], node0->RankOut()[1], node0->RankIn()});
    multi_array<double, 2> S00_comparison({node0->RankOut()[0], node0->RankOut()[0]});

    std::fill(std::begin(Q0_comparison), std::end(Q0_comparison), 0.0);

    Q0_comparison(0, 0, 0) = 2.0 * std::exp(-0.75) * sqrt(2.0 + std::exp(-4.0));

    std::function<double(double *, double *)> ip;
    ip = inner_product_from_const_weight(1.0, node0->RankIn() * node0->RankOut()[1]);

    // Compute QR decomposition C^n = (S^(n+id))^T * G^n
    Matrix::Matricize(Q0_comparison, Qmat, 0);
    gs(Qmat, S00_comparison, ip);
    Matrix::Tensorize(Qmat, G0_comparison, 0);

    double tau = 1.0;
    integrator.SubflowPhi<0>(root, tau);

    std::vector<multi_array<double, 2>> A0_comparison(node0->grid.n_reactions);
    std::vector<multi_array<double, 2>> B0_comparison(node0->grid.n_reactions);
    std::vector<multi_array<double, 2>> A00_comparison(node00->grid.n_reactions);
    std::vector<multi_array<double, 2>> B00_comparison(node00->grid.n_reactions);
    std::vector<multi_array<double, 2>> A1_bar_comparison(node1->grid.n_reactions);
    std::vector<multi_array<double, 2>> B1_bar_comparison(node1->grid.n_reactions);
    std::vector<multi_array<double, 2>> A01_bar_comparison(node01->grid.n_reactions);
    std::vector<multi_array<double, 2>> B01_bar_comparison(node01->grid.n_reactions);

    for (Index mu = 0; mu < root->grid.n_reactions; ++mu)
    {
        A0_comparison[mu].resize({node0->RankIn(), node0->RankIn()});
        B0_comparison[mu].resize({node0->RankIn(), node0->RankIn()});
        A00_comparison[mu].resize({node00->RankIn(), node00->RankIn()});
        B00_comparison[mu].resize({node00->RankIn(), node00->RankIn()});
        A1_bar_comparison[mu].resize({node1->RankIn(), node1->RankIn()});
        B1_bar_comparison[mu].resize({node1->RankIn(), node1->RankIn()});
        A01_bar_comparison[mu].resize({node01->RankIn(), node01->RankIn()});
        B01_bar_comparison[mu].resize({node01->RankIn(), node01->RankIn()});

        set_zero(A0_comparison[mu]);
        set_zero(B0_comparison[mu]);
        set_zero(A00_comparison[mu]);
        set_zero(B00_comparison[mu]);
        set_zero(A1_bar_comparison[mu]);
        set_zero(B1_bar_comparison[mu]);
        set_zero(A01_bar_comparison[mu]);
        set_zero(B01_bar_comparison[mu]);
    }

    // Calculate A1_comparison and B1_comparison
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

    A1_bar_comparison[2](0, 0) = 1.0 / (2.0 + std::exp(-4.0)) * (1.0 + 2.0 * std::exp(-2.0));
    A1_bar_comparison[2](0, 1) =-1.0 / (2.0 + std::exp(-4.0)) * sqrt(1.0 + 0.5 * std::exp(-4.0));
    A1_bar_comparison[2](0, 2) = 1.0 / (2.0 + std::exp(-4.0)) * (std::exp(-2.0) - 4.0) / sqrt(2.0);
    A1_bar_comparison[2](1, 0) = 1.0 / (2.0 + std::exp(-4.0)) * sqrt(1.0 + 0.5 * std::exp(-4.0)) * (1.0 - 2.0 * std::exp(-2.0));
    A1_bar_comparison[2](1, 1) =-0.5;
    A1_bar_comparison[2](1, 2) = 1.0 / (2.0 + std::exp(-4.0)) * sqrt(1.0 + 0.5 * std::exp(-4.0)) * (std::exp(-2.0) + 4.0) / sqrt(2.0);
    A1_bar_comparison[2](2, 0) = 1.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) * (1.0 + 2.0 * std::exp(-2.0)) / sqrt(2.0);
    A1_bar_comparison[2](2, 1) =-1.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) * sqrt(1.0 + 0.5 * std::exp(-4.0)) / sqrt(2.0);
    A1_bar_comparison[2](2, 2) = 1.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) * 0.5 * (std::exp(-2.0) - 4.0);

    A1_bar_comparison[5](0, 0) = 1.0 / (2.0 + std::exp(-4.0)) * (1.0 + 0.5 * std::exp(-2.0));
    A1_bar_comparison[5](0, 1) = 1.0 / (2.0 + std::exp(-4.0)) * sqrt(1.0 + 0.5 * std::exp(-4.0)) * (1.0 - 0.5 * std::exp(-2.0));
    A1_bar_comparison[5](0, 2) = 1.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) * (1.0 + 0.5 * std::exp(-2.0)) / sqrt(2.0);
    A1_bar_comparison[5](1, 0) =-1.0 / (2.0 + std::exp(-4.0)) * sqrt(1.0 + 0.5 * std::exp(-4.0));
    A1_bar_comparison[5](1, 1) =-0.5;
    A1_bar_comparison[5](1, 2) =-1.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) * sqrt(1.0 + 0.5 * std::exp(-4.0)) / sqrt(2.0);
    A1_bar_comparison[5](2, 0) = 1.0 / (2.0 + std::exp(-4.0)) * (std::exp(-2.0) - 1.0) / sqrt(2.0);
    A1_bar_comparison[5](2, 1) = 1.0 / (2.0 + std::exp(-4.0)) * sqrt(1.0 + 0.5 * std::exp(-4.0)) * (1.0 + std::exp(-2.0)) / sqrt(2.0);
    A1_bar_comparison[5](2, 2) = 1.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) * 0.5 * (std::exp(-2.0) - 1.0);

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

    B1_bar_comparison[2](0, 0) = 1.0 / (2.0 + std::exp(-4.0)) * (1.0 + 2.0 * std::exp(-4.0));
    B1_bar_comparison[2](0, 1) =-1.0 / (2.0 + std::exp(-4.0)) * sqrt(1.0 + 0.5 * std::exp(-4.0));
    B1_bar_comparison[2](1, 0) = B1_bar_comparison[2](0, 1);
    B1_bar_comparison[2](0, 2) =-3.0 / (2.0 + std::exp(-4.0)) * std::exp(-2.0) / sqrt(2);
    B1_bar_comparison[2](2, 0) = B1_bar_comparison[2](0, 2);
    B1_bar_comparison[2](1, 1) = 0.5;
    B1_bar_comparison[2](1, 2) =-1.0 / (2.0 + std::exp(-4.0)) * sqrt(1.0 + 0.5 * std::exp(-4.0)) * std::exp(-2.0) / sqrt(2.0);
    B1_bar_comparison[2](2, 1) = B1_bar_comparison[2](1, 2);
    B1_bar_comparison[2](2, 2) = 1.0 / (2.0 + std::exp(-4.0)) * (0.5 * std::exp(-4.0) + 4.0);

    B1_bar_comparison[5](0, 0) = 1.0 / (2.0 + std::exp(-4.0)) * (1.5 + std::exp(-4.0) / 3.0);
    B1_bar_comparison[5](0, 1) = 1.0 / (2.0 + std::exp(-4.0)) * 0.5 * sqrt(1.0 + 0.5 * std::exp(-4.0));
    B1_bar_comparison[5](1, 0) = B1_bar_comparison[5](0, 1);
    B1_bar_comparison[5](0, 2) = 1.0 / (2.0 + std::exp(-4.0)) * 5.0 * std::exp(-2.0) / (6.0 * sqrt(2));
    B1_bar_comparison[5](2, 0) = B1_bar_comparison[5](0, 2);
    B1_bar_comparison[5](1, 1) = 0.75;
    B1_bar_comparison[5](1, 2) = 1.0 / (2.0 + std::exp(-4.0)) * 0.5 * sqrt(1.0 + 0.5 * std::exp(-4.0)) * std::exp(-2.0) / sqrt(2.0);
    B1_bar_comparison[5](2, 1) = B1_bar_comparison[5](1, 2);
    B1_bar_comparison[5](2, 2) = 1.0 / (2.0 + std::exp(-4.0)) * (0.75 * std::exp(-4.0) + 2.0 / 3.0);

    for (Index mu = 0; mu < root->grid.n_reactions; ++mu)
    {
        for (Index i0 = 0; i0 < root->RankOut()[0]; ++i0)
        {
            for (Index j0 = 0; j0 < root->RankOut()[0]; ++j0)
            {
                for (Index i1 = 0; i1 < root->RankOut()[1]; ++i1)
                {
                    for (Index j1 = 0; j1 < root->RankOut()[1]; ++j1)
                    {
                        A0_comparison[mu](i0, j0) += root->G(i0, i1, 0) * root->G(j0, j1, 0) * A1_bar_comparison[mu](i1, j1);
                        B0_comparison[mu](i0, j0) += root->G(i0, i1, 0) * root->G(j0, j1, 0) * B1_bar_comparison[mu](i1, j1);
                    }
                }
            }
        }
    }

    // Calculate A01_comparison
    A01_bar_comparison[0](0, 0) = 1.0;
    A01_bar_comparison[0](1, 1) = 1.0;
    A01_bar_comparison[2](0, 0) = 1.0;
    A01_bar_comparison[2](1, 1) = 1.0;
    A01_bar_comparison[3](0, 0) = 1.0;
    A01_bar_comparison[3](1, 1) = 1.0;
    A01_bar_comparison[5](0, 0) = 1.0;
    A01_bar_comparison[5](1, 1) = 1.0;

    A01_bar_comparison[1](0, 0) = 0.5;
    A01_bar_comparison[1](0, 1) =-0.5;
    A01_bar_comparison[1](1, 0) = 0.5;
    A01_bar_comparison[1](1, 1) =-0.5;

    A01_bar_comparison[4](0, 0) = 0.5;
    A01_bar_comparison[4](0, 1) = 0.5;
    A01_bar_comparison[4](1, 0) =-0.5;
    A01_bar_comparison[4](1, 1) =-0.5;

    B01_bar_comparison[0](0, 0) = 1.0;
    B01_bar_comparison[0](1, 1) = 1.0;
    B01_bar_comparison[2](0, 0) = 1.0;
    B01_bar_comparison[2](1, 1) = 1.0;
    B01_bar_comparison[3](0, 0) = 1.0;
    B01_bar_comparison[3](1, 1) = 1.0;
    B01_bar_comparison[5](0, 0) = 1.0;
    B01_bar_comparison[5](1, 1) = 1.0;

    B01_bar_comparison[1](0, 0) = 0.5;
    B01_bar_comparison[1](0, 1) =-0.5;
    B01_bar_comparison[1](1, 0) = B01_bar_comparison[1](0, 1);
    B01_bar_comparison[1](1, 1) = 0.5;

    B01_bar_comparison[4](0, 0) = 0.75;
    B01_bar_comparison[4](0, 1) = 0.25;
    B01_bar_comparison[4](1, 0) = B01_bar_comparison[4](0, 1);
    B01_bar_comparison[4](1, 1) = 0.75;

    for (Index mu = 0; mu < node0->grid.n_reactions; ++mu)
    {
        for (Index i0 = 0; i0 < node0->RankIn(); ++i0)
        {
            for (Index j0 = 0; j0 < node0->RankIn(); ++j0)
            {
                for (Index i00 = 0; i00 < node0->RankOut()[0]; ++i00)
                {
                    for (Index j00 = 0; j00 < node0->RankOut()[0]; ++j00)
                    {
                        for (Index i01 = 0; i01 < node0->RankOut()[1]; ++i01)
                        {
                            for (Index j01 = 0; j01 < node0->RankOut()[1]; ++j01)
                            {
                                A00_comparison[mu](i00, j00) += G0_comparison(i00, i01, i0) * G0_comparison(j00, j01, j0) * A01_bar_comparison[mu](i01, j01) * A0_comparison[mu](i0, j0);
                                B00_comparison[mu](i00, j00) += G0_comparison(i00, i01, i0) * G0_comparison(j00, j01, j0) * B01_bar_comparison[mu](i01, j01) * B0_comparison[mu](i0, j0);
                            }
                        }
                    }
                }
            }
        }
    }

    for (Index mu = 0; mu < root->grid.n_reactions; ++mu)
    {
        REQUIRE(bool(node0->coefficients.A[mu] == A0_comparison[mu]));
        REQUIRE(bool(node0->coefficients.B[mu] == B0_comparison[mu]));
        REQUIRE(bool(node00->coefficients.A[mu] == A00_comparison[mu]));
        REQUIRE(bool(node00->coefficients.B[mu] == B00_comparison[mu]));
    }

    multi_array<double, 2> K00(X00.shape());
    multi_array<double, 2> K00C(X00.shape());
    multi_array<double, 2> K00D(X00.shape());
    multi_array<double, 2> K00C_shift(X00.shape());
    multi_array<double, 2> K_dot(X00.shape());
    set_zero(K_dot);
    blas.matmul(X00, S00_comparison, K00);

    std::vector<Index> vec_index(node00->grid.d);
    for (Index mu = 0; mu < node00->grid.n_reactions; ++mu)
    {
        set_zero(K00C);
        set_zero(K00D);
        std::fill(std::begin(vec_index), std::end(vec_index), 0);
        for (Index x = 0; x < node00->grid.dx; ++x)
        {
            Index x_dep = IndexFunction::VecIndexToDepCombIndex(std::begin(vec_index), std::begin(node00->grid.n_dep[mu]), std::begin(node00->grid.idx_dep[mu]), std::end(node00->grid.idx_dep[mu]));
            for (Index i = 0; i < node00->RankIn(); ++i)
            {
                for (Index j = 0; j < node00->RankIn(); ++j)
                {
                    K00C(x, i) += tau * A00_comparison[mu](i, j) * node00->external_coefficients.propensity[mu][x_dep] * K00(x, j);
                    K00D(x, i) += tau * B00_comparison[mu](i, j) * node00->external_coefficients.propensity[mu][x_dep] * K00(x, j);
                }
            }
            IndexFunction::IncrVecIndex(std::begin(node00->grid.n), std::begin(vec_index), std::end(vec_index));
        }
        Matrix::ShiftRows<1>(K00C_shift, K00C, node00->grid, mu);
        K_dot += K00C_shift;
        K_dot -= K00D;
    }
    K00 += K_dot;

    std::function<double(double *, double *)> ip00;
    ip00 = inner_product_from_const_weight(node00->grid.h_mult, node00->grid.dx);
    gs(K00, node00->S, ip00);

    REQUIRE(bool(K00 == node00->X));
}