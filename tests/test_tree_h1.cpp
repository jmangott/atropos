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

TEST_CASE("tree_h1", "[tree_h1]")
{

    Index r = 2;
    Index n_basisfunctions = 1;
    Index n_reactions = 4;

    Index val_n0 = 2, val_n1 = 2;
    double val_liml0 = 0.0, val_liml1 = 0.0;
    Index val_binsize0 = 1, val_binsize1 = 1;

    vector<Index> n{val_n0, val_n1};
    vector<Index> n0{val_n0};
    vector<Index> n1{val_n1};

    vector<double> liml{val_liml0, val_liml1};
    vector<double> liml0{val_liml0};
    vector<double> liml1{val_liml1};

    vector<Index> binsize{val_binsize0, val_binsize1};
    vector<Index> binsize0{val_binsize0};
    vector<Index> binsize1{val_binsize1};

    multi_array<bool, 2> dep({n_reactions, (Index)n.size()});
    multi_array<bool, 2> dep0({n_reactions, (Index)n0.size()});
    multi_array<bool, 2> dep1({n_reactions, (Index)n1.size()});

    std::fill(std::begin(dep), std::end(dep), false);
    std::fill(std::begin(dep0), std::end(dep0), false);
    std::fill(std::begin(dep1), std::end(dep1), false);

    dep(0, 0) = true;
    dep(1, 1) = true;
    dep(2, 1) = true;
    dep(3, 0) = true;

    dep0(0, 0) = true;
    dep0(3, 0) = true;

    dep1(1, 0) = true;
    dep1(2, 0) = true;

    multi_array<Index, 2> nu({n_reactions, (Index)n.size()});
    multi_array<Index, 2> nu0({n_reactions, (Index)n0.size()});
    multi_array<Index, 2> nu1({n_reactions, (Index)n1.size()});

    std::fill(std::begin(nu), std::end(nu), 0);
    std::fill(std::begin(nu0), std::end(nu0), 0);
    std::fill(std::begin(nu1), std::end(nu1), 0);

    nu(0, 0) = -1;
    nu(1, 1) = -1;
    nu(2, 0) = 1;
    nu(3, 1) = 1;

    nu0(0, 0) = -1;
    nu0(2, 0) = 1;

    nu1(1, 0) = -1;
    nu1(3, 0) = 1;

    grid_parms grid(n, binsize, liml, dep, nu);
    grid_parms grid0(n0, binsize0, liml0, dep0, nu0);
    grid_parms grid1(n1, binsize1, liml1, dep1, nu1);
    grid.Initialize();
    grid0.Initialize();
    grid1.Initialize();

    std::vector<std::vector<double>> propensity(grid.n_reactions);
    std::vector<std::vector<double>> propensity0(grid.n_reactions);
    std::vector<std::vector<double>> propensity1(grid.n_reactions);

    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        propensity[mu].resize(grid.dx_dep[mu]);
        propensity0[mu].resize(grid.dx_dep[mu]);
        propensity1[mu].resize(grid.dx_dep[mu]);
    }

    propensity[0] = {0.0, 1.0};
    propensity[1] = {0.0, 1.0};
    propensity[2] = {1.0, 0.5};
    propensity[3] = {1.0, 0.5};

    propensity0[0] = {0.0, 1.0};
    propensity0[1] = {1.0};
    propensity0[2] = {1.0};
    propensity0[3] = {1.0, 0.5};

    propensity1[0] = {1.0};
    propensity1[1] = {0.0, 1.0};
    propensity1[2] = {1.0, 0.5};
    propensity1[3] = {1.0};

    multi_array<double, 3> Q({r, r, 1});
    multi_array<double, 2> X0({val_n0, r}), X1({val_n1, r});

    // Initialize Q and Q0
    std::fill(std::begin(Q), std::end(Q), 0.0);
    Q(0, 0, 0) = 1.0;

    // Initialize and normalize X0, X1
    set_zero(X0);
    set_zero(X1);

    X0(0, 0) = 1.0;
    X0(1, 0) = 1.0;
    X0 *= std::exp(-0.25);

    X1(0, 0) = 1.0;
    X1(1, 0) = 1.0;
    X1 *= std::exp(-0.25);

    // Construct cme_lr_tree
    cme_internal_node *root = new cme_internal_node("", nullptr, grid, 1, {r, r}, 1);
    cme_external_node *node0 = new cme_external_node("0", root, grid0, r, n_basisfunctions);
    cme_external_node *node1 = new cme_external_node("1", root, grid1, r, n_basisfunctions);

    root->Q = Q;
    node0->X = X0;
    node1->X = X1;

    root->coefficients.propensity = propensity;
    node0->coefficients.propensity = propensity0;
    node1->coefficients.propensity = propensity1;

    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        root->coefficients.A(mu, 0, 0) = 1.0;
        root->coefficients.B(mu, 0, 0) = 1.0;
    }
    root->coefficients.E(0, 0, 0, 0) = 1.0;
    root->coefficients.F(0, 0, 0, 0) = 1.0;

    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        node0->external_coefficients.C[mu].resize({node0->grid.dx_dep[mu], node0->RankIn(), node0->RankIn()});
        node0->external_coefficients.D[mu].resize({node0->grid.dx_dep[mu], node0->RankIn(), node0->RankIn()});

        node1->external_coefficients.C[mu].resize({node1->grid.dx_dep[mu], node1->RankIn(), node1->RankIn()});
        node1->external_coefficients.D[mu].resize({node1->grid.dx_dep[mu], node1->RankIn(), node1->RankIn()});
    }

    root->child[0] = node0;
    root->child[1] = node1;

    cme_lr_tree tree(root);

    // Calculate probability distribution
    multi_array<double, 2> p({val_n0, val_n1}), p_ortho({val_n0, val_n1});
    std::fill(std::begin(p), std::end(p), 0.0);

    for (Index i = 0; i < r; ++i)
    {
        for (Index j = 0; j < r; ++j)
        {
            for (Index x0 = 0; x0 < val_n0; ++x0)
            {
                for (Index x1 = 0; x1 < val_n1; ++x1)
                {
                    p(x0, x1) += Q(i, j, 0) * X0(x0, i) * X1(x1, j);
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
            for (Index x0 = 0; x0 < val_n0; ++x0)
            {
                for (Index x1 = 0; x1 < val_n1; ++x1)
                {
                    p_ortho(x0, x1) += tree.root->Q(i, j, 0) * ((cme_external_node *)tree.root->child[0])->X(x0, i) * ((cme_external_node *)tree.root->child[1])->X(x1, j);
                }
            }
        }
    }

    REQUIRE(bool(p == p_ortho));

    // Orthogonalize Q and X manually
    std::fill(std::begin(Q), std::end(Q), 0.0);

    Q(0, 0, 0) = 2.0 * std::exp(-0.5);

    set_zero(X0);
    set_zero(X1);
    
    X0(0, 0) = 1.0;
    X0(1, 0) = 1.0;
    X0(0, 1) = 1.0;
    X0(1, 1) =-1.0;
    X0 /= sqrt(2.0);

    X1(0, 0) = 1.0;
    X1(1, 0) = 1.0;
    X1(0, 1) = 1.0;
    X1(1, 1) =-1.0;
    X1 /= sqrt(2.0);

    root->Q = Q;
    node0->X = X0;
    node1->X = X1;

    // Check if this yields the same probability distribution
    std::fill(std::begin(p), std::end(p), 0.0);

    for (Index i = 0; i < r; ++i)
    {
        for (Index j = 0; j < r; ++j)
        {
            for (Index x0 = 0; x0 < val_n0; ++x0)
            {
                for (Index x1 = 0; x1 < val_n1; ++x1)
                {
                    p(x0, x1) += Q(i, j, 0) * X0(x0, i) * X1(x1, j);
                }
            }
        }
    }

    REQUIRE(bool(p == p_ortho));

    blas_ops blas;
    gram_schmidt gs(&blas);

    double tau = 1.0;
    SubflowPhi<0>(root, blas, tau);

    multi_array<double, 3> A0_comparison({node0->grid.n_reactions, node0->RankIn(), node0->RankIn()});
    multi_array<double, 3> B0_comparison({node0->grid.n_reactions, node0->RankIn(), node0->RankIn()});
    multi_array<double, 3> A1_bar_comparison({node1->grid.n_reactions, node1->RankIn(), node1->RankIn()});
    multi_array<double, 3> B1_bar_comparison({node1->grid.n_reactions, node1->RankIn(), node1->RankIn()});

    // Calculate A1_comparison
    std::fill(std::begin(A1_bar_comparison), std::end(A1_bar_comparison), 0.0);
    std::fill(std::begin(B1_bar_comparison), std::end(B1_bar_comparison), 0.0);

    A1_bar_comparison(0, 0, 0) = 1.0;
    A1_bar_comparison(0, 1, 1) = 1.0;

    A1_bar_comparison(1, 0, 0) = 0.5;
    A1_bar_comparison(1, 0, 1) =-0.5;
    A1_bar_comparison(1, 1, 0) = 0.5;
    A1_bar_comparison(1, 1, 1) =-0.5;

    A1_bar_comparison(2, 0, 0) = 0.75;
    A1_bar_comparison(2, 0, 1) = 0.25;
    A1_bar_comparison(2, 1, 0) = A1_bar_comparison(2, 0, 1);
    A1_bar_comparison(2, 1, 1) = 0.75;

    A1_bar_comparison(3, 0, 0) = 0.5;
    A1_bar_comparison(3, 0, 1) = 0.5;
    A1_bar_comparison(3, 1, 0) =-0.5;
    A1_bar_comparison(3, 1, 1) =-0.5;

    B1_bar_comparison(0, 0, 0) = 1.0;
    B1_bar_comparison(0, 1, 1) = 1.0;
    B1_bar_comparison(3, 0, 0) = 1.0;
    B1_bar_comparison(3, 1, 1) = 1.0;

    B1_bar_comparison(1, 0, 0) = 0.5;
    B1_bar_comparison(1, 0, 1) =-0.5;
    B1_bar_comparison(1, 1, 0) = B1_bar_comparison(1, 0, 1);
    B1_bar_comparison(1, 1, 1) = 0.5;

    B1_bar_comparison(2, 0, 0) = 0.75;
    B1_bar_comparison(2, 0, 1) = 0.25;
    B1_bar_comparison(2, 1, 0) = B1_bar_comparison(2, 0, 1);
    B1_bar_comparison(2, 1, 1) = 0.75;

    std::fill(std::begin(A0_comparison), std::end(A0_comparison), 0.0);
    std::fill(std::begin(B0_comparison), std::end(B0_comparison), 0.0);

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
                        A0_comparison(mu, i0, j0) += root->G(i0, i1, 0) * root->G(j0, j1, 0) * A1_bar_comparison(mu, i1, j1);
                        B0_comparison(mu, i0, j0) += root->G(i0, i1, 0) * root->G(j0, j1, 0) * B1_bar_comparison(mu, i1, j1);
                    }
                }
            }
        }
    }

    REQUIRE(bool(node0->coefficients.A == A0_comparison));
    REQUIRE(bool(node0->coefficients.B == B0_comparison));

    std::vector<multi_array<double, 3>> C_comparison(node0->grid.n_reactions), D_comparison(node0->grid.n_reactions);
    for (Index mu = 0; mu < node0->grid.n_reactions; ++mu)
    {
        C_comparison[mu].resize({node0->grid.dx_dep[mu], node0->RankIn(), node0->RankIn()});
        D_comparison[mu].resize({node0->grid.dx_dep[mu], node0->RankIn(), node0->RankIn()});
    }

    C_comparison[0](0, 0, 0) = 0.0;
    C_comparison[0](0, 0, 1) = 0.0;
    C_comparison[0](0, 1, 0) = 0.0;
    C_comparison[0](0, 1, 1) = 0.0;

    C_comparison[0](1, 0, 0) = 1.0;
    C_comparison[0](1, 0, 1) = 0.0;
    C_comparison[0](1, 1, 0) = 0.0;
    C_comparison[0](1, 1, 1) = 1.0;

    C_comparison[1](0, 0, 0) = 0.5;
    C_comparison[1](0, 0, 1) = -0.5;
    C_comparison[1](0, 1, 0) = 0.5;
    C_comparison[1](0, 1, 1) = -0.5;

    C_comparison[2](0, 0, 0) = 0.75;
    C_comparison[2](0, 0, 1) = 0.25;
    C_comparison[2](0, 1, 0) = 0.25;
    C_comparison[2](0, 1, 1) = 0.75;

    C_comparison[3](0, 0, 0) = 0.5;
    C_comparison[3](0, 0, 1) = 0.5;
    C_comparison[3](0, 1, 0) = -0.5;
    C_comparison[3](0, 1, 1) = -0.5;

    C_comparison[3](1, 0, 0) = 0.25;
    C_comparison[3](1, 0, 1) = 0.25;
    C_comparison[3](1, 1, 0) = -0.25;
    C_comparison[3](1, 1, 1) = -0.25;

    D_comparison[0](0, 0, 0) = 0.0;
    D_comparison[0](0, 0, 1) = 0.0;
    D_comparison[0](0, 1, 0) = 0.0;
    D_comparison[0](0, 1, 1) = 0.0;

    D_comparison[0](1, 0, 0) = 1.0;
    D_comparison[0](1, 0, 1) = 0.0;
    D_comparison[0](1, 1, 0) = 0.0;
    D_comparison[0](1, 1, 1) = 1.0;

    D_comparison[1](0, 0, 0) = 0.5;
    D_comparison[1](0, 0, 1) = -0.5;
    D_comparison[1](0, 1, 0) = -0.5;
    D_comparison[1](0, 1, 1) = 0.5;

    D_comparison[2](0, 0, 0) = 0.75;
    D_comparison[2](0, 0, 1) = 0.25;
    D_comparison[2](0, 1, 0) = 0.25;
    D_comparison[2](0, 1, 1) = 0.75;

    D_comparison[3](0, 0, 0) = 1.0;
    D_comparison[3](0, 0, 1) = 0.0;
    D_comparison[3](0, 1, 0) = 0.0;
    D_comparison[3](0, 1, 1) = 1.0;

    D_comparison[3](1, 0, 0) = 0.5;
    D_comparison[3](1, 0, 1) = 0.0;
    D_comparison[3](1, 1, 0) = 0.0;
    D_comparison[3](1, 1, 1) = 0.5;

    // Reset S, G and X0 to get reproducable results
    node0->S(0, 0) = 2.0 * std::exp(-0.5);
    node0->S(0, 1) = 0.0;
    node0->S(1, 0) = 0.0;
    node0->S(1, 1) = 0.0;
    root->G(0, 0, 0) = 1.0;
    root->G(1, 0, 0) = 0.0;
    root->G(0, 1, 0) = 0.0;
    root->G(1, 1, 0) = 1.0;
    node0->X = X0;

    // Recalculate the coefficients A, B, C and D
    root->CalculateAB<0>(blas);
    node0->CalculateCD(blas);

    for (Index mu = 0; mu < node0->grid.n_reactions; ++mu)
    {
        REQUIRE(bool(node0->external_coefficients.C[mu] == C_comparison[mu]));
        REQUIRE(bool(node0->external_coefficients.D[mu] == D_comparison[mu]));
    }

    // Test K step
    multi_array<double, 2> K0_comparison(node0->X.shape());
    double norm_2e = std::sqrt(2.0) * std::exp(-0.5);

    K0_comparison(0, 0) = (1.0 - 0.25 * tau) * norm_2e;
    K0_comparison(0, 1) = 0.25 * tau * norm_2e;
    K0_comparison(1, 0) = (1.0 - 1.25 * tau) * norm_2e;
    K0_comparison(1, 1) = 0.75 * tau * norm_2e;

    multi_array<double, 2> tmp(node0->X);
    blas.matmul(tmp, node0->S, node0->X);
    node0->CalculateK(blas, tau);

    REQUIRE(bool(K0_comparison == node0->X));

    std::function<double(double *, double *)> ip0;
    ip0 = inner_product_from_const_weight(node0->grid.h_mult, node0->grid.dx);
    gs(node0->X, node0->S, ip0);

    // Test S step
    multi_array<double, 2> S0_comparison(node0->S.shape());

    // Reset S and X0 to get reproducable results
    node0->X = X0;
    node0->S(0, 0) = 2.0 * std::exp(-0.5);
    node0->S(0, 1) = 0.0;
    node0->S(1, 0) = 0.0;
    node0->S(1, 1) = 0.0;
    multi_array<double, 2> S0_old(node0->S);

    S0_comparison(0, 0) = tau * 1.5 * std::exp(-0.5);
    S0_comparison(0, 1) = -tau * std::exp(-0.5);
    S0_comparison(1, 0) = -tau * std::exp(-0.5);
    S0_comparison(1, 1) = tau * 0.5 * std::exp(-0.5);

    node0->CalculateEF(blas);

    multi_array<double, 4> E_comparison({node0->RankIn(), node0->RankIn(), node0->RankIn(), node0->RankIn()});
    multi_array<double, 4> F_comparison({node0->RankIn(), node0->RankIn(), node0->RankIn(), node0->RankIn()});

    // E_comparison, mu = 0
    E_comparison(0, 0, 0, 0) = 0.5;
    E_comparison(0, 0, 1, 0) = -0.5;
    E_comparison(1, 0, 0, 0) = 0.5;
    E_comparison(1, 0, 1, 0) = -0.5;

    E_comparison(0, 0, 0, 1) = 0.0;
    E_comparison(0, 0, 1, 1) = 0.0;
    E_comparison(1, 0, 0, 1) = 0.0;
    E_comparison(1, 0, 1, 1) = 0.0;

    E_comparison(0, 1, 0, 0) = 0.0;
    E_comparison(0, 1, 1, 0) = 0.0;
    E_comparison(1, 1, 0, 0) = 0.0;
    E_comparison(1, 1, 1, 0) = 0.0;

    E_comparison(0, 1, 0, 1) = 0.5;
    E_comparison(0, 1, 1, 1) = -0.5;
    E_comparison(1, 1, 0, 1) = 0.5;
    E_comparison(1, 1, 1, 1) = -0.5;

    // F_comparison, mu = 0
    F_comparison(0, 0, 0, 0) = 0.5;
    F_comparison(0, 0, 1, 0) = -0.5;
    F_comparison(1, 0, 0, 0) = -0.5;
    F_comparison(1, 0, 1, 0) = 0.5;

    F_comparison(0, 0, 0, 1) = 0.0;
    F_comparison(0, 0, 1, 1) = 0.0;
    F_comparison(1, 0, 0, 1) = 0.0;
    F_comparison(1, 0, 1, 1) = 0.0;

    F_comparison(0, 1, 0, 0) = 0.0;
    F_comparison(0, 1, 1, 0) = 0.0;
    F_comparison(1, 1, 0, 0) = 0.0;
    F_comparison(1, 1, 1, 0) = 0.0;

    F_comparison(0, 1, 0, 1) = 0.5;
    F_comparison(0, 1, 1, 1) = -0.5;
    F_comparison(1, 1, 0, 1) = -0.5;
    F_comparison(1, 1, 1, 1) = 0.5;

    // E_comparison, mu = 1
    E_comparison(0, 0, 0, 0) += 0.5;
    E_comparison(0, 0, 1, 0) += 0.0;
    E_comparison(1, 0, 0, 0) += 0.0;
    E_comparison(1, 0, 1, 0) += 0.5;

    E_comparison(0, 0, 0, 1) += -0.5;
    E_comparison(0, 0, 1, 1) += -0.0;
    E_comparison(1, 0, 0, 1) += -0.0;
    E_comparison(1, 0, 1, 1) += -0.5;

    E_comparison(0, 1, 0, 0) += 0.5;
    E_comparison(0, 1, 1, 0) += 0.0;
    E_comparison(1, 1, 0, 0) += 0.0;
    E_comparison(1, 1, 1, 0) += 0.5;

    E_comparison(0, 1, 0, 1) += -0.5;
    E_comparison(0, 1, 1, 1) += -0.0;
    E_comparison(1, 1, 0, 1) += -0.0;
    E_comparison(1, 1, 1, 1) += -0.5;

    // F_comparison, mu = 1
    F_comparison(0, 0, 0, 0) += 0.5;
    F_comparison(0, 0, 1, 0) += 0.0;
    F_comparison(1, 0, 0, 0) += 0.0;
    F_comparison(1, 0, 1, 0) += 0.5;

    F_comparison(0, 0, 0, 1) += -0.5;
    F_comparison(0, 0, 1, 1) += -0.0;
    F_comparison(1, 0, 0, 1) += -0.0;
    F_comparison(1, 0, 1, 1) += -0.5;

    F_comparison(0, 1, 0, 0) += -0.5;
    F_comparison(0, 1, 1, 0) += -0.0;
    F_comparison(1, 1, 0, 0) += -0.0;
    F_comparison(1, 1, 1, 0) += -0.5;

    F_comparison(0, 1, 0, 1) += 0.5;
    F_comparison(0, 1, 1, 1) += 0.0;
    F_comparison(1, 1, 0, 1) += 0.0;
    F_comparison(1, 1, 1, 1) += 0.5;

    // E_comparison, mu = 2
    E_comparison(0, 0, 0, 0) += 0.375;
    E_comparison(0, 0, 1, 0) += 0.375;
    E_comparison(1, 0, 0, 0) += -0.375;
    E_comparison(1, 0, 1, 0) += -0.375;

    E_comparison(0, 0, 0, 1) += 0.125;
    E_comparison(0, 0, 1, 1) += 0.125;
    E_comparison(1, 0, 0, 1) += -0.125;
    E_comparison(1, 0, 1, 1) += -0.125;

    E_comparison(0, 1, 0, 0) += 0.125;
    E_comparison(0, 1, 1, 0) += 0.125;
    E_comparison(1, 1, 0, 0) += -0.125;
    E_comparison(1, 1, 1, 0) += -0.125;

    E_comparison(0, 1, 0, 1) += 0.375;
    E_comparison(0, 1, 1, 1) += 0.375;
    E_comparison(1, 1, 0, 1) += -0.375;
    E_comparison(1, 1, 1, 1) += -0.375;

    // F_comparison, mu = 2
    F_comparison(0, 0, 0, 0) += 0.75;
    F_comparison(0, 0, 1, 0) += 0.0;
    F_comparison(1, 0, 0, 0) += 0.0;
    F_comparison(1, 0, 1, 0) += 0.75;

    F_comparison(0, 0, 0, 1) += 0.25;
    F_comparison(0, 0, 1, 1) += 0.0;
    F_comparison(1, 0, 0, 1) += 0.0;
    F_comparison(1, 0, 1, 1) += 0.25;

    F_comparison(0, 1, 0, 0) += 0.25;
    F_comparison(0, 1, 1, 0) += 0.0;
    F_comparison(1, 1, 0, 0) += 0.0;
    F_comparison(1, 1, 1, 0) += 0.25;

    F_comparison(0, 1, 0, 1) += 0.75;
    F_comparison(0, 1, 1, 1) += 0.0;
    F_comparison(1, 1, 0, 1) += 0.0;
    F_comparison(1, 1, 1, 1) += 0.75;

    // E_comparison, mu = 3
    E_comparison(0, 0, 0, 0) += 0.375;
    E_comparison(0, 0, 1, 0) += 0.125;
    E_comparison(1, 0, 0, 0) += 0.125;
    E_comparison(1, 0, 1, 0) += 0.375;

    E_comparison(0, 0, 0, 1) += 0.375;
    E_comparison(0, 0, 1, 1) += 0.125;
    E_comparison(1, 0, 0, 1) += 0.125;
    E_comparison(1, 0, 1, 1) += 0.375;

    E_comparison(0, 1, 0, 1) += -0.375;
    E_comparison(0, 1, 1, 1) += -0.125;
    E_comparison(1, 1, 0, 1) += -0.125;
    E_comparison(1, 1, 1, 1) += -0.375;

    E_comparison(0, 1, 0, 0) += -0.375;
    E_comparison(0, 1, 1, 0) += -0.125;
    E_comparison(1, 1, 0, 0) += -0.125;
    E_comparison(1, 1, 1, 0) += -0.375;

    // F_comparison, mu = 3
    F_comparison(0, 0, 0, 0) += 0.75;
    F_comparison(0, 0, 1, 0) += 0.25;
    F_comparison(1, 0, 0, 0) += 0.25;
    F_comparison(1, 0, 1, 0) += 0.75;

    F_comparison(0, 0, 0, 1) += 0.0;
    F_comparison(0, 0, 1, 1) += 0.0;
    F_comparison(1, 0, 0, 1) += 0.0;
    F_comparison(1, 0, 1, 1) += 0.0;

    F_comparison(0, 1, 0, 0) += 0.0;
    F_comparison(0, 1, 1, 0) += 0.0;
    F_comparison(1, 1, 0, 0) += 0.0;
    F_comparison(1, 1, 1, 0) += 0.0;

    F_comparison(0, 1, 0, 1) += 0.75;
    F_comparison(0, 1, 1, 1) += 0.25;
    F_comparison(1, 1, 0, 1) += 0.25;
    F_comparison(1, 1, 1, 1) += 0.75;

    REQUIRE(bool(node0->coefficients.E == E_comparison));
    REQUIRE(bool(node0->coefficients.F == F_comparison));

    node0->CalculateS(blas, tau);
    node0->S -= S0_old;
    REQUIRE(bool(S0_comparison == node0->S));
}