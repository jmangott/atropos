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

TEST_CASE("tree_h1", "[tree_h1]")
{
    unsigned int substeps = 1;
    std::map<std::string, integration_method *> integrations_methods;
    integrations_methods["K"] = new explicit_euler(substeps);
    integrations_methods["S"] = new explicit_euler(substeps);
    integrations_methods["Q"] = new explicit_euler(substeps);

    blas_ops blas;
    ttn_integrator integrator(blas, integrations_methods);

    Index r = 2;
    Index n_basisfunctions = 1;
    Index n_reactions = 4;

    Index val_n0 = 2, val_n1 = 2;
    double val_liml0 = 0.0, val_liml1 = 0.0;
    Index val_binsize0 = 1, val_binsize1 = 1;
    int val_species0 = 0, val_species1 = 1;

    std::vector<Index> n{val_n0, val_n1};
    std::vector<Index> n0{val_n0};
    std::vector<Index> n1{val_n1};

    std::vector<double> liml{val_liml0, val_liml1};
    std::vector<double> liml0{val_liml0};
    std::vector<double> liml1{val_liml1};

    std::vector<Index> binsize{val_binsize0, val_binsize1};
    std::vector<Index> binsize0{val_binsize0};
    std::vector<Index> binsize1{val_binsize1};

    std::vector<int> species{val_species0, val_species1};
    std::vector<int> species0{val_species0};
    std::vector<int> species1{val_species1};

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

    grid_parms grid(n, binsize, liml, dep, nu, species);
    grid_parms grid0(n0, binsize0, liml0, dep0, nu0, species0);
    grid_parms grid1(n1, binsize1, liml1, dep1, nu1, species1);
    grid.Initialize();
    grid0.Initialize();
    grid1.Initialize();

    std::vector<std::vector<double>> propensity0(grid.n_reactions);
    std::vector<std::vector<double>> propensity1(grid.n_reactions);

    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        propensity0[mu].resize(grid.dx_dep[mu]);
        propensity1[mu].resize(grid.dx_dep[mu]);
    }

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

    node0->external_coefficients.propensity = propensity0;
    node1->external_coefficients.propensity = propensity1;

    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        root->coefficients.A[mu](0, 0) = 1.0;
        root->coefficients.B[mu](0, 0) = 1.0;
    }
    root->coefficients.E(0, 0, 0, 0) = 1.0;
    root->coefficients.F(0, 0, 0, 0) = 1.0;

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
    tree.Orthogonalize(blas);
    tree.InitializeAB_bar(blas);

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

    gram_schmidt gs(&blas);

    double tau = 1.0;
    integrator.SubflowPhi<0>(root, tau);

    std::vector<multi_array<double, 2>> A0_comparison(node0->grid.n_reactions);
    std::vector<multi_array<double, 2>> B0_comparison(node0->grid.n_reactions);
    std::vector<multi_array<double, 2>> A1_bar_comparison(node1->grid.n_reactions);
    std::vector<multi_array<double, 2>> B1_bar_comparison(node1->grid.n_reactions);

    // Calculate A1_comparison
    for (Index mu = 0; mu < root->grid.n_reactions; ++mu)
    {
        A0_comparison[mu].resize({node0->RankIn(), node0->RankIn()});
        B0_comparison[mu].resize({node0->RankIn(), node0->RankIn()});
        A1_bar_comparison[mu].resize({node1->RankIn(), node1->RankIn()});
        B1_bar_comparison[mu].resize({node1->RankIn(), node1->RankIn()});

        set_zero(A0_comparison[mu]);
        set_zero(B0_comparison[mu]);
        set_zero(A1_bar_comparison[mu]);
        set_zero(B1_bar_comparison[mu]);
    }

    A1_bar_comparison[0](0, 0) = 1.0;
    A1_bar_comparison[0](1, 1) = 1.0;

    A1_bar_comparison[1](0, 0) = 0.5;
    A1_bar_comparison[1](0, 1) =-0.5;
    A1_bar_comparison[1](1, 0) = 0.5;
    A1_bar_comparison[1](1, 1) =-0.5;

    A1_bar_comparison[2](0, 0) = 0.75;
    A1_bar_comparison[2](0, 1) = 0.25;
    A1_bar_comparison[2](1, 0) = A1_bar_comparison[2](0, 1);
    A1_bar_comparison[2](1, 1) = 0.75;

    A1_bar_comparison[3](0, 0) = 0.5;
    A1_bar_comparison[3](0, 1) = 0.5;
    A1_bar_comparison[3](1, 0) =-0.5;
    A1_bar_comparison[3](1, 1) =-0.5;

    B1_bar_comparison[0](0, 0) = 1.0;
    B1_bar_comparison[0](1, 1) = 1.0;
    B1_bar_comparison[3](0, 0) = 1.0;
    B1_bar_comparison[3](1, 1) = 1.0;

    B1_bar_comparison[1](0, 0) = 0.5;
    B1_bar_comparison[1](0, 1) =-0.5;
    B1_bar_comparison[1](1, 0) = B1_bar_comparison[1](0, 1);
    B1_bar_comparison[1](1, 1) = 0.5;

    B1_bar_comparison[2](0, 0) = 0.75;
    B1_bar_comparison[2](0, 1) = 0.25;
    B1_bar_comparison[2](1, 0) = B1_bar_comparison[2](0, 1);
    B1_bar_comparison[2](1, 1) = 0.75;

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

    for (Index mu = 0; mu < root->grid.n_reactions; ++mu)
    {
        REQUIRE(bool(node0->coefficients.A[mu] == A0_comparison[mu]));
        REQUIRE(bool(node0->coefficients.B[mu] == B0_comparison[mu]));
    }

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

    // Recalculate the coefficients A and B
    root->CalculateAB<0>(blas);

    // Test K step
    multi_array<double, 2> K0_comparison(node0->X.shape());
    double norm_2e = std::sqrt(2.0) * std::exp(-0.5);

    K0_comparison(0, 0) = (1.0 - 0.25 * tau) * norm_2e;
    K0_comparison(0, 1) = 0.25 * tau * norm_2e;
    K0_comparison(1, 0) = (1.0 - 1.25 * tau) * norm_2e;
    K0_comparison(1, 1) = 0.75 * tau * norm_2e;

    multi_array<double, 2> tmp(node0->X);
    blas.matmul(tmp, node0->S, node0->X);
    const auto K_step_rhs0 = [node0, blas] (const multi_array<double, 2> &K) { return CalculateKDot(K, node0, blas); };
    integrator.integration_methods.at("K")->integrate(node0->X, K_step_rhs0, tau);

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
    node0->CalculateAB_bar(blas);
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

    const auto S_step_rhs0 = [node0](const multi_array<double, 2> &S)
    { return CalculateSDot(S, node0); };
    integrator.integration_methods.at("S")->integrate(node0->S, S_step_rhs0, -1.0 * tau);
    node0->S -= S0_old;
    REQUIRE(bool(S0_comparison == node0->S));

    root->Q = Q;
    node0->X = X0;
    node1->X = X1;

    integrator(tree.root, tau);

    std::fill(std::begin(E_comparison), std::end(E_comparison), 0.0);
    std::fill(std::begin(F_comparison), std::end(F_comparison), 0.0);

    // Test the relation G(i0, k1, i) * G(j0, l1, j) * g(i, i0, i1, j, j0, j1) = e(i1, k1, j1, l1)
    for (Index i = 0; i < root->RankIn(); ++i)
    {
        for (Index j = 0; j < root->RankIn(); ++j)
        {
            for (Index i0 = 0; i0 < root->RankOut()[0]; ++i0)
            {
                for (Index j0 = 0; j0 < root->RankOut()[0]; ++j0)
                {
                    for (Index i1 = 0; i1 < root->RankOut()[1]; ++i1)
                    {
                        for (Index j1 = 0; j1 < root->RankOut()[1]; ++j1)
                        {
                            for (Index k1 = 0; k1 < root->RankOut()[1]; ++k1)
                            {
                                for (Index l1 = 0; l1 < root->RankOut()[1]; ++l1)
                                {
                                    E_comparison(i1, k1, j1, l1) += root->internal_coefficients.G(i, i0 + root->RankOut()[0] * i1, j, j0 + root->RankOut()[0] * j1) * root->G(i0, k1, i) * root->G(j0, l1, j);
                                    F_comparison(i1, k1, j1, l1) += root->internal_coefficients.H(i, i0 + root->RankOut()[0] * i1, j, j0 + root->RankOut()[0] * j1) * root->G(i0, k1, i) * root->G(j0, l1, j);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    REQUIRE(bool(node1->coefficients.E == E_comparison));
    REQUIRE(bool(node1->coefficients.F == F_comparison));

    // Test the relation S_dot(i1, j1) = Q_dot(i0, i1, i) * G(i0, j1, i)
    root->Q = Q;
    node0->X = X0;
    node1->X = X1;
    tau = 0.001;
    node0->CalculateAB_bar(blas);
    node1->CalculateAB_bar(blas);

    integrator.SubflowPhi<0>(tree.root, tau);

    // Perform the L step manually
    multi_array<double, 2> Qmat({root->RankIn() * root->RankOut()[0], root->RankOut()[1]});
    std::function<double(double *, double *)> ip;
    ip = inner_product_from_const_weight(1.0, root->RankIn() * root->RankOut()[0]);

    // Compute QR decomposition C^n = G^n * (S^(n+id))^T
    Matrix::Matricize(root->Q, Qmat, 1);
    gs(Qmat, root->child[1]->S, ip);
    Matrix::Tensorize(Qmat, root->G, 1);
    transpose_inplace(root->child[1]->S);

    root->CalculateAB<1>(blas);

    // Compute K = X * S
    multi_array<double, 2> tmp_x(node1->X);
    blas.matmul(tmp_x, node1->S, node1->X);

    // K step
    const auto K_step_rhs1 = [node1, blas](const multi_array<double, 2> &K) { return CalculateKDot(K, node1, blas); };
    integrator.integration_methods.at("K")->integrate(node1->X, K_step_rhs1, tau);

    // Perform the QR decomposition K = X * S
    std::function<double(double *, double *)> ip_x;
    ip_x = inner_product_from_const_weight(node1->grid.h_mult, node1->grid.dx);
    gs(node1->X, node1->S, ip_x);

    node1->CalculateAB_bar(blas);

    // Integrate S
    multi_array<double, 2> deltaS(node1->S);

    node1->CalculateEF(blas);
    const auto S_step_rhs1 = [node1](const multi_array<double, 2> &S)
    { return CalculateSDot(S, node1); };
    integrator.integration_methods.at("S")->integrate(node1->S, S_step_rhs1, -1.0 * tau);

    // Set C^n = G^n * (S^(n+id))^T
    multi_array<double, 2> Gmat({root->RankIn() * root->RankOut()[0], root->RankOut()[1]});
    Matrix::Matricize(root->G, Gmat, 1);
    set_zero(Qmat);
    blas.matmul_transb(Gmat, deltaS, Qmat);
    Matrix::Tensorize(Qmat, root->Q, 1);

    multi_array<double, 3> deltaQ(root->Q);

    integrator.SubflowPsi(tree.root, tau);
    deltaS -= node1->S;
    deltaQ -= root->Q;

    multi_array<double, 2> deltaQG({node1->RankIn(), node1->RankIn()});
    set_zero(deltaQG);

    for (Index i = 0; i < root->RankIn(); ++i)
    {
        for (Index i0 = 0; i0 < root->RankOut()[0]; ++i0)
        {
            for (Index i1 = 0; i1 < root->RankOut()[1]; ++i1)
            {
                for (Index j1 = 0; j1 < root->RankOut()[1]; ++j1)
                {
                    deltaQG(i1, j1) -= deltaQ(i0, i1, i) * root->G(i0, j1, i);
                }
            }
        }
    }

    REQUIRE(bool(deltaQG == deltaS));
}