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

TEST_CASE("tree_h1", "[tree_h1]")
{
    Ensign::Matrix::blas_ops blas;

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

    Ensign::multi_array<bool, 2> dep({n_reactions, (Index)n.size()});
    Ensign::multi_array<bool, 2> dep0({n_reactions, (Index)n0.size()});
    Ensign::multi_array<bool, 2> dep1({n_reactions, (Index)n1.size()});

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

    Ensign::multi_array<Index, 2> nu({n_reactions, (Index)n.size()});
    Ensign::multi_array<Index, 2> nu0({n_reactions, (Index)n0.size()});
    Ensign::multi_array<Index, 2> nu1({n_reactions, (Index)n1.size()});

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

    for (Index mu = 0; mu < grid.n_reactions; ++mu) {
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

    Ensign::multi_array<double, 3> Q({r, r, 1});
    Ensign::multi_array<double, 2> X0({val_n0, r}), X1({val_n1, r});

    // Initialize Q and Q0
    std::fill(std::begin(Q), std::end(Q), 0.0);
    Q(0, 0, 0) = 1.0;

    // Initialize and normalize X0, X1
    Ensign::Matrix::set_zero(X0);
    Ensign::Matrix::set_zero(X1);

    X0(0, 0) = 1.0;
    X0(1, 0) = 1.0;
    X0 *= std::exp(-0.25);

    X1(0, 0) = 1.0;
    X1(1, 0) = 1.0;
    X1 *= std::exp(-0.25);

    // Construct cme_lr_tree
    cme_internal_node* root = new cme_internal_node("", nullptr, grid, 1, {r, r}, 1);
    cme_external_node* node0 =
        new cme_external_node("0", root, grid0, r, n_basisfunctions);
    cme_external_node* node1 =
        new cme_external_node("1", root, grid1, r, n_basisfunctions);

    root->Q = Q;
    node0->X = X0;
    node1->X = X1;

    node0->propensity = propensity0;
    node1->propensity = propensity1;

    for (Index mu = 0; mu < grid.n_reactions; ++mu) {
        root->coefficients.A[mu](0, 0) = 1.0;
        root->coefficients.B[mu](0, 0) = 1.0;
    }

    root->child[0] = node0;
    root->child[1] = node1;

    cme_lr_tree tree(root);

    // Calculate probability distribution
    Ensign::multi_array<double, 2> p({val_n0, val_n1}), p_ortho({val_n0, val_n1});
    std::fill(std::begin(p), std::end(p), 0.0);

    for (Index i = 0; i < r; ++i) {
        for (Index j = 0; j < r; ++j) {
            for (Index x0 = 0; x0 < val_n0; ++x0) {
                for (Index x1 = 0; x1 < val_n1; ++x1) {
                    p(x0, x1) += Q(i, j, 0) * X0(x0, i) * X1(x1, j);
                }
            }
        }
    }

    // Check if the probability distribution remains the same under orthogonalization
    std::fill(std::begin(p_ortho), std::end(p_ortho), 0.0);
    tree.Orthogonalize(blas);
    tree.InitializeAB_bar(blas);

    for (Index i = 0; i < r; ++i) {
        for (Index j = 0; j < r; ++j) {
            for (Index x0 = 0; x0 < val_n0; ++x0) {
                for (Index x1 = 0; x1 < val_n1; ++x1) {
                    p_ortho(x0, x1) +=
                        tree.root->Q(i, j, 0) *
                        ((cme_external_node*)tree.root->child[0])->X(x0, i) *
                        ((cme_external_node*)tree.root->child[1])->X(x1, j);
                }
            }
        }
    }

    REQUIRE(bool(p == p_ortho));

    // Orthogonalize Q and X manually
    std::fill(std::begin(Q), std::end(Q), 0.0);

    Q(0, 0, 0) = 2.0 * std::exp(-0.5);

    Ensign::Matrix::set_zero(X0);
    Ensign::Matrix::set_zero(X1);

    X0(0, 0) = 1.0;
    X0(1, 0) = 1.0;
    X0(0, 1) = 1.0;
    X0(1, 1) = -1.0;
    X0 /= sqrt(2.0);

    X1(0, 0) = 1.0;
    X1(1, 0) = 1.0;
    X1(0, 1) = 1.0;
    X1(1, 1) = -1.0;
    X1 /= sqrt(2.0);

    root->Q = Q;
    node0->X = X0;
    node1->X = X1;

    // Check if this yields the same probability distribution
    std::fill(std::begin(p), std::end(p), 0.0);

    for (Index i = 0; i < r; ++i) {
        for (Index j = 0; j < r; ++j) {
            for (Index x0 = 0; x0 < val_n0; ++x0) {
                for (Index x1 = 0; x1 < val_n1; ++x1) {
                    p(x0, x1) += Q(i, j, 0) * X0(x0, i) * X1(x1, j);
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

    Q(0, 0, 0) = 2.0 * std::exp(-0.5);
    G(0, 0, 0) = 1.0;
    G(1, 1, 0) = 1.0;
    S0(0, 0) = 2.0 * std::exp(-0.5);

    root->Q = Q;
    root->G = G;
    node0->S = S0;

    tree.root->CalculateAB<0>(blas);

    std::vector<Ensign::multi_array<double, 2>> A0_comparison(node0->grid.n_reactions);
    std::vector<Ensign::multi_array<double, 2>> B0_comparison(node0->grid.n_reactions);
    std::vector<Ensign::multi_array<double, 2>> A1_bar_comparison(
        node1->grid.n_reactions);
    std::vector<Ensign::multi_array<double, 2>> B1_bar_comparison(
        node1->grid.n_reactions);

    // Calculate A1_comparison
    for (Index mu = 0; mu < root->grid.n_reactions; ++mu) {
        A0_comparison[mu].resize({node0->RankIn(), node0->RankIn()});
        B0_comparison[mu].resize({node0->RankIn(), node0->RankIn()});
        A1_bar_comparison[mu].resize({node1->RankIn(), node1->RankIn()});
        B1_bar_comparison[mu].resize({node1->RankIn(), node1->RankIn()});

        Ensign::Matrix::set_zero(A0_comparison[mu]);
        Ensign::Matrix::set_zero(B0_comparison[mu]);
        Ensign::Matrix::set_zero(A1_bar_comparison[mu]);
        Ensign::Matrix::set_zero(B1_bar_comparison[mu]);
    }

    A1_bar_comparison[0](0, 0) = 1.0;
    A1_bar_comparison[0](1, 1) = 1.0;

    A1_bar_comparison[1](0, 0) = 0.5;
    A1_bar_comparison[1](0, 1) = -0.5;
    A1_bar_comparison[1](1, 0) = 0.5;
    A1_bar_comparison[1](1, 1) = -0.5;

    A1_bar_comparison[2](0, 0) = 0.75;
    A1_bar_comparison[2](0, 1) = 0.25;
    A1_bar_comparison[2](1, 0) = A1_bar_comparison[2](0, 1);
    A1_bar_comparison[2](1, 1) = 0.75;

    A1_bar_comparison[3](0, 0) = 0.5;
    A1_bar_comparison[3](0, 1) = 0.5;
    A1_bar_comparison[3](1, 0) = -0.5;
    A1_bar_comparison[3](1, 1) = -0.5;

    B1_bar_comparison[0](0, 0) = 1.0;
    B1_bar_comparison[0](1, 1) = 1.0;
    B1_bar_comparison[3](0, 0) = 1.0;
    B1_bar_comparison[3](1, 1) = 1.0;

    B1_bar_comparison[1](0, 0) = 0.5;
    B1_bar_comparison[1](0, 1) = -0.5;
    B1_bar_comparison[1](1, 0) = B1_bar_comparison[1](0, 1);
    B1_bar_comparison[1](1, 1) = 0.5;

    B1_bar_comparison[2](0, 0) = 0.75;
    B1_bar_comparison[2](0, 1) = 0.25;
    B1_bar_comparison[2](1, 0) = B1_bar_comparison[2](0, 1);
    B1_bar_comparison[2](1, 1) = 0.75;

    for (Index mu = 0; mu < root->grid.n_reactions; ++mu) {
        for (Index i0 = 0; i0 < root->RankOut()[0]; ++i0) {
            for (Index j0 = 0; j0 < root->RankOut()[0]; ++j0) {
                for (Index i1 = 0; i1 < root->RankOut()[1]; ++i1) {
                    for (Index j1 = 0; j1 < root->RankOut()[1]; ++j1) {
                        A0_comparison[mu](i0, j0) += root->G(i0, i1, 0) *
                                                     root->G(j0, j1, 0) *
                                                     A1_bar_comparison[mu](i1, j1);
                        B0_comparison[mu](i0, j0) += root->G(i0, i1, 0) *
                                                     root->G(j0, j1, 0) *
                                                     B1_bar_comparison[mu](i1, j1);
                    }
                }
            }
        }
    }

    for (Index mu = 0; mu < root->grid.n_reactions; ++mu) {
        REQUIRE(bool(node0->coefficients.A[mu] == A0_comparison[mu]));
        REQUIRE(bool(node0->coefficients.B[mu] == B0_comparison[mu]));
    }

    // Test K step
    Ensign::multi_array<double, 2> K0_dot_comparison(node0->X.shape());
    double norm_2e = std::sqrt(2.0) * std::exp(-0.5);

    K0_dot_comparison(0, 0) = -0.25 * norm_2e;
    K0_dot_comparison(0, 1) = 0.25 * norm_2e;
    K0_dot_comparison(1, 0) = -1.25 * norm_2e;
    K0_dot_comparison(1, 1) = 0.75 * norm_2e;

    Ensign::multi_array<double, 2> K0(node0->X.shape());
    Ensign::multi_array<double, 2> K0_dot(node0->X.shape());

    blas.matmul(node0->X, node0->S, K0);
    K0_dot = CalculateKDot(K0, node0, blas);

    REQUIRE(bool(K0_dot_comparison == K0_dot));

    // Test S step
    Ensign::multi_array<double, 2> S0_dot(node0->S.shape());
    Ensign::multi_array<double, 2> S0_dot_comparison(node0->S.shape());

    // Reset S and X0 to get reproducable results
    node0->X = X0;
    node0->S = S0;

    S0_dot_comparison(0, 0) = -1.5 * std::exp(-0.5);
    S0_dot_comparison(0, 1) = std::exp(-0.5);
    S0_dot_comparison(1, 0) = std::exp(-0.5);
    S0_dot_comparison(1, 1) = -0.5 * std::exp(-0.5);

    S0_dot = CalculateSDot(node0->S, node0, blas);

    REQUIRE(bool(S0_dot_comparison == S0_dot));

    // Test the relation S1_dot(i1, j1) = Q_dot(i0, i1, i) * G(i0, j1, i)
    Ensign::multi_array<double, 2> S1_dot(
        {tree.root->RankOut()[1], tree.root->RankOut()[1]});
    Ensign::multi_array<double, 3> Q_dot(root->Q.shape());
    Ensign::multi_array<double, 2> Qmat1(
        {root->RankIn() * root->RankOut()[0], root->RankOut()[1]});
    Ensign::multi_array<double, 2> Qmat2(
        {Ensign::prod(root->RankOut()), root->RankIn()});
    Ensign::multi_array<double, 2> Qmat2_dot(Qmat2.shape());
    Ensign::multi_array<double, 2> Gmat1(Qmat1.shape());
    Ensign::Matrix::set_zero(Qmat1);

    Ensign::multi_array<double, 2> S1(
        {tree.root->RankOut()[1], tree.root->RankOut()[1]});
    Ensign::Matrix::set_zero(S1);
    S1(0, 0) = 2.0 * std::exp(-0.5);
    node1->S = S1;

    root->CalculateAB<1>(blas);
    S1_dot = CalculateSDot(node1->S, node1, blas);

    Ensign::Tensor::matricize<1>(root->G, Gmat1);
    blas.matmul_transb(Gmat1, node1->S, Qmat1);
    Ensign::Tensor::tensorize<1>(Qmat1, root->Q);

    Ensign::Tensor::matricize<2>(root->Q, Qmat2);
    Qmat2_dot = CalculateQDot(Qmat2, root, blas);
    Ensign::Tensor::tensorize<2>(Qmat2_dot, Q_dot);

    Ensign::multi_array<double, 2> Q_dotG(S1_dot.shape());
    Ensign::Matrix::set_zero(Q_dotG);

    for (Index i = 0; i < root->RankIn(); ++i) {
        for (Index i0 = 0; i0 < root->RankOut()[0]; ++i0) {
            for (Index i1 = 0; i1 < root->RankOut()[1]; ++i1) {
                for (Index j1 = 0; j1 < root->RankOut()[1]; ++j1) {
                    Q_dotG(i1, j1) += Q_dot(i0, i1, i) * root->G(i0, j1, i);
                }
            }
        }
    }

    REQUIRE(bool(Q_dotG == S1_dot));
}
