#include "bug_integrator.hpp"

template <Index id>
void bug_integrator::SubflowPhi(cme_internal_node* const node, const double tau) const
{
    Index id_c = (id == 0) ? 1 : 0;

    orthogonalize gs(&blas);
    multi_array<double, 2> Qmat(
        {node->RankIn() * node->RankOut()[id_c], node->RankOut()[id]});
    // std::function<double(double *, double *)> ip;

    // Compute QR decomposition C^n = G^n * (S^(n+id))^T
    gt::start("Mat/Ten");
    Matrix::Matricize<id>(node->Q, Qmat);
    gt::stop("Mat/Ten");
    gt::start("gs");
    gs(Qmat, node->child[id]->S, 1.0);
    gt::stop("gs");
    gt::start("Mat/Ten");
    Matrix::Tensorize<id>(Qmat, node->G);
    gt::stop("Mat/Ten");
    transpose_inplace(node->child[id]->S);

    gt::start("CalculateAB");
    node->CalculateAB<id>(blas);
    gt::stop("CalculateAB");

    if (node->child[id]->IsExternal()) {
        gt::start("External");
        cme_external_node* child_node = (cme_external_node*)node->child[id];

        // Compute K = X * S
        multi_array<double, 2> tmp_x(child_node->X);
        blas.matmul(tmp_x, child_node->S, child_node->X);

        // K step
        const auto K_step_rhs = [child_node, this](const multi_array<double, 2>& K) {
            return CalculateKDot(K, child_node, this->blas);
        };
        gt::start("Integrate K");
        integration_methods.at("K")->integrate(child_node->X, K_step_rhs, tau);
        gt::stop("Integrate K");

        // Perform the QR decomposition K = X * S
        // std::function<double(double *, double *)> ip_x;
        // ip_x = inner_product_from_const_weight(child_node->grid.h_mult,
        // child_node->grid.dx);
        gt::start("gs");
        gs(child_node->X, child_node->S, child_node->grid.h_mult);
        gt::stop("gs");
        gt::stop("External");
    }
    else {
        gt::start("Internal");
        cme_internal_node* child_node = (cme_internal_node*)node->child[id];

        // Set C^(n+i) = Q^(n+id) * S^(n+id)
        multi_array<double, 2> Cmat_child(
            {prod(child_node->RankOut()), child_node->RankIn()});
        multi_array<double, 2> Qmat_child(
            {prod(child_node->RankOut()), child_node->RankIn()});
        gt::start("Mat/Ten");
        Matrix::Matricize<2>(child_node->Q, Qmat_child);
        gt::stop("Mat/Ten");
        set_zero(Cmat_child);
        blas.matmul(Qmat_child, child_node->S, Cmat_child);
        gt::start("Mat/Ten");
        Matrix::Tensorize<2>(Cmat_child, child_node->Q);
        gt::stop("Mat/Ten");
        gt::stop("Internal");

        bug_integrator::operator()(child_node, tau);

        gt::start("Internal");
        // Compute QR decomposition C^(n+id) = Q^(n+id) * S^(n+id)
        // std::function<double(double *, double *)> ip_child;
        // ip_child = inner_product_from_const_weight(1.0, prod(child_node->RankOut()));
        Matrix::Matricize<2>(child_node->Q, Cmat_child);
        gt::start("gs");
        gs(Cmat_child, child_node->S, 1.0);
        gt::stop("gs");
        Matrix::Tensorize<2>(Cmat_child, child_node->Q);
        gt::stop("Internal");
    }
    gt::start("CalculateAB_bar");
    node->child[id]->CalculateAB_bar(blas);
    gt::stop("CalculateAB_bar");

    // Integrate S
    gt::start("S");
    const auto S_step_rhs = [node, this](const multi_array<double, 2>& S) {
        return CalculateSDot(S, node->child[id], this->blas);
    };
    gt::start("Integrate S");
    integration_methods.at("S")->integrate(node->child[id]->S, S_step_rhs, -1.0 * tau);
    gt::stop("Integrate S");

    // Set C^n = G^n * (S^(n+id))^T
    multi_array<double, 2> Gmat(
        {node->RankIn() * node->RankOut()[id_c], node->RankOut()[id]});
    Matrix::Matricize<id>(node->G, Gmat);
    set_zero(Qmat);
    blas.matmul_transb(Gmat, node->child[id]->S, Qmat);
    Matrix::Tensorize<id>(Qmat, node->Q);
    gt::stop("S");
}

template void bug_integrator::SubflowPhi<0>(cme_internal_node* const node,
                                            const double tau) const;

template void bug_integrator::SubflowPhi<1>(cme_internal_node* const node,
                                            const double tau) const;

void bug_integrator::SubflowPsi(cme_internal_node* const node, const double tau) const
{
    multi_array<double, 2> Qmat({prod(node->RankOut()), node->RankIn()});

    gt::start("Mat/Ten");
    Matrix::Matricize<2>(node->Q, Qmat);
    gt::stop("Mat/Ten");

    const auto Q_step_rhs = [node, this](const multi_array<double, 2>& Qmat) {
        return CalculateQDot(Qmat, node, this->blas);
    };
    gt::start("Integrate Q");
    integration_methods.at("Q")->integrate(Qmat, Q_step_rhs, tau);
    gt::stop("Integrate Q");

    gt::start("Mat/Ten");
    Matrix::Tensorize<2>(Qmat, node->Q);
    gt::stop("Mat/Ten");
}

void bug_integrator::SubflowTheta(cme_internal_node* const node) const {}

void bug_integrator::operator()(cme_internal_node* const node, const double tau) const
{
    SubflowPhi<0>(node, tau);
    SubflowPhi<1>(node, tau);
    SubflowPsi(node, tau);
    SubflowTheta(node);
}
