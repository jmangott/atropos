#include "bug_integrator.hpp"

template <Index id>
void bug_integrator::SubflowPhi(cme_internal_node* const node, const double tau) const
{
    Index id_c = (id == 0) ? 1 : 0;

    Ensign::orthogonalize gs(&blas);
    Ensign::multi_array<double, 2> Qmat(
        {node->RankIn() * node->RankOut()[id_c], node->RankOut()[id]});
    // std::function<double(double *, double *)> ip;

    // Compute QR decomposition C^n = G^n * (S^(n+id))^T
    Ensign::gt::start("Mat/Ten");
    Ensign::Tensor::matricize<id>(node->Q, Qmat);
    Ensign::gt::stop("Mat/Ten");
    Ensign::gt::start("gs");
    gs(Qmat, node->child[id]->S, 1.0);
    Ensign::gt::stop("gs");
    Ensign::gt::start("Mat/Ten");
    Ensign::Tensor::tensorize<id>(Qmat, node->G);
    Ensign::gt::stop("Mat/Ten");
    Ensign::Matrix::transpose_inplace(node->child[id]->S);

    Ensign::gt::start("CalculateAB");
    node->CalculateAB<id>(blas);
    Ensign::gt::stop("CalculateAB");

    if (node->child[id]->IsExternal()) {
        Ensign::gt::start("External");
        cme_external_node* child_node = (cme_external_node*)node->child[id];

        // Compute K = X * S
        Ensign::multi_array<double, 2> tmp_x(child_node->X);
        blas.matmul(tmp_x, child_node->S, child_node->X);

        // K step
        const auto K_step_rhs = [child_node,
                                 this](const Ensign::multi_array<double, 2>& K) {
            return CalculateKDot(K, child_node, this->blas);
        };
        Ensign::gt::start("Integrate K");
        integration_methods.at("K")->integrate(child_node->X, K_step_rhs, tau);
        Ensign::gt::stop("Integrate K");

        // Perform the QR decomposition K = X * S
        // std::function<double(double *, double *)> ip_x;
        // ip_x = Ensign::inner_product_from_const_weight(child_node->grid.h_mult,
        // child_node->grid.dx);
        Ensign::gt::start("gs");
        gs(child_node->X, child_node->S, child_node->grid.h_mult);
        Ensign::gt::stop("gs");
        Ensign::gt::stop("External");
    }
    else {
        Ensign::gt::start("Internal");
        cme_internal_node* child_node = (cme_internal_node*)node->child[id];

        // Set C^(n+i) = Q^(n+id) * S^(n+id)
        Ensign::multi_array<double, 2> Cmat_child(
            {Ensign::prod(child_node->RankOut()), child_node->RankIn()});
        Ensign::multi_array<double, 2> Qmat_child(
            {Ensign::prod(child_node->RankOut()), child_node->RankIn()});
        Ensign::gt::start("Mat/Ten");
        Ensign::Tensor::matricize<2>(child_node->Q, Qmat_child);
        Ensign::gt::stop("Mat/Ten");
        Ensign::Matrix::set_zero(Cmat_child);
        blas.matmul(Qmat_child, child_node->S, Cmat_child);
        Ensign::gt::start("Mat/Ten");
        Ensign::Tensor::tensorize<2>(Cmat_child, child_node->Q);
        Ensign::gt::stop("Mat/Ten");
        Ensign::gt::stop("Internal");

        bug_integrator::operator()(child_node, tau);

        Ensign::gt::start("Internal");
        // Compute QR decomposition C^(n+id) = Q^(n+id) * S^(n+id)
        // std::function<double(double *, double *)> ip_child;
        // ip_child = Ensign::inner_product_from_const_weight(1.0,
        // Ensign::prod(child_node->RankOut()));
        Ensign::Tensor::matricize<2>(child_node->Q, Cmat_child);
        Ensign::gt::start("gs");
        gs(Cmat_child, child_node->S, 1.0);
        Ensign::gt::stop("gs");
        Ensign::Tensor::tensorize<2>(Cmat_child, child_node->Q);
        Ensign::gt::stop("Internal");
    }
    Ensign::gt::start("CalculateAB_bar");
    node->child[id]->CalculateAB_bar(blas);
    Ensign::gt::stop("CalculateAB_bar");

    // Integrate S
    Ensign::gt::start("S");
    const auto S_step_rhs = [node, this](const Ensign::multi_array<double, 2>& S) {
        return CalculateSDot(S, node->child[id], this->blas);
    };
    Ensign::gt::start("Integrate S");
    integration_methods.at("S")->integrate(node->child[id]->S, S_step_rhs, -1.0 * tau);
    Ensign::gt::stop("Integrate S");

    // Set C^n = G^n * (S^(n+id))^T
    Ensign::multi_array<double, 2> Gmat(
        {node->RankIn() * node->RankOut()[id_c], node->RankOut()[id]});
    Ensign::Tensor::matricize<id>(node->G, Gmat);
    Ensign::Matrix::set_zero(Qmat);
    blas.matmul_transb(Gmat, node->child[id]->S, Qmat);
    Ensign::Tensor::tensorize<id>(Qmat, node->Q);
    Ensign::gt::stop("S");
}

template void bug_integrator::SubflowPhi<0>(cme_internal_node* const node,
                                            const double tau) const;

template void bug_integrator::SubflowPhi<1>(cme_internal_node* const node,
                                            const double tau) const;

void bug_integrator::SubflowPsi(cme_internal_node* const node, const double tau) const
{
    Ensign::multi_array<double, 2> Qmat(
        {Ensign::prod(node->RankOut()), node->RankIn()});

    Ensign::gt::start("Mat/Ten");
    Ensign::Tensor::matricize<2>(node->Q, Qmat);
    Ensign::gt::stop("Mat/Ten");

    const auto Q_step_rhs = [node, this](const Ensign::multi_array<double, 2>& Qmat) {
        return CalculateQDot(Qmat, node, this->blas);
    };
    Ensign::gt::start("Integrate Q");
    integration_methods.at("Q")->integrate(Qmat, Q_step_rhs, tau);
    Ensign::gt::stop("Integrate Q");

    Ensign::gt::start("Mat/Ten");
    Ensign::Tensor::tensorize<2>(Qmat, node->Q);
    Ensign::gt::stop("Mat/Ten");
}

void bug_integrator::SubflowTheta(cme_internal_node* const node) const {}

void bug_integrator::operator()(cme_internal_node* const node, const double tau) const
{
    SubflowPhi<0>(node, tau);
    SubflowPhi<1>(node, tau);
    SubflowPsi(node, tau);
    SubflowTheta(node);
}
