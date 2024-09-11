#include "ps_integrator.hpp"

template <Index id>
void ps_integrator::SubflowPhi(cme_internal_node * const node, const double tau) const
{
    Index id_c = (id == 0) ? 1 : 0;

    orthogonalize gs(&blas);
    multi_array<double, 2> Qmat({node->RankIn() * node->RankOut()[id_c], node->RankOut()[id]});
    // std::function<double(double *, double *)> ip;

    // Compute QR decomposition C^n = G^n * (S^(n+id))^T
    get_time::start("Mat/Ten");
    Matrix::Matricize<id>(node->Q, Qmat);
    get_time::stop("Mat/Ten");
    get_time::start("gs");
    gs(Qmat, node->child[id]->S, 1.0);
    get_time::stop("gs");
    get_time::start("Mat/Ten");
    Matrix::Tensorize<id>(Qmat, node->G);
    get_time::stop("Mat/Ten");
    transpose_inplace(node->child[id]->S);

    get_time::start("CalculateAB");
    node->CalculateAB<id>(blas);
    get_time::stop("CalculateAB");

    if (node->child[id]->IsExternal())
    {
        get_time::start("External");
        cme_external_node *child_node = (cme_external_node *)node->child[id];

        // Compute K = X * S
        multi_array<double, 2> tmp_x(child_node->X);
        blas.matmul(tmp_x, child_node->S, child_node->X);

        // K step
        const auto K_step_rhs = [child_node, this](const multi_array<double, 2> &K) { return CalculateKDot(K, child_node, this->blas); };
        get_time::start("Integrate K");
        integration_methods.at("K")->integrate(child_node->X, K_step_rhs, tau);
        get_time::stop("Integrate K");

        // Perform the QR decomposition K = X * S
        // std::function<double(double *, double *)> ip_x;
        // ip_x = inner_product_from_const_weight(child_node->grid.h_mult, child_node->grid.dx);
        get_time::start("gs");
        gs(child_node->X, child_node->S, child_node->grid.h_mult);
        get_time::stop("gs");
        get_time::stop("External");
    }
    else
    {
        get_time::start("Internal");
        cme_internal_node *child_node = (cme_internal_node *)node->child[id];

        // Set C^(n+i) = Q^(n+id) * S^(n+id)
        multi_array<double, 2> Cmat_child({prod(child_node->RankOut()), child_node->RankIn()});
        multi_array<double, 2> Qmat_child({prod(child_node->RankOut()), child_node->RankIn()});
        get_time::start("Mat/Ten");
        Matrix::Matricize<2>(child_node->Q, Qmat_child);
        get_time::stop("Mat/Ten");
        set_zero(Cmat_child);
        blas.matmul(Qmat_child, child_node->S, Cmat_child);
        get_time::start("Mat/Ten");
        Matrix::Tensorize<2>(Cmat_child, child_node->Q);
        get_time::stop("Mat/Ten");
        get_time::stop("Internal");

        ps_integrator::operator()(child_node, tau);

        get_time::start("Internal");
        // Compute QR decomposition C^(n+id) = Q^(n+id) * S^(n+id)
        // std::function<double(double *, double *)> ip_child;
        // ip_child = inner_product_from_const_weight(1.0, prod(child_node->RankOut()));
        Matrix::Matricize<2>(child_node->Q, Cmat_child);
        get_time::start("gs");
        gs(Cmat_child, child_node->S, 1.0);
        get_time::stop("gs");
        Matrix::Tensorize<2>(Cmat_child, child_node->Q);
        get_time::stop("Internal");
    }
    get_time::start("CalculateAB_bar");
    node->child[id]->CalculateAB_bar(blas);
    get_time::stop("CalculateAB_bar");

    // Integrate S
    get_time::start("S");
    const auto S_step_rhs = [node, this](const multi_array<double, 2> &S) { return CalculateSDot(S, node->child[id], this->blas); };
    get_time::start("Integrate S");
    integration_methods.at("S")->integrate(node->child[id]->S, S_step_rhs, -1.0 * tau);
    get_time::stop("Integrate S");

    // Set C^n = G^n * (S^(n+id))^T
    multi_array<double, 2> Gmat({node->RankIn() * node->RankOut()[id_c], node->RankOut()[id]});
    Matrix::Matricize<id>(node->G, Gmat);
    set_zero(Qmat);
    blas.matmul_transb(Gmat, node->child[id]->S, Qmat);
    Matrix::Tensorize<id>(Qmat, node->Q);
    get_time::stop("S");
}

template void ps_integrator::SubflowPhi<0>(cme_internal_node * const node, const double tau) const;

template void ps_integrator::SubflowPhi<1>(cme_internal_node * const node, const double tau) const;

void ps_integrator::SubflowPsi(cme_internal_node * const node, const double tau) const
{
    multi_array<double, 2> Qmat({prod(node->RankOut()), node->RankIn()});

    get_time::start("Mat/Ten");
    Matrix::Matricize<2>(node->Q, Qmat);
    get_time::stop("Mat/Ten");

    const auto Q_step_rhs = [node, this](const multi_array<double, 2> &Qmat) { return CalculateQDot(Qmat, node, this->blas); };
    get_time::start("Integrate Q");
    integration_methods.at("Q")->integrate(Qmat, Q_step_rhs, tau);
    get_time::stop("Integrate Q");

    get_time::start("Mat/Ten");
    Matrix::Tensorize<2>(Qmat, node->Q);
    get_time::stop("Mat/Ten");
}

void ps_integrator::operator()(cme_internal_node * const node, const double tau) const
{
    SubflowPhi<0>(node, tau);
    SubflowPhi<1>(node, tau);
    SubflowPsi(node, tau);
}