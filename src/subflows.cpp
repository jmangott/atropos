#include "subflows.hpp"

template <Index id>
void SubflowPhi(cme_internal_node *node, const blas_ops &blas, const double tau)
{
    Index id_c = (id == 0) ? 1 : 0;

    gram_schmidt gs(&blas);
    multi_array<double, 2> Qmat({node->RankIn() * node->RankOut()[id_c], node->RankOut()[id]});
    std::function<double(double *, double *)> ip;
    ip = inner_product_from_const_weight(1.0, node->RankIn() * node->RankOut()[id_c]);

    // Compute QR decomposition C^n = (S^(n+id))^T * G^n
    Matrix::Matricize(node->Q, Qmat, id);
    gs(Qmat, node->child[id]->S, ip);
    Matrix::Tensorize(Qmat, node->G, id);
    transpose_inplace(node->child[id]->S);

    node->CalculateAB<id>(blas);
    if (node->child[id]->IsExternal())
    {
        cme_external_node *child_node = (cme_external_node *)node->child[id];

        // Compute coefficients C and D
        child_node->CalculateCD(blas);

        // Compute K = X * S
        multi_array<double, 2> tmp_x(child_node->X);
        blas.matmul(tmp_x, child_node->S, child_node->X);

        // K step
        // TODO:`SubflowK` and `SubflowS` should be functions depending on the current node, not on the child node
        SubflowK(child_node, blas, tau);

        // Perform the QR decomposition K = X * S
        std::function<double(double *, double *)> ip_x;
        ip_x = inner_product_from_const_weight(child_node->grid.h_mult, child_node->grid.dx);
        gs(child_node->X, child_node->S, ip_x);
    }
    else
    {
        cme_internal_node *child_node = (cme_internal_node *)node->child[id];

        // Set C^(n+i) = Q^(n+id) * S^(n+id)
        multi_array<double, 2> Cmat_child({prod(child_node->RankOut()), child_node->RankIn()});
        multi_array<double, 2> Qmat_child({prod(child_node->RankOut()), child_node->RankIn()});
        Matrix::Matricize(child_node->Q, Qmat_child, 2);
        set_zero(Cmat_child);
        blas.matmul(Qmat_child, child_node->S, Cmat_child);
        Matrix::Tensorize(Cmat_child, child_node->Q, 2);

        TTNIntegrator(child_node, blas);

        // Compute QR decomposition C^(n+id) = Q^(n+id) * S^(n+id)
        std::function<double(double *, double *)> ip_child;
        ip = inner_product_from_const_weight(1.0, prod(child_node->RankOut()));
        Matrix::Matricize(child_node->Q, Cmat_child, 2);
        gs(Cmat_child, child_node->S, ip_child);
    }

    // Integrate S
    node->CalculateEF<id>(blas);
    SubflowS(node->child[id], blas, tau);

    // Set C^n = (S^(n+id))^T * G^n
    multi_array<double, 2> Gmat({node->RankIn() * node->RankOut()[id_c], node->RankOut()[id]});
    Matrix::Matricize(node->G, Gmat, id);
    set_zero(Qmat);
    blas.matmul_transb(Gmat, node->child[id]->S, Qmat);
    Matrix::Tensorize(Qmat, node->Q, id);
}

template void SubflowPhi<0>(cme_internal_node * const node, const blas_ops &blas, const double tau);

template void SubflowPhi<1>(cme_internal_node * const node, const blas_ops &blas, const double tau);

// TODO:
void SubflowPsi(cme_internal_node * const node, const blas_ops &blas, const double tau)
{
    // Compute coefficients g and h
    // Integrate C
}

void SubflowK(cme_external_node* const node, const blas_ops &blas, const double tau)
{
    multi_array<double, 2> prod_KC(node->X.shape());
    multi_array<double, 2> prod_KC_shift(node->X.shape());
    multi_array<double, 2> prod_KD(node->X.shape());
    multi_array<double, 2> K_dot(node->X.shape());
    set_zero(K_dot);

    std::vector<Index> vec_index(node->grid.d);

    for (Index mu = 0; mu < node->grid.n_reactions; ++mu)
    {
        set_zero(prod_KC);
        set_zero(prod_KD);
        std::fill(vec_index.begin(), vec_index.end(), 0);

#ifdef __OPENMP__
#pragma omp parallel firstprivate(vec_index)
#endif
        {
            Index alpha;
#ifdef __OPENMP__
            Index chunk_size = SetVecIndex(std::begin(vec_index), std::end(vec_index), std::begin(node->grid.n), node->grid.dx());
#endif

#ifdef __OPENMP__
#pragma omp for schedule(static, chunk_size)
#endif
            for (Index i = 0; i < node->grid.dx; ++i)
            {
                alpha = IndexFunction::VecIndexToDepCombIndex(std::begin(vec_index), std::begin(node->grid.n_dep[mu]), std::begin(node->grid.idx_dep[mu]), std::end(node->grid.idx_dep[mu]));
                IndexFunction::IncrVecIndex(std::begin(node->grid.n), std::begin(vec_index), std::end(vec_index));

                // Calculate matrix-vector multiplication of C2*K and D2*K
                for (Index j = 0; j < node->RankIn(); ++j)
                {
                    for (Index l = 0; l < node->RankIn(); ++l)
                    {
                        prod_KC(i, j) += tau * node->X(i, l) * node->external_coefficients.C[mu](alpha, j, l);
                        prod_KD(i, j) += tau * node->X(i, l) * node->external_coefficients.D[mu](alpha, j, l);
                    }
                }
            }
        }
        // Shift prod_KC
        Matrix::ShiftRows<1>(prod_KC_shift, prod_KC, node->grid, mu);

        // Calculate k_dot = shift(C * K) - D * K
        K_dot += prod_KC_shift;
        K_dot -= prod_KD;
    }
    node->X += K_dot;
}

void SubflowS(cme_node* node, const blas_ops &blas, const double tau)
{
    multi_array<double, 2> S_dot(node->S.shape());
    set_zero(S_dot);

#ifdef __OPENMP__
#pragma omp parallel for collapse(2)
#endif
    for (Index i = 0; i < node->RankIn(); i++)
    {
        for (Index j = 0; j < node->RankIn(); j++)
        {
            for (Index k = 0; k < node->RankIn(); k++)
            {
                for (Index l = 0; l < node->RankIn(); l++)
                {
                    S_dot(i, j) += tau * node->S(k, l) * node->coefficients.E(i, j, k, l);
                    S_dot(i, j) -= tau * node->S(k, l) * node->coefficients.F(i, j, k, l);
                }
            }
        }
    }
    node->S += S_dot;
}