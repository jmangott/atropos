#ifndef TREE_CLASS_HPP
#define TREE_CLASS_HPP

#include <iostream>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include <netcdf.h>

#include "coeff_class.hpp"
#include "grid_class.hpp"
#include "netcdf_check.hpp"
#include "matrix.hpp"
#include "timer_class.hpp"

#ifdef __OPENMP__
#pragma omp declare reduction(+: multi_array<double, 2>: omp_out += omp_in) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig))
#endif

#ifdef __OPENMP__
#pragma omp declare reduction(+: multi_array<double, 4>: omp_out += omp_in) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig))
#endif

// General classes for the hierarchical DLR approximation
// TODO: introduce a template parameter `N` for arbitrary many outgoing legs
template<class T>
struct node
{
    const std::string id;

    node* const parent;
    std::array<node*, 2> child;
    const Index n_basisfunctions;

    multi_array<T, 2> S;

    node() = default;

    node(const std::string _id, node * const _parent, std::array<node*, 2> _child, const Index _r_in, const Index _n_basisfunctions) 
    : id(_id)
    , parent(_parent)
    , child(_child)
    , n_basisfunctions(_n_basisfunctions)
    , S({_r_in, _r_in})
    {
        assert(n_basisfunctions <= _r_in);
    }

    virtual ~node() = default;

    virtual bool IsInternal() const = 0;
    virtual bool IsExternal() const = 0;
    virtual void Initialize(int ncid) = 0;

    Index RankIn() const
    {
        return S.shape()[0];
    }
};

template<class T>
struct internal_node : virtual node<T>
{
    multi_array<T, 3> Q;
    multi_array<T, 3> G;

    internal_node(const std::string _id, internal_node * const _parent, const Index _r_in, const std::array<Index, 2> _r_out, const Index _n_basisfunctions)
    : node<T>(_id, _parent, {nullptr, nullptr}, _r_in, _n_basisfunctions)
    , Q({_r_out[0], _r_out[1], _r_in})
    , G({_r_out[0], _r_out[1], _r_in})
    {}
    
    bool IsInternal() const override
    {
        return true;
    }

    bool IsExternal() const override
    {
        return false;
    }

    std::array<Index, 2> RankOut() const
    {
        return array<Index, 2>({Q.shape()[0], Q.shape()[1]});
    }

    void Initialize(int ncid) override;

    void Write(int ncid, int id_r_in, std::array<int, 2> id_r_out) const;

    multi_array<T, 2> Orthogonalize(std::function<T(T *, T *)> inner_product, const blas_ops &blas);
};

template<class T>
struct external_node : virtual node<T>
{
    multi_array<T, 2> X;

    external_node(const std::string _id, internal_node<T> * const _parent, const Index _dx, const Index _r_in, const Index _n_basisfunctions) 
    : node<T>(_id, _parent, {nullptr, nullptr}, _r_in, _n_basisfunctions)
    , X({_dx, _r_in})
    {}
    
    bool IsInternal() const override
    {
        return false;
    }

    bool IsExternal() const override
    {
        return true;
    }

    Index ProblemSize() const
    {
        return X.shape()[0];
    }

    void Initialize(int ncid) override;

    void Write(int ncid, int id_r_in, int id_dx) const;

    multi_array<T, 2> Orthogonalize(std::function<T(T *, T *)> inner_product, const blas_ops &blas);
};

struct cme_node : virtual node<double>
{
    std::array<cme_node*, 2> child;
    const grid_parms grid;
    cme_coeff coefficients;

    cme_node(const std::string _id, cme_node * const _parent, const grid_parms _grid, const Index _r_in, const Index _n_basisfunctions)
    : node<double>(_id, _parent, {nullptr, nullptr}, _r_in, _n_basisfunctions)
    , child({nullptr, nullptr})
    , grid(_grid)
    , coefficients(_grid.n_reactions, _r_in)
    {}
    void CalculateAB_bar(const blas_ops &blas);
    void CalculateEF(const blas_ops& blas);
};

struct cme_internal_node : cme_node, internal_node<double>
{
    cme_internal_coeff internal_coefficients;

    cme_internal_node(const std::string _id, cme_internal_node * const _parent, const grid_parms _grid, const Index _r_in, const std::array<Index, 2> _r_out, const Index _n_basisfunctions) 
    : node(_id, _parent, {nullptr, nullptr}, _r_in, _n_basisfunctions)
    , cme_node(_id, _parent, _grid, _r_in, _n_basisfunctions)
    , internal_node<double>(_id, _parent, _r_in, _r_out, _n_basisfunctions)
    , internal_coefficients(_r_in, _r_out)
    {}
    void Initialize(int ncid) override;

    template <Index id>
    void CalculateAB(const blas_ops &blas);
    void CalculateGH(const blas_ops &blas);
    void CalculateQ(const double tau);
};

struct cme_external_node : cme_node, external_node<double>
{
    cme_external_coeff external_coefficients;

    cme_external_node(const std::string _id, cme_internal_node * const _parent, const grid_parms _grid, const Index _r_in, const Index _n_basisfunctions)
    : node(_id, _parent, {nullptr, nullptr}, _r_in, _n_basisfunctions)
    , cme_node(_id, _parent, _grid, _r_in, _n_basisfunctions)
    , external_node<double>(_id, _parent, _grid.dx, _r_in, _n_basisfunctions)
    , external_coefficients(_grid.n_reactions)
    {}
    void Initialize(int ncid) override;
};

multi_array<double, 2> CalculateKDot(const multi_array<double, 2> &K, const cme_external_node* const node, const blas_ops &blas);

multi_array<double, 2> CalculateSDot(const multi_array<double, 2> &S, const cme_node *const node);

multi_array<double, 2> CalculateQDot(const multi_array<double, 2> &Qmat, const cme_internal_node* const node);

struct cme_lr_tree
{
    cme_internal_node * root;
    std::string partition_str;

    friend std::ostream &operator<<(std::ostream &os, cme_lr_tree const& tree)
    {
        tree.PrintHelper(os, tree.root);
        return os;
    }

    private:
        void PrintHelper(std::ostream &os, cme_node const * const node) const;
        void OrthogonalizeHelper(cme_internal_node * const node, const blas_ops &blas) const;
        void InitializeAB_barHelper(cme_node * const node, const blas_ops &blas) const;
        std::vector<double> NormalizeHelper(cme_node const * const node) const;

    public:
        void Read(const std::string fn);
        void Write(const std::string, const double t, const double tau, const double dm) const;
        void Orthogonalize(const blas_ops &blas) const;
        void InitializeAB_bar(const blas_ops &blas) const;
        double Normalize() const;
};

namespace WriteHelpers
{
    void WritePartitionStr(int ncid, const std::string partition_str);
    void WriteGridParms(int ncid, const grid_parms grid);
    void WriteNode(int ncid, cme_node const * const node);
}

namespace ReadHelpers
{
    std::string ReadPartitionStr(int ncid);
    grid_parms ReadGridParms(int ncid);
    std::array<Index, 2> ReadRankOut(int ncid);
    Index ReadNBasisfunctions(int ncid);
    std::vector<std::vector<double>> ReadPropensity(int ncid, const Index n_reactions);
    cme_node *ReadNode(int ncid, const std::string id, cme_internal_node * const parent_node, const Index r_in);
}

template <class T>
multi_array<T, 2> internal_node<T>::Orthogonalize(std::function<T(T *, T *)> inner_product, const blas_ops &blas)
{
    multi_array<T, 2> Qmat({prod(RankOut()), node<T>::RankIn()});
    multi_array<T, 2> Q_R({node<T>::RankIn(), node<T>::RankIn()});
    Matrix::Matricize(Q, Qmat, 2);
    Q_R = Matrix::Orthogonalize(Qmat, node<T>::n_basisfunctions, inner_product, blas);
    Matrix::Tensorize(Qmat, Q, 2);

    return Q_R;
};

template <class T>
multi_array<T, 2> external_node<T>::Orthogonalize(std::function<T(T *, T *)> inner_product, const blas_ops &blas)
{
    multi_array<T, 2> X_R({node<T>::RankIn(), node<T>::RankIn()});
    X_R = Matrix::Orthogonalize(X, node<T>::n_basisfunctions, inner_product, blas);

    return X_R;
};

template <Index id>
void cme_internal_node::CalculateAB(const blas_ops &blas)
{
    const Index id_c = (id == 0) ? 1 : 0;
    Index rank_out_c = RankOut()[id_c];

    std::vector<multi_array<double, 2>> AA_bar(grid.n_reactions);
    std::vector<multi_array<double, 2>> BB_bar(grid.n_reactions);

    multi_array<double, 2> Gmat({RankIn() * rank_out_c, RankOut()[id]});

    multi_array<double, 2> A_tmp1({RankIn() * rank_out_c, RankOut()[id]});
    multi_array<double, 2> B_tmp1({RankIn() * rank_out_c, RankOut()[id]});
    multi_array<double, 2> A_tmp2({RankOut()[id], RankOut()[id]});
    multi_array<double, 2> B_tmp2({RankOut()[id], RankOut()[id]});

#ifdef __OPENMP__
#pragma omp parallel for
#endif
    for (Index mu = 0; mu < child[id]->grid.n_reactions; ++mu)
    {
        set_zero(child[id]->coefficients.A[mu]);
        set_zero(child[id]->coefficients.B[mu]);
        AA_bar[mu].resize({RankIn() * rank_out_c, RankIn() * rank_out_c});
        BB_bar[mu].resize({RankIn() * rank_out_c, RankIn() * rank_out_c});
    }

#ifdef __OPENMP__
#pragma omp parallel for
#endif
    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        for (Index i1 = 0; i1 < rank_out_c; ++i1)
        {
            for (Index j1 = 0; j1 < rank_out_c; ++j1)
            {
                for (Index i = 0; i < RankIn(); ++i)
                {
                    for (Index j = 0; j < RankIn(); ++j)
                    {
                        AA_bar[mu](i1 + rank_out_c * i, j1 + rank_out_c * j) = coefficients.A[mu](i, j) * child[id_c]->coefficients.A_bar[mu](i1, j1);
                        BB_bar[mu](i1 + rank_out_c * i, j1 + rank_out_c * j) = coefficients.B[mu](i, j) * child[id_c]->coefficients.B_bar[mu](i1, j1);
                    }
                }
            }
        }
    }

    Matrix::Matricize(G, Gmat, id);

#ifdef __OPENMP__
#pragma omp parallel for
#endif
    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        blas.matmul(AA_bar[mu], Gmat, A_tmp1);
        blas.matmul_transa(Gmat, A_tmp1, A_tmp2);

        blas.matmul(BB_bar[mu], Gmat, B_tmp1);
        blas.matmul_transa(Gmat, B_tmp1, B_tmp2);

        child[id]->coefficients.A[mu] += A_tmp2;
        child[id]->coefficients.B[mu] += B_tmp2;
    }
};

#endif