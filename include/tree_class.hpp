#ifndef TREE_CLASS_HPP
#define TREE_CLASS_HPP

#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include <netcdf.h>

#include "coeff_class.hpp"
#include "grid_class.hpp"
#include "io_functions.hpp"
#include "netcdf_error.hpp"
#include "matrix.hpp"

// General classes for the hierarchical DLR approximation
// TODO: introduce a template parameter `N` for arbitrary many outgoing legs
struct node
{
    std::string id;

    node* parent;
    std::array<node*, 2> child;

    node(std::string _id, node* _parent, std::array<node*, 2> _child) 
    : id(_id)
    , parent(_parent)
    , child(_child) {}

    virtual ~node() {}

    virtual bool IsInternal() const = 0;
    virtual bool IsExternal() const = 0;
    virtual Index RankIn() const = 0;
    virtual void Initialize(int ncid) = 0;
};

template<class T>
struct internal_node : node
{
    multi_array<T, 3> Q;
    multi_array<T, 3> G;
    multi_array<T, 2> S;
    Index n_basisfunctions;

    internal_node(const std::string _id, internal_node * const _parent, const Index _r_in, const std::array<Index, 2> _r_out, const Index _n_basisfunctions)
    : node(_id, _parent, {nullptr, nullptr})
    , Q((_parent == nullptr) ? std::array<Index, 3>({_r_out[0], _r_out[1], 1}) : std::array<Index, 3>({_r_out[0], _r_out[1], _r_in}))
    , G((_parent == nullptr) ? std::array<Index, 3>({_r_out[0], _r_out[1], 1}) : std::array<Index, 3>({_r_out[0], _r_out[1], _r_in}))
    , S(_r_out)
    , n_basisfunctions(_n_basisfunctions)
    {
        assert(n_basisfunctions <= _r_in);
    }
    
    bool IsInternal() const
    {
        return true;
    }

    bool IsExternal() const
    {
        return false;
    }

    std::array<Index, 2> RankOut() const
    {
        return S.shape();
    }

    Index RankIn() const
    {
        return Q.shape()[2];
    }

    void Initialize(int ncid);

    multi_array<T, 2> Orthogonalize(std::function<T(T *, T *)> inner_product, const blas_ops &blas);
};

template<class T>
struct external_node : node
{
    multi_array<T, 2> X;
    Index n_basisfunctions;

    external_node(std::string _id, internal_node<T> * const _parent, const Index _dx, const Index _r_in, const Index _n_basisfunctions) 
    : node(_id, _parent, {nullptr, nullptr})
    , X({_dx, _r_in})
    , n_basisfunctions(_n_basisfunctions)
    {
        assert(n_basisfunctions < _r_in);
    }
    
    bool IsInternal() const
    {
        return false;
    }

    bool IsExternal() const
    {
        return true;
    }

    Index ProblemSize() const
    {
        return X.shape()[0];
    }

    Index RankIn() const
    {
        return X.shape()[1];
    }

    void Initialize(int ncid);

    multi_array<T, 2> Orthogonalize(std::function<T(T *, T *)> inner_product, const blas_ops &blas);
};

struct cme_internal_node : internal_node<double>
{
    grid_parms grid;
    cme_internal_coeff coefficients;

    cme_internal_node(std::string _id, cme_internal_node *_parent, grid_parms _grid, Index _r_in, std::array<Index, 2> _r_out, Index _n_basisfunctions) 
    : internal_node<double>(_id, _parent, _r_in, _r_out, _n_basisfunctions)
    , grid(_grid), coefficients(_grid.n_reactions, _r_in, _r_out) {}
};

struct cme_external_node : external_node<double>
{
    grid_parms grid;
    cme_external_coeff coefficients;

    cme_external_node(std::string _id, cme_internal_node *_parent, grid_parms _grid, Index _r_in, Index _n_basisfunctions) 
    : external_node<double>(_id, _parent, _grid.dx, _r_in, _n_basisfunctions)
    , grid(_grid), coefficients(_grid.n_reactions) {}
};

struct cme_lr_tree
{
    cme_internal_node *root;

    private:
        void PrintHelper(node* node);
        void OrthogonalizeHelper(cme_internal_node *node);

    public:
        void Read(std::string fn);
        void Print();
        void Orthogonalize();
};

namespace ReadHelpers
{
    grid_parms ReadGridParms(int ncid);

    std::array<Index, 2> ReadRankOut(int ncid);

    Index ReadNBasisfunctions(int ncid);

    node *ReadNode(int ncid, std::string id, cme_internal_node *parent_node, Index r_in);
}

template <class T>
multi_array<T, 2> internal_node<T>::Orthogonalize(std::function<T(T *, T *)> inner_product, const blas_ops &blas)
{
    multi_array<T, 2> Qmat({prod(RankOut()), RankIn()});
    multi_array<T, 2> Q_R({RankIn(), RankIn()});
    Matrix::Matricize(Q, Qmat, 2);
    Q_R = Matrix::Orthogonalize(Qmat, n_basisfunctions, inner_product, blas);
    Matrix::Tensorize(Qmat, Q, 2);

    return Q_R;
};

template <class T>
multi_array<T, 2> external_node<T>::Orthogonalize(std::function<T(T *, T *)> inner_product, const blas_ops &blas)
{
    multi_array<T, 2> X_R({RankIn(), RankIn()});
    X_R = Matrix::Orthogonalize(X, n_basisfunctions, inner_product, blas);

    return X_R;
};

#endif