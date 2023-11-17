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
struct node
{
    std::string id;

    node* parent;
    node* left;
    node* right;

    node(std::string _id, node* _parent, node* _left, node* _right) : id(_id), parent(_parent), left(_left), right(_right) {};

    virtual ~node() {};

    // virtual bool IsRoot() = 0;
    virtual bool IsInternal() const = 0;
    virtual bool IsExternal() const = 0;
    virtual Index ParentRank() const = 0;
    virtual void Initialize(int ncid) = 0;
};


// template<class T>
// struct root_node : node
// {
//     multi_array<T, 2> S;

//     root_node(Index _r) : node("root", nullptr, nullptr, nullptr), S({_r, _r}) {};

//     bool IsRoot()
//     {
//         return true;
//     }
    
//     bool IsInternal()
//     {
//         return false;
//     }
    
//     bool IsExternal()
//     {
//         return false;
//     }

//     Index rank() const
//     {
//         return S.shape()[0];
//     }

//     void InitializeNode(int ncid);
// };

template<class T>
struct internal_node : node
{
    multi_array<T, 3> Q;
    multi_array<T, 3> G;
    multi_array<T, 2> S;
    Index n_basisfunctions;

    // internal_node(std::string _id, root_node<T>* _parent, Index _r, Index _n_basisfunctions) : node(_id, _parent, nullptr, nullptr), Q({_r, _r, _parent->rank()}), G({_r, _r, _parent->rank()}), S({_r, _r}), n_basisfunctions(_n_basisfunctions) {};
    internal_node(std::string _id, internal_node *_parent, Index _r, Index _n_basisfunctions)
        : node(_id, _parent, nullptr, nullptr), Q((_parent == nullptr) ? std::array<Index, 3>({_r, _r, 1}) : std::array<Index, 3>({_r, _r, _parent->Rank()})), G((_parent == nullptr) ? std::array<Index, 3>({_r, _r, 1}) : std::array<Index, 3>({_r, _r, _parent->Rank()})), S({_r, _r}), n_basisfunctions(_n_basisfunctions){};

    // bool IsRoot()
    // {
    //     return false;
    // }
    
    bool IsInternal() const
    {
        return true;
    }

    bool IsExternal() const
    {
        return false;
    }

    Index Rank() const
    {
        return S.shape()[0];
    }

    Index ParentRank() const
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

    // external_node(std::string _id, root_node<T> *_parent, Index _dx, Index _n_basisfunctions) : node(_id, _parent, nullptr, nullptr), X({_dx, _parent->rank()}), n_basisfunctions(_n_basisfunctions){};
    external_node(std::string _id, internal_node<T> *_parent, Index _dx, Index _n_basisfunctions) : node(_id, _parent, nullptr, nullptr), X({_dx, _parent->Rank()}), n_basisfunctions(_n_basisfunctions){};

    // bool IsRoot()
    // {
    //     return false;
    // }
    
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

    Index ParentRank() const
    {
        return X.shape()[1];
    }

    void Initialize(int ncid);

    multi_array<T, 2> Orthogonalize(std::function<T(T *, T *)> inner_product, const blas_ops &blas);
};

// Derived, CME-specific classes
// struct cme_root_node : root_node<double>
// {
//     grid_parms grid;
//     cme_root_coeff coefficients;

//     cme_root_node(grid_parms _grid, Index _r) : root_node<double>(_r), grid(_grid) {};
// };

struct cme_internal_node : internal_node<double>
{
    grid_parms grid;
    cme_internal_coeff coefficients;

    // cme_internal_node(std::string _id, cme_root_node *_parent, grid_parms _grid, Index _r, Index _n_basisfunctions) : internal_node<double> (_id, _parent, _r, _n_basisfunctions), grid(_grid) {};
    cme_internal_node(std::string _id, cme_internal_node *_parent, grid_parms _grid, Index _r, Index _n_basisfunctions) : internal_node<double>(_id, _parent, _r, _n_basisfunctions), grid(_grid) {};

    // multi_array<double, 2> Orthogonalize(std::function<double(double *, double *)> inner_product, const blas_ops &blas);
};

struct cme_external_node : external_node<double>
{
    grid_parms grid;
    cme_external_coeff coefficients;

    // cme_external_node(std::string _id, cme_root_node* _parent, grid_parms _grid, Index _n_basisfunctions) : external_node<double>(_id, _parent, _grid.dx(), _n_basisfunctions), grid(_grid) {};
    cme_external_node(std::string _id, cme_internal_node* _parent, grid_parms _grid, Index _n_basisfunctions) : external_node<double>(_id, _parent, _grid.dx(), _n_basisfunctions), grid(_grid) {};

    // multi_array<double, 2> Orthogonalize(std::function<double(double *, double *)> inner_product, const blas_ops &blas);
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

    Index ReadRank(int ncid);

    Index ReadNBasisfunctions(int ncid);

    node *ReadNode(int ncid, std::string id, cme_internal_node *parent_node);
}

template <class T>
multi_array<T, 2> internal_node<T>::Orthogonalize(std::function<T(T *, T *)> inner_product, const blas_ops &blas)
{
    multi_array<T, 2> Qmat({Rank() * Rank(), ParentRank()});
    multi_array<T, 2> Q_R({ParentRank(), ParentRank()});
    Matrix::Matricize(Q, Qmat, 2);
    Q_R = Matrix::Orthogonalize(Qmat, n_basisfunctions, inner_product, blas);
    Matrix::Tensorize(Qmat, Q, 0);

    return Q_R;
};

template <class T>
multi_array<T, 2> external_node<T>::Orthogonalize(std::function<T(T *, T *)> inner_product, const blas_ops &blas)
{
    multi_array<T, 2> X_R({ParentRank(), ParentRank()});
    X_R = Matrix::Orthogonalize(X, n_basisfunctions, inner_product, blas);

    return X_R;
};

#endif