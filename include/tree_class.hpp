#ifndef TREE_CLASS_HPP
#define TREE_CLASS_HPP

#include <variant>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include <netcdf.h>

#include "coeff_class.hpp"
#include "grid_class.hpp"
#include "netcdf_error.hpp"

// General classes for the hierarchical DLR approximation
struct node
{
    std::string id;

    node* parent;
    node* left;
    node* right;

    node(std::string _id, node* _parent, node* _left, node* _right) : id(_id), parent(_parent), left(_left), right(_right) {};

    virtual ~node() {};

    virtual bool IsRoot() = 0;
    virtual bool IsExternal() = 0;
    virtual void InitializeNode(int ncid) = 0;
};

struct root_node : node
{
    multi_array<double, 2> S;

    root_node(Index _r) : node("root", nullptr, nullptr, nullptr), S({_r, _r}) {};

    bool IsRoot()
    {
        return true;
    }
    
    bool IsExternal()
    {
        return false;
    }

    void InitializeNode(int ncid);

    Index rank() const
    {
        return S.shape()[0];
    }
};

struct internal_node : node
{
    multi_array<double, 3> Q;
    multi_array<double, 3> G;
    multi_array<double, 2> S;

    internal_node(std::string _id, root_node* _parent, Index _r) : node(_id, _parent, nullptr, nullptr), Q({_parent->rank(), _r, _r}), G({_parent->rank(), _r, _r}), S({_r, _r}) {};
    internal_node(std::string _id, internal_node *_parent, Index _r) : node(_id, _parent, nullptr, nullptr), Q({_parent->rank(), _r, _r}), G({_parent->rank(), _r, _r}), S({_r, _r}) {};

    bool IsRoot()
    {
        return false;
    }

    bool IsExternal()
    {
        return false;
    }

    void InitializeNode(int ncid);

    Index rank() const
    {
        return S.shape()[0];
    }
};

struct external_node : node
{
    multi_array<double, 2> X;

    external_node(std::string _id, root_node* _parent, Index _dx) : node(_id, _parent, nullptr, nullptr), X({_parent->rank(), _dx}) {};
    external_node(std::string _id, internal_node *_parent, Index _dx) : node(_id, _parent, nullptr, nullptr), X({_parent->rank(), _dx}) {};

    bool IsRoot()
    {
        return false;
    }

    bool IsExternal()
    {
        return true;
    }

    void InitializeNode(int ncid);

    Index problem_size() const
    {
        return X.shape()[0];
    }

// TODO: Ask Lukas how orthogonalization should be done, inplace or with a temporary array?
    Index n_basisfunctions() const
    {
        return X.shape()[1];
    }
};

// Derived, CME-specific classes
struct cme_root_node : root_node
{
    grid_parms grid;
    root_coeff coefficients;

    cme_root_node(grid_parms _grid, Index _r) : root_node(_r), grid(_grid) {};
};

struct cme_internal_node : internal_node
{
    grid_parms grid;
    internal_coeff coefficients;

    cme_internal_node(std::string _id, cme_root_node *_parent, grid_parms _grid, Index _r) : internal_node
    (_id, _parent, _r), grid(_grid) {};
    cme_internal_node(std::string _id, cme_internal_node *_parent, grid_parms _grid, Index _r) : internal_node(_id, _parent, _r), grid(_grid) {};
};

struct cme_external_node : external_node
{
    grid_parms grid;
    external_coeff coefficients;

    cme_external_node(std::string _id, cme_root_node* _parent, grid_parms _grid) : external_node(_id, _parent, _grid.dx()), grid(_grid) {};
    cme_external_node(std::string _id, cme_internal_node* _parent, grid_parms _grid) : external_node(_id, _parent, _grid.dx()), grid(_grid) {};
};

struct cme_lr_tree
{
    cme_root_node *root;

    private:
        struct ReadTreeHelpers;
        ReadTreeHelpers* pReadTreeHelpers;
        void PrintTreeHelper(node* node);

    public:
        void ReadTree(std::string fn);
        void PrintTree();
};

#endif