#ifndef TREE_CLASS_HPP
#define TREE_CLASS_HPP

#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "coeff_class.hpp"
#include "grid_class.hpp"


// General classes for the hierarchical DLR approximation
template<class T>
struct node
{
    node<T>* left = nullptr;
    node<T>* right = nullptr;

    virtual ~node() {};
};

template<class T>
struct root : node<T>
{
    T root_coeff;
    multi_array<double, 2> S;

    root(Index _r) : S({_r, _r}) {};

    Index rank() const
    {
        return S.shape()[0];
    }
};

template<class T>
struct internal_node : node<T>
{
    node<T>* parent;
    T internal_coeff;
    multi_array<double, 3> Q;
    multi_array<double, 3> G;
    multi_array<double, 2> S;

    internal_node(root<T>* _parent, Index _r) : Q({_parent->rank(), _r, _r}), G({_parent->rank(), _r, _r}), S({_r, _r}) {};
    internal_node(internal_node<T>* _parent, Index _r) : Q({_parent->rank(), _r, _r}), G({_parent->rank(), _r, _r}), S({_r, _r}) {};

    Index rank() const
    {
        return S.shape()[0];
    }
};

template<class T>
struct external_node : node<T>
{
    node<T>* parent;
    T external_coeff;
    multi_array<double, 2> X;

    external_node(root<T>* _parent, Index _dx) : X({_parent->rank(), _dx}) {};
    external_node(internal_node<T>* _parent, Index _dx) : X({_parent->rank(), _dx}) {};

    Index problem_size() const
    {
        return X.shape()[0];
    }
};

template<class T>
struct lr_tree
{
    root<T>* root_node;

    // lr_tree(root<T>* _root_node) : root_node(_root_node) {};
};


// Dervied, CME-specific classes
template<class T>
struct cme_root : root<T>
{
    grid_parms grid;

    cme_root(grid_parms _grid, Index _r) : root<T>(_r), grid(_grid) {};
};

template<class T>
struct cme_internal_node : internal_node<T>
{
    grid_parms grid;

    cme_internal_node(root<T> *_parent, grid_parms _grid, Index _r) : internal_node<T>(_parent, _r), grid(_grid){};
    cme_internal_node(internal_node<T> *_parent, grid_parms _grid, Index _r) : internal_node<T>(_parent, _r), grid(_grid) {};
};

template<class T>
struct cme_external_node : external_node<T>
{
    grid_parms grid;

    cme_external_node(root<T> *_parent, grid_parms _grid) : external_node<T>(_parent, _grid.dx()), grid(_grid) {};
    cme_external_node(internal_node<T> *_parent, grid_parms _grid) : external_node<T>(_parent, _grid.dx()), grid(_grid) {};
};

#endif