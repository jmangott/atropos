#ifndef TREE_CLASS_HPP
#define TREE_CLASS_HPP

#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "coeff_class.hpp"
#include "grid_class.hpp"


// struct for storing low-rank factors and coefficients
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
    grid_parms grid;
    T root_coeff;
    multi_array<double, 2> S;

    root(grid_parms _grid, Index _r) : grid(_grid), S({_r, _r}) {};

    Index rank() const
    {
        return S.shape()[0];
    }
};

template<class T>
struct internal_node : node<T>
{
    node<T>* parent;
    grid_parms grid;
    T internal_coeff;
    multi_array<double, 3> Q;
    multi_array<double, 3> G;
    multi_array<double, 2> S;

    internal_node(root<T>* _parent, grid_parms _grid, Index _r) : grid(_grid), Q({_parent->rank(), _r, _r}), G({_parent->rank(), _r, _r}), S({_r, _r}) {};
    internal_node(internal_node<T>* _parent, grid_parms _grid, Index _r) : grid(_grid), Q({_parent->rank(), _r, _r}), G({_parent->rank(), _r, _r}), S({_r, _r}) {};

    Index rank() const
    {
        return S.shape()[0];
    }
};

template<class T>
struct external_node : node<T>
{
    node<T>* parent;
    grid_parms grid;
    T external_coeff;
    multi_array<double, 2> X;

    external_node(root<T>* _parent, grid_parms _grid) : grid(_grid), X({_parent->rank(), _grid.dx}) {};
    external_node(internal_node<T>* _parent, grid_parms _grid) : grid(_grid), X({_parent->rank(), _grid.dx}) {};

    Index problem_size() const
    {
        return X.shape()[0];
    }
};

template<class T>
struct tree
{
    root<T>* root_node;
};

#endif