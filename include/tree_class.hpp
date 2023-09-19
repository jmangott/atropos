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
    node<T>* left;
    node<T>* right;
};

template<class T>
struct root : node<T>
{
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

    // internal_node(node<T>* _parent, Index _r) : Q({_parent->rank(), _r, _r}), G({_parent->rank(), _r, _r}), S({_r, _r}) {};
    internal_node(node<T>* _parent, Index _r) : Q({_r, _r, _r}), G({_r, _r, _r}), S({_r, _r}) {};

    Index rank() const
    {
        return S.shape()[0];
    }    
};

template<class T>
struct external_node : node<T>
{
    node<T>* parent;
    node<T>* left = nullptr;
    node<T>* right = nullptr;
    grid_parms grid;
    T external_coeff;
    multi_array<double, 2> X;

    // external_node<T>(node<T>* _parent, grid_parms _grid) : grid(_grid), X({_parent->rank(), _grid.dx}) {};
    external_node<T>(node<T>* _parent, grid_parms _grid) : grid(_grid), X({1, _grid.dx}) {};

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