#ifndef TREE_CLASS_HPP
#define TREE_CLASS_HPP

#include<memory>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "grid_class.hpp"


// struct for storing low-rank factors and coefficients
template<class T>
struct node
{
    std::unique_ptr<node<T>> left;
    std::unique_ptr<node<T>> right;
};

template<class T>
struct root : node<T>
{
    Index r;
    multi_array<double, 2> s_matrix;
};

template<class T>
struct internal_node : node<T>
{
    Index r;
    T internal_coeff;
    multi_array<double, 2> s_matrix;
    multi_array<double, 3> q_tensor;
    multi_array<double, 3> g_tensor;
};

template<class T>
struct external_node : node<T>
{
    std::unique_ptr<node<T>> left = nullptr;
    std::unique_ptr<node<T>> right = nullptr;
    grid_info grid;
    T external_coeff;
    multi_array<double, 2> s_matrix;
};

template<class T>
struct tree
{
    std::unique_ptr<root<T>> root_node;
};

#endif