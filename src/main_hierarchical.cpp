#include "tree_class.hpp"
#include "coeff_class.hpp"

int main()
{
    Index d = 3;
    multi_array<Index, 1> n({d});
    n(0) = 1;
    n(1) = 2;
    n(2) = 3;
    multi_array<Index, 1> binsize({d});
    binsize(0) = 1;
    binsize(1) = 2;
    binsize(2) = 3;
    multi_array<double, 1> liml({d});
    liml(0) = 0.0;
    liml(1) = 0.0;
    liml(2) = 0.0;
    grid_parms grid(d, n, binsize, liml);
    root<coeff> root_node(1);
    internal_node<coeff> child_node(&root_node, 2);
    external_node<coeff> child_child_node(&child_node, grid);
    root_node.left = &child_node;
    cout << root_node.rank() << " " << child_node.rank() << endl;
    root_node.left->left = &child_child_node;
    // cout << root_node.left->left->grid.d << endl;

    // cout << child_node.parent->rank() << endl;
    // cout << root_node.left->rank() << endl;
    return 0;
}