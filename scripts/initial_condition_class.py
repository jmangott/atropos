"""Contains the `InitialCondition` class for setting up initial conditions."""

import numpy as np
import numpy.typing as npt

from scripts.tree_class import Tree, Node, ExternalNode, InternalNode

# TODO: Q should have shape (n_basisfunctions, child(n_basisfunctions), child(n_basisfunctions))


class InitialCondition:
    """
    Provides the Q tensors and the low-rank factors as an array according to the following (recursive) ordering convention (OC):

        1. left node (or root node)
        2. apply OC
        3. right node

    These arrays can be conveniently set to obey the initial conditions. For setting up the low-rank factors, the external nodes are also accessible via the array `external_nodes`.
    """

    def __setNodeData(self, node: Node, nb):
        if nb > node.rankIn():
            raise Exception(
                "Number of basisfunctions must be smaller or equal to the incoming rank"
            )

        if isinstance(node, ExternalNode):
            node.X.resize((node.grid.dx(), nb), refcheck=False)

        elif isinstance(node, InternalNode):
            node.Q.resize((node.rankOut()[0], node.rankOut()[1], nb), refcheck=False)

            if isinstance(self.n_basisfunctions, dict):
                current_id = node.id
                nb0 = self.n_basisfunctions[current_id + 0]
                nb1 = self.n_basisfunctions[current_id + 1]
                self.__setNodeData(node.child[0], nb0)
                self.__setNodeData(node.child[1], nb1)
            elif isinstance(self.n_basisfunctions, np.ndarray):
                next_nb = next(self.n_basisfunctions_iter)
                self.__setNodeData(node.child[0], next_nb)
                self.__setNodeData(node.child[1], next_nb)

    def __getNodeData(self, node: Node):
        if isinstance(node, ExternalNode):
            self.external_nodes.append(node)
            self.X.append(node.X)
        elif isinstance(node, InternalNode):
            self.Q.append(node.Q)
            self.__getNodeData(node.child[0])
            self.__getNodeData(node.child[1])

    def __init__(self, _tree: Tree, _n_basisfunctions: npt.NDArray[np.int_]):
        if isinstance(_n_basisfunctions, np.ndarray):
            if _n_basisfunctions.size != _tree.n_internal_nodes:
                raise Exception(
                    "`_n_basisfunctions.size` must be equal to the number of internal nodes"
                )

        if isinstance(_n_basisfunctions, dict):
            self.n_basisfunctions = _n_basisfunctions
        elif isinstance(_n_basisfunctions, np.ndarray):
            self.n_basisfunctions = _n_basisfunctions
            self.n_basisfunctions_iter = iter(_n_basisfunctions)
        else:
            raise Exception(
                "`_n_basisfunctions` must be either a `dict` or a `np.ndarray`"
            )
        self.external_nodes = []
        self.Q = []
        self.X = []

        if isinstance(self.n_basisfunctions, dict):
            current_id = _tree.root.id
            nb0 = self.n_basisfunctions[current_id + 0]
            nb1 = self.n_basisfunctions[current_id + 1]
            self.__setNodeData(_tree.root.child[0], nb0)
            self.__setNodeData(_tree.root.child[1], nb1)
        elif isinstance(self.n_basisfunctions, np.ndarray):
            nb = next(self.n_basisfunctions_iter)
            self.__setNodeData(_tree.root.child[0], nb)
            self.__setNodeData(_tree.root.child[1], nb)
        self.__getNodeData(_tree.root)
