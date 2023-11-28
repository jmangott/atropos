"""Contains the `InitialCondition` class for setting up initial conditions."""
from tree_class import *

# TODO: Q should have shape (n_basisfunctions, child(n_basisfunctions), child(n_basisfunctions))

class InitialCondition:
    """
    Provides all S matrices, the Q tensors and the low-rank factors as an array according to the following (recursive) ordering convention (OC):

        1. left node (or root node)
        2. apply OC
        3. right node

    These arrays can be conveniently set to obey the initial conditions. For setting up the low-rank factors, the grid objects of the external nodes are also provided as an array (`grid_external_nodes`).
    """

    def __setNodeData(self, node: Node, n_basisfunctions_iter):
        nb = next(n_basisfunctions_iter)

        if isinstance(node.child[0], ExternalNode) and isinstance(node.child[1], ExternalNode):
            node.child[0].X.resize((nb, node.child[0].grid.dx), refcheck=False)
            node.child[1].X.resize((nb, node.child[1].grid.dx), refcheck=False)

        elif isinstance(node.child[0], ExternalNode) and isinstance(node.child[1], InternalNode):
            node.child[0].X.resize((nb, node.child[0].grid.dx), refcheck=False)
            node.child[1].Q.resize((nb, node.child[1].r_out, node.child[1].r_out), refcheck=False)

            self.__setNodeData(node.child[1], n_basisfunctions_iter)

        elif isinstance(node.child[0], InternalNode) and isinstance(node.child[1], ExternalNode):
            node.child[0].Q.resize((nb, node.child[0].r_out, node.child[0].r_out), refcheck=False)
            node.child[1].X.resize((nb, node.child[1].grid.dx), refcheck=False)

            self.__setNodeData(node.child[0], n_basisfunctions_iter)

        else:
            node.child[0].Q.resize((nb, node.child[0].r_out, node.child[0].r_out), refcheck=False)
            node.child[1].Q.resize((nb, node.child[1].r_out, node.child[1].r_out), refcheck=False)

            self.__setNodeData(node.child[0], n_basisfunctions_iter)
            self.__setNodeData(node.child[1], n_basisfunctions_iter)

    def __getNodeData(self, node: Node):
        if isinstance(node, ExternalNode):
            self.external_nodes.append(node)
            self.X.append(node.X)
        else:
            self.Q.append(node.Q)
            self.S.append(node.S)
            self.__getNodeData(node.child[0])
            self.__getNodeData(node.child[1])

    def __init__(self, _tree: Tree, _n_basisfunctions: npt.NDArray[np.int_]):
        if (_n_basisfunctions.size != _tree.n_internal_nodes):
            raise ValueError(
                "`_n_basisfunctions.size` must be equal to the number of internal nodes")
        
        if (np.any(_n_basisfunctions > _tree.r_out)):
            raise ValueError(
                "`_n_basisfunctions` must be smaller or equal than the incoming rank")

        self.n_basisfunctions = _n_basisfunctions
        self.external_nodes = []
        self.S = []
        self.Q = []
        self.X = []

        n_basisfunctions_iter = iter(self.n_basisfunctions)
        self.__setNodeData(_tree.root, n_basisfunctions_iter)
        self.__getNodeData(_tree.root)

if __name__ == "__main__":
    import models.lambda_phage as model

    # Partition string
    partition_str = "((0 1)(2 3))(4)"

    # Rank
    r_out = np.array([5, 4])

    # Grid parameters
    n = np.array([16, 41, 11, 11, 11])
    d = n.size
    binsize = np.ones(d, dtype=int)
    liml = np.zeros(d)
    grid = GridParms(n, binsize, liml)

    # Set up the partition tree
    tree = Tree(model.reaction_system, partition_str, grid, r_out)
    tree.buildTree()

    # Initial distribution
    def gaussian(x: np.ndarray, mu: np.ndarray, C) -> float:
        Cinv = 1 / C
        return np.exp(-0.5 * np.dot(np.transpose(x - mu), np.dot(Cinv, (x - mu))))

    # Number of basisfunctions
    n_basisfunctions = np.array([2, 1])

    # Low-rank initial conditions
    initial_conditions = InitialCondition(tree, n_basisfunctions)
    for S_el in initial_conditions.S:
        S_el[0, 0] = 1

    for Q_el in initial_conditions.Q:
        Q_el[0, 0, 0] = 1

    for k, node in enumerate(initial_conditions.external_nodes):
        for i in range(node.X.shape[0]):
            X_function = (lambda x: gaussian(x, np.zeros((node.grid.d)), 0.5))
            vec_index = np.zeros(node.grid.d)
            for j in range(node.X.shape[1]):
                state_x = vecIndexToState(vec_index, node.grid.liml, node.grid.binsize)
                # TODO: Check whether state_x + grid.binsize * 0.5 is correct
                initial_conditions.X[k][i, j] = X_function(
                    state_x + node.grid.binsize * 0.5)
                incrVecIndex(vec_index, node.grid.n, node.grid.d)
            # Normalization
            initial_conditions.X[k][i] /= np.linalg.norm(initial_conditions.X[k][i])

    # Print tree and write it to a netCDF file
    tree.printTree()
    tree.writeTree()