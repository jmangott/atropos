import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial

from scripts.hierarchical.grid_class import GridParms
from scripts.hierarchical.initial_condition_class import InitialCondition
from scripts.hierarchical.tree_class import Tree
from scripts.index_functions import incrVecIndex

import scripts.hierarchical.models.bax as model

def tensorUnfold(tensor, mode):
    """
    Cf. https://stackoverflow.com/questions/49970141/using-numpy-reshape-to-perform-3rd-rank-tensor-unfold-operation
    """
    return np.reshape(np.moveaxis(tensor, mode, 0), (tensor.shape[mode], -1), order="F")

# Partition string
# partition_str = "(0 1 2 3 4)(5 6 7 8 9 10)"
# partition_str = "((0 1)(2 3 4))((5 6 7 8)(9 10))"
partition_str = "(0 1 2)(((3 4 6 7)(5 8))(9 10))"

# Rank
# r_out = np.array([5])
r_out = np.array([5, 3, 3])

# Grid parameters
n = np.array([46, 16, 16, 11, 11, 11, 4, 4, 4, 56, 56])
d = n.size
binsize = np.ones(d, dtype=int)
liml = np.zeros(d)
grid = GridParms(n, binsize, liml)

# Set up the partition tree
tree = Tree(model.reaction_system, partition_str, grid, r_out)
tree.buildTree()

C = 0.2
Cinv = 1 / C
mu = np.array([40, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0])

def eval_x(x: np.ndarray, mu: np.ndarray):
    return np.exp(-0.5 * Cinv * np.dot(np.transpose(x - mu), (x - mu)))

# Number of basisfunctions
# n_basisfunctions = np.array([1])
n_basisfunctions = np.array([1, 1, 1])

# Low-rank initial conditions
initial_conditions = InitialCondition(tree, n_basisfunctions)

for Q in initial_conditions.Q:
    Q[0, 0, 0] = 1.0

idx = 0
for node in range(tree.n_external_nodes):
    vec_index = np.zeros(initial_conditions.external_nodes[node].grid.d())
    for i in range(initial_conditions.external_nodes[node].grid.dx()):
        initial_conditions.X[node][i, :] = eval_x(vec_index, mu[idx : idx+len(vec_index)])
        incrVecIndex(vec_index, initial_conditions.external_nodes[node].grid.n, initial_conditions.external_nodes[node].grid.d())
    idx += len(vec_index)

# x0_sum = np.sum(initial_conditions.X[0], axis=0)
# x1_sum = np.sum(initial_conditions.X[1], axis=0)

x0_sum = np.sum(initial_conditions.X[0], axis=0)
# x01_sum = np.sum(initial_conditions.X[1], axis=0)
x100_sum = np.sum(initial_conditions.X[1], axis=0)
x101_sum = np.sum(initial_conditions.X[2], axis=0)
x11_sum = np.sum(initial_conditions.X[3], axis=0)

# x0_sum = x00_sum * x01_sum
x10_sum = x100_sum * x101_sum
x1_sum = x10_sum * x11_sum

norm = x0_sum @ x1_sum
print(norm)
initial_conditions.Q[0] = 1.0 / norm

# Print tree and write it to a netCDF file
print(tree)
tree.writeTree()

dx = np.prod(initial_conditions.external_nodes[0].grid.n[:2])
P_sum = np.array([np.sum(initial_conditions.X[0][i::dx,0]) for i in range(dx)]) * x1_sum

# X0 = initial_conditions.X[0] * x01_sum
# P_sum = X0 * x1_sum

x = np.arange(0, 46)
y = np.arange(0, 16)
X, Y = np.meshgrid(x, y, indexing="ij")

plt.figure(1)
plt.contour(X, Y, np.reshape(P_sum, (46, 16), order="F"))
plt.show()