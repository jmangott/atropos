import numpy as np 

from scripts.grid_class import GridParms
from scripts.initial_condition_class import InitialCondition
from scripts.tree_class import Tree
from scripts.index_functions import incrVecIndex, vecIndexToState

import scripts.models.toggle_switch as model

# Partition string
partition_str = "(0)(1)"

# Rank
r_out = np.array([5])

# Grid parameters
n = np.array([51, 51])
d = n.size
binsize = np.ones(d, dtype=int)
liml = np.zeros(d)
grid = GridParms(n, binsize, liml)

# Set up the partition tree
tree = Tree(partition_str, grid)
tree.initialize(model.reaction_system, r_out)

# Initial distribution
def gaussian(x: np.ndarray) -> float:
    C = 0.5 * np.array([[75, -15], [-15, 75]])
    Cinv = np.linalg.inv(C)
    mu = np.array([30, 5])
    return np.exp(-0.5 * np.dot(np.transpose(x - mu), np.dot(Cinv, (x - mu))))

p0 = np.zeros(grid.dx())
vec_index = np.zeros(grid.d())
for i in range(grid.dx()):
    state = vecIndexToState(vec_index, grid.liml, grid.binsize)
    p0[i] = gaussian(state + (grid.binsize - 1.0) * 0.5)
    incrVecIndex(vec_index, grid.n, grid.d())

# Normalize p0
p0 = p0 / np.sum(p0)
p0_mat = p0.reshape((tree.root.child[0].grid.dx(), tree.root.child[1].grid.dx()), order="F")

# SVD of p0
u, s, vh = np.linalg.svd(p0_mat, full_matrices=False)

# Use only the first `r` singular values
x0 = u[:, :tree.root.rankOut()]
s = np.diag(s[:tree.root.rankOut()])
x1 = vh[:tree.root.rankOut(), :].T

# Number of basisfunctions
n_basisfunctions = r_out

# Low-rank initial conditions
initial_conditions = InitialCondition(tree, n_basisfunctions)

initial_conditions.Q[0][:, :, 0] = s

initial_conditions.X[0][:] = x0
initial_conditions.X[1][:] = x1

# Print tree and write it to a netCDF file
print(tree)
tree.write()