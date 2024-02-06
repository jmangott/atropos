import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial

from scripts.grid_class import GridParms
from scripts.initial_condition_class import InitialCondition
from scripts.tree_class import Tree
from scripts.index_functions import incrVecIndex, vecIndexToState, tensorUnfold

import scripts.hierarchical.models.lambda_phage as model

# Partition string
partition_str = "((0 1)(2 4))(3)"

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
def multinomial(x):
    abs_x = np.sum(x)
    if (abs_x <= 3):
        p0 = factorial(3) * (0.05 ** abs_x) * \
            ((1.0 - 5 * 0.05) ** (3 - abs_x)) / (np.prod(factorial(x)) * factorial(3 - abs_x))
    else:
        p0 = 0.0
    return p0

p0 = np.zeros(grid.dx())
vec_index = np.zeros(grid.d())
for i in range(grid.dx()):
    state = vecIndexToState(vec_index, grid.liml, grid.binsize)
    p0[i] = multinomial(state + (grid.binsize - 1.0) * 0.5)
    incrVecIndex(vec_index, grid.n, grid.d())

p0_mat = p0.reshape((tree.root.child[0].grid.dx(), tree.root.child[1].grid.dx()), order="F")

# SVD of p0
u, s, vh = np.linalg.svd(p0_mat, full_matrices=False)

# Use only the first `r` singular values
x0 = u[:, :tree.root.r_out]
q = np.diag(s[:tree.root.r_out])
x1 = vh[:tree.root.r_out, :].T

# SVD of x0
q0 = np.zeros((r_out[1], r_out[1], r_out[0]))
x0_tensor = np.zeros((tree.root.child[0].child[0].grid.dx(), tree.root.child[0].child[1].grid.dx(), tree.root.r_out))

for i in range(tree.root.r_out):
    x0_tensor[:, :, i] = x0[:, i].reshape((tree.root.child[0].child[0].grid.dx(), tree.root.child[0].child[1].grid.dx()), order="F")

x0_mat = tensorUnfold(x0_tensor, 0)
u, _, _ = np.linalg.svd(x0_mat, full_matrices=False)
x00 = u[:, :tree.root.child[0].r_out]

x0_mat = tensorUnfold(x0_tensor, 1)
u, _, _ = np.linalg.svd(x0_mat, full_matrices=False)
x01 = u[:, :tree.root.child[0].r_out]

q0 = np.einsum('ik,jl,ijm', x00, x01, x0_tensor)

# Number of basisfunctions
n_basisfunctions = r_out

# Low-rank initial conditions
initial_conditions = InitialCondition(tree, n_basisfunctions)

initial_conditions.Q[0][:, :, 0] = q
initial_conditions.Q[1][:] = q0

initial_conditions.X[0][:] = x00
initial_conditions.X[1][:] = x01
initial_conditions.X[2][:] = x1

# Print tree and write it to a netCDF file
print(tree)
tree.writeTree()

x00_sum = np.sum(x00, axis=0)
x01_sum = np.sum(x01, axis=0)
x1_sum = np.sum(x1, axis=0)
x0_sum = np.array([x00_sum @ q0[:, :, i] @ x01_sum.T for i in range(r_out[0])])
norm = x0_sum @ q @ x1_sum.T
print(norm)

X0 = np.array([x00 @ q0[:, :, i] @ x01_sum.T for i in range(q0.shape[-1])]).T
P_sum = X0 @ q @ x1_sum.T

x = np.arange(0, 16)
y = np.arange(0, 41)
X, Y = np.meshgrid(x, y, indexing="ij")

plt.figure(1)
plt.contour(X, Y, np.reshape(P_sum, (16, 41), order="F"))
plt.show()