import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial

from initial_condition_class import *

import models.lambda_phage as model

def tensorUnfold(tensor, mode):
    """
    Cf. https://stackoverflow.com/questions/49970141/using-numpy-reshape-to-perform-3rd-rank-tensor-unfold-operation
    """
    return np.reshape(np.moveaxis(tensor, mode, 0), (tensor.shape[mode], -1), order="F")

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
def multinomial(x):
    abs_x = np.sum(x)
    if (abs_x <= 3):
        p0 = factorial(3) * (0.05 ** abs_x) * \
            ((1.0 - 5 * 0.05) ** (3 - abs_x)) / (np.prod(factorial(x)) * factorial(3 - abs_x))
    else:
        p0 = 0.0
    return p0

p0 = np.zeros(grid.dx)
vec_index = np.zeros(grid.d)
for i in range(grid.dx):
    state = vecIndexToState(vec_index, grid.liml, grid.binsize)
    p0[i] = multinomial(state + (grid.binsize - 1.0) * 0.5)
    incrVecIndex(vec_index, grid.n, grid.d)

p0_mat = p0.reshape((tree.root.child[0].grid.dx, tree.root.child[1].grid.dx), order="F")

# SVD of p0
u, s, vh = np.linalg.svd(p0_mat, full_matrices=False)

# Use only the first `r` singular values
# Transpose x1 and x2, as Ensign works with column-major order arrays
x0 = u[:, :tree.root.r_out].T
q = np.diag(s[:tree.root.r_out])
x1 = vh[:tree.root.r_out, :]

# SVD of x0
q0 = np.zeros((r_out[0], r_out[1], r_out[1]))
x0_tensor = np.zeros((tree.root.r_out, tree.root.child[0].child[0].grid.dx, tree.root.child[0].child[1].grid.dx))

for i in range(tree.root.r_out):
    x0_tensor[i, :, :] = x0[i].reshape((tree.root.child[0].child[0].grid.dx, tree.root.child[0].child[1].grid.dx), order="F")

x0_mat = tensorUnfold(x0_tensor, 1)
u, _, _ = np.linalg.svd(x0_mat, full_matrices=False)
x00 = u[:, :tree.root.child[0].r_out].T

x0_mat = tensorUnfold(x0_tensor, 2)
u, _, _ = np.linalg.svd(x0_mat, full_matrices=False)
x01 = u[:, :tree.root.child[0].r_out].T

q0 = np.einsum('ijk,lj,mk', x0_tensor, x00, x01)

# Number of basisfunctions
n_basisfunctions = r_out

# Low-rank initial conditions
initial_conditions = InitialCondition(tree, n_basisfunctions)

initial_conditions.Q[0][0, :, :] = q

q0T = np.zeros(q0.shape)
for i in range(tree.root.r_out):
    q0T[i, :, :] = np.transpose(q0[i, :, :])

initial_conditions.Q[1][:] = q0T

initial_conditions.X[0][:] = x00
initial_conditions.X[1][:] = x01
initial_conditions.X[2][:] = x1

# Print tree and write it to a netCDF file
print(tree)
tree.writeTree()

x00_sum = np.sum(x00, axis=1)
x01_sum = np.sum(x01, axis=1)
x1_sum = np.sum(x1, axis=1)
x0_sum = np.array([x00_sum @ q0[i, :, :] @ x01_sum.T for i in range(r_out[0])])
norm = x0_sum @ q @ x1_sum.T
print(norm)
print(q0)
print(x1)

X0 = np.array([x00.T @ q0[i, :, :] @ x01_sum for i in range(q0.shape[0])])
P_sum = X0.T @ q @ x1_sum

x = np.arange(0, 16)
y = np.arange(0, 41)
X, Y = np.meshgrid(x, y, indexing="ij")

plt.figure(1)
plt.contour(X, Y, np.reshape(P_sum, (16, 41), order="F"))
plt.show()