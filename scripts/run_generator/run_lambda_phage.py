from scripts.generator_class import Model, Partitioning, run
import numpy as np
import sympy as sp

from scipy.special import factorial
from scripts.index_functions import incrVecIndex, vecIndexToState, tensorUnfold

# Partitition p0, rank 5 4
# Snapshot 100, tau 1e-3, final_time 1, rest default

"""
Define variables
"""
S1, S2, S3, S4, S5 = sp.symbols("S1, S2, S3, S4, S5")


"""
Generate model
"""
model = Model((S1, S2, S3, S4, S5))

kA0 = 0.5
kA1 = 1.0
kA2 = 0.15
kA3 = 0.3
kA4 = 0.3
kB0 = 0.12
kB1 = 0.6
kB2 = 1.0
kB3 = 1.0
kB4 = 1.0
kC0 = 0.0025
kC1 = 0.0007
kC2 = 0.0231
kC3 = 0.01
kC4 = 0.01

### Full reaction generation:

#model.add_reaction(0, S1, {S2: kA0 * kB0 / (kB0 + S2)})
#model.add_reaction(0, S2, {S1: kB1 / (kB1 + S1), S5: kA1 + S5})
#model.add_reaction(0, S3, {S2: kA2 * kB2 * S2 / (kB2 * S2 + 1.0)})
#model.add_reaction(0, S4, {S3: kA3 * kB3 * S3 / (kB3 * S3 + 1.0)})
#model.add_reaction(0, S5, {S3: kA4 * kB4 * S3 / (kB4 * S3 + 1.0)})
#model.add_reaction(S1, 0, {S1: kC0 * S1})
#model.add_reaction(S2, 0, {S2: kC1 * S2})
#model.add_reaction(S3, 0, {S3: kC2 * S3})
#model.add_reaction(S4, 0, {S4: kC3 * S4})
#model.add_reaction(S5, 0, {S5: kC4 * S5})


### Shorter reaction generation:

model.add_reaction(0, S1, {S2: kA0 * kB0 / (kB0 + S2)})
model.add_reaction(0, S2, {S1: kB1 / (kB1 + S1), S5: kA1 + S5})
model.add_reaction(0, S3, {S2: kA2 * kB2 * S2 / (kB2 * S2 + 1.0)})
model.add_reaction(0, S4, {S3: kA3 * kB3 * S3 / (kB3 * S3 + 1.0)})
model.add_reaction(0, S5, {S3: kA4 * kB4 * S3 / (kB4 * S3 + 1.0)})
model.add_reaction(S1, 0, kC0)
model.add_reaction(S2, 0, kC1)
model.add_reaction(S3, 0, kC2)
model.add_reaction(S4, 0, kC3)
model.add_reaction(S5, 0, kC4)

model.generate_reaction_system()



"""
Generate tree and initial condition
"""
r = np.array([5,4])
p0 = '(S1 S2)((S3 S4)(S5))'
partitioning = Partitioning(p0, r, model)

n = np.array([16, 41, 11, 11, 11])
d = n.size
binsize = np.ones(d, dtype=int)
liml = np.zeros(d)
partitioning.add_grid_params(n, binsize, liml)

partitioning.generate_tree()



# Initial distribution
def multinomial(x):
    abs_x = np.sum(x)
    if (abs_x <= 3):
        p0 = factorial(3) * (0.05 ** abs_x) * \
            ((1.0 - 5 * 0.05) ** (3 - abs_x)) / (np.prod(factorial(x)) * factorial(3 - abs_x))
    else:
        p0 = 0.0
    return p0

# Helper function for factorizing a low-rank factor
def factorizeLRFactor(x, node):
    x_tensor = np.zeros((node.child[0].grid.dx(), node.child[1].grid.dx(), node.rankIn()))
    for i in range(node.rankIn()):
        x_tensor[:, :, i] = x[:, i].reshape((node.child[0].grid.dx(), node.child[1].grid.dx()), order="F")

    x_mat = tensorUnfold(x_tensor, 0)
    u, _, _ = np.linalg.svd(x_mat, full_matrices=False)
    x0 = u[:, :node.child[0].rankIn()]

    x_mat = tensorUnfold(x_tensor, 1)
    u, _, _ = np.linalg.svd(x_mat, full_matrices=False)
    x1 = u[:, :node.child[1].rankIn()]

    q = np.einsum('ik,jl,ijm', x0, x1, x_tensor)

    return q, x0, x1

p0 = np.zeros(partitioning.grid.dx())
vec_index = np.zeros(partitioning.grid.d())
for i in range(partitioning.grid.dx()):
    state = vecIndexToState(vec_index, partitioning.grid.liml, partitioning.grid.binsize)
    p0[i] = multinomial(state + (partitioning.grid.binsize - 1.0) * 0.5)
    incrVecIndex(vec_index, partitioning.grid.n, partitioning.grid.d())

p0_mat = p0.reshape((partitioning.tree.root.child[0].grid.dx(), partitioning.tree.root.child[1].grid.dx()), order="F")

# SVD of p0
u, s, vh = np.linalg.svd(p0_mat, full_matrices=False)

# Use only the first `r` singular values
x0 = u[:, :partitioning.tree.root.rankOut()]
q = np.diag(s[:partitioning.tree.root.rankOut()])
x1 = vh[:partitioning.tree.root.rankOut(), :].T

# SVD of x0
q1, x10, x11 = factorizeLRFactor(x1, partitioning.tree.root.child[1])       # This because p0

# Number of basisfunctions
n_basisfunctions = r

# Low-rank initial conditions
partitioning.generate_initial_condition(n_basisfunctions)

partitioning.initial_conditions.Q[0][:, :, 0] = q
partitioning.initial_conditions.Q[1][:] = q1
partitioning.initial_conditions.X[0][:] = x0
partitioning.initial_conditions.X[1][:] = x10
partitioning.initial_conditions.X[2][:] = x11

x10_sum = np.sum(x10, axis=0)
x11_sum = np.sum(x11, axis=0)
x0_sum = np.sum(x0, axis=0)
x1_sum = np.array([x10_sum @ q1[:, :, i] @ x11_sum.T for i in range(r[0])])

norm = x0_sum @ q @ x1_sum.T
print("norm:", norm)



print(partitioning.tree)



"""
write input file and run
"""
run(partitioning, 'lambda_phage_Stefan', 10, 1e-3, 1, 1, "implicit_Euler")
