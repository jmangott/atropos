"""Script for setting the initial conditions for the enzymatic futile cycle model."""
import argparse
import numpy as np
import sys

from scripts.grid_class import GridParms
from scripts.initial_condition_class import InitialCondition
from scripts.tree_class import Tree
from scripts.index_functions import incrVecIndex

import scripts.models.enzymatic_futile_cycle as model

parser = argparse.ArgumentParser(
                    prog='set_enzymatic_futile_cycle',
                    usage='python3 scripts/input/set_enzymatic_futile_cycle.py --rank 5',
                    description='This script sets the initial conditions for the enzymatic futile cycle model.')

parser.add_argument('-r',
                    '--rank',
                    nargs='+',
                    type=int,
                    required=True,
                    help="Specify the ranks of the internal nodes",
                    )
args = parser.parse_args()

partition_str = "(0 1 2)(3 4 5)"
r_out = np.array(args.rank)
n_basisfunctions = np.ones(r_out.size, dtype="int")

# Grid parameters
n = np.array([128, 4, 4, 128, 4, 4])
d = n.size
binsize = np.ones(d, dtype=int)
liml = np.zeros(d)
grid = GridParms(n, binsize, liml)

# Set up the partition tree
tree = Tree(partition_str, grid)
tree.initialize(model.reaction_system, r_out)

mu = np.array([30, 2, 0, 90, 2, 0])

def eval_x(x: np.ndarray, mu: np.ndarray):
    return 1.0 if np.all(x == mu) else 0.0

# Low-rank initial conditions
initial_conditions = InitialCondition(tree, n_basisfunctions)

for Q in initial_conditions.Q:
    Q[0, 0, 0] = 1.0

idx = 0
mu_perm = mu[tree.species]
for node in range(tree.n_external_nodes):
    vec_index = np.zeros(initial_conditions.external_nodes[node].grid.d())
    for i in range(initial_conditions.external_nodes[node].grid.dx()):
        initial_conditions.X[node][i, :] = eval_x(vec_index, mu_perm[idx : idx+len(vec_index)])
        incrVecIndex(vec_index, initial_conditions.external_nodes[node].grid.n, initial_conditions.external_nodes[node].grid.d())
    idx += len(vec_index)

# Calculate norm
_, marginal_distribution = tree.calculateObservables(np.zeros(tree.root.grid.d(), dtype="int"))
norm = np.sum(marginal_distribution[0])

# import matplotlib.pyplot as plt
# plt.plot(np.arange(tree.grid.n[0]), marginal_distribution[0], label="$x_0$")
# plt.plot(np.arange(tree.grid.n[1]), marginal_distribution[1], label="$x_1$")
# plt.plot(np.arange(tree.grid.n[2]), marginal_distribution[2], label="$x_2$")
# plt.plot(np.arange(tree.grid.n[3]), marginal_distribution[3], label="$x_3$")
# plt.plot(np.arange(tree.grid.n[4]), marginal_distribution[4], label="$x_4$")
# plt.plot(np.arange(tree.grid.n[5]), marginal_distribution[5], label="$x_5$")
# plt.legend()
# plt.show()

print("norm:", norm)
tree.root.Q[0, 0, 0] /= norm

# Print tree and write it to a netCDF file
print(tree)
tree.write()