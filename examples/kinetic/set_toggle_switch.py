"""Script for setting the initial conditions for the toggle switch model."""

import argparse
import numpy as np

from src.grid import GridParms
from src.initial_condition import InitialCondition
from src.tree import Tree
from src.index_functions import incrVecIndex

import examples.models.kinetic.toggle_switch as model

parser = argparse.ArgumentParser(
    prog="set_toggle_switch",
    usage="python3 examples/kinetic/set_toggle_switch.py --rank 5",
    description="This script sets the initial conditions for the toggle switch model.",
)

parser.add_argument(
    "-r",
    "--rank",
    type=int,
    required=True,
    help="Specify the ranks of the internal nodes",
)

args = parser.parse_args()

partition_str = "(0)(1)"
r_out = np.array([args.rank])
n_basisfunctions = np.ones(r_out.size, dtype="int")

# Grid parameters
n = np.array([51, 51])
d = n.size
binsize = np.ones(d, dtype=int)
liml = np.zeros(d)
grid = GridParms(n, binsize, liml)

# Set up the partition tree
tree = Tree(partition_str, grid)
tree.initialize(model.reaction_system, r_out)

C = 0.2
Cinv = 1 / C
mu = np.array([30, 5])


def eval_x(x: np.ndarray, mu: np.ndarray):
    return np.exp(-0.5 * Cinv * np.dot(np.transpose(x - mu), (x - mu)))


# Low-rank initial conditions
initial_conditions = InitialCondition(tree, n_basisfunctions)

for Q in initial_conditions.Q:
    Q[0, 0, 0] = 1.0

idx = 0
mu_perm = mu[tree.species]
for node in range(tree.n_external_nodes):
    vec_index = np.zeros(initial_conditions.external_nodes[node].grid.d())
    for i in range(initial_conditions.external_nodes[node].grid.dx()):
        initial_conditions.X[node][i, :] = eval_x(
            vec_index, mu_perm[idx : idx + len(vec_index)]
        )
        incrVecIndex(
            vec_index,
            initial_conditions.external_nodes[node].grid.n,
            initial_conditions.external_nodes[node].grid.d(),
        )
    idx += len(vec_index)

# Calculate norm
_, marginal_distribution = tree.calculateObservables(
    np.zeros(tree.root.grid.d(), dtype="int")
)
norm = np.sum(marginal_distribution[tree.species_names[0]])
print("norm:", norm)
tree.root.Q[0, 0, 0] /= norm

# Print tree and write it to a netCDF file
print(tree)
tree.write()
