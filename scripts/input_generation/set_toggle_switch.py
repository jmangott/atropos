"""Script for setting the initial conditions for the toggle switch model."""

import argparse
import numpy as np

from scripts.grid_class import GridParms
from scripts.initial_condition_class import InitialCondition
from scripts.tree_class import Tree
from scripts.index_functions import incrVecIndex

import scripts.models.toggle_switch as model


def constructP0(eval_P0: callable, interval: np.ndarray) -> np.ndarray:
    """Set up the initial probability distribution according to a given function `eval_P0`."""
    dx = np.prod(interval)
    m = interval.size
    P0 = np.zeros(dx)
    vec_index = np.zeros(m)
    for i in range(dx):
        P0[i] = eval_P0(vec_index)
        incrVecIndex(vec_index, interval, m)
    return P0 / np.sum(P0)


parser = argparse.ArgumentParser(
    prog="set_toggle_switch",
    usage="python3 scripts/input_generation/set_toggle_switch.py --rank 5",
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
n_basisfunctions = r_out

# Grid parameters
n = np.array([51, 51])
d = n.size
binsize = np.ones(d, dtype=int)
liml = np.zeros(d)
grid = GridParms(n, binsize, liml)

# Set up the partition tree
tree = Tree(partition_str, grid)
tree.initialize(model.reaction_system, r_out)

# Set up the initial condition
C = 0.5 * np.array([[75, -15], [-15, 75]])
Cinv = np.linalg.inv(C)
mu = np.array([30, 5])


def eval_P0(x: np.ndarray) -> float:
    return np.exp(-0.5 * np.dot(np.transpose(x - mu), np.dot(Cinv, (x - mu))))


P0 = constructP0(eval_P0, n)

u, s, vh = np.linalg.svd(np.reshape(P0, (n[0], n[1]), order="F"), full_matrices=False)

# Low-rank initial conditions
initial_conditions = InitialCondition(tree, n_basisfunctions)

tree.root.child[0].X = u[:, : r_out[0]]
tree.root.child[1].X = vh[: r_out[0], :].T
tree.root.Q[:, :, 0] = np.diag(s[: r_out[0]])

# Print tree and write it to a netCDF file
print(tree)
tree.write()
