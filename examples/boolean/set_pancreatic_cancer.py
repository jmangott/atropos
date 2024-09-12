"""Script for setting the initial conditions for the Boolean pancreatic cancer model."""

import argparse
import numpy as np
import sys

import src.boolean_helper
from src.grid import GridParms
from src.initial_condition import InitialCondition
from src.tree import Tree
from src.index_functions import incrVecIndex

reaction_system = src.boolean_helper.convertRulesToReactions(
    "examples/models/boolean/pancreatic_cancer.hpp"
)

p_best = "((0 1 2 3 4 5 6 9 12)(7 8 10 11 17 21 23 26))((13 14 19 20 22 27 28 31 33)(15 16 18 24 25 29 30 32))"
p_worst = "((0 1 2 4 8 11 17 26)(3 5 6 7 9 10 12 21 23))((13 18 19 20 22 25 31 33)(14 15 16 24 27 28 29 30 32))"
p_reasonable = "((0 1 2 3 4 5 7 9)(13 14 19 20 25 27 29 30 32))((6 10 12 16 18 21 24 26 31)(8 11 15 17 22 23 28 33))"

parser = argparse.ArgumentParser(
    prog="set_pancreatic",
    usage="python3 examples/boolean/set_pancreatic_cancer.py --partition_best --rank 5",
    description="This script sets the initial conditions for the Boolean pancreatic cancer model.",
)

parser.add_argument(
    "-pb",
    "--partition_best",
    action="store_const",
    const=p_best,
    required=False,
    help="Set the partition string to the best partition in terms of entropy",
    dest="partition",
)

parser.add_argument(
    "-pw",
    "--partition_worst",
    action="store_const",
    const=p_worst,
    required=False,
    help="Set the partition string to the worst partition in terms of entropy",
    dest="partition",
)

parser.add_argument(
    "-pr",
    "--partition_reasonable",
    action="store_const",
    const=p_reasonable,
    required=False,
    help="Set the partition string to the best partition in terms of Kerninghan-Lin counts",
    dest="partition",
)

parser.add_argument(
    "-p",
    "--partition",
    type=str,
    required=False,
    help="Specify a general partition string",
    dest="partition",
)

parser.add_argument(
    "-r",
    "--rank",
    type=int,
    required=True,
    help="Specify the ranks of the internal nodes",
)
args = parser.parse_args()

if args.partition is None:
    print("usage:", parser.usage)
    print(
        parser.prog + ":",
        """
          error: one of the following arguments is required:
          -p/--partition`,
          -pb/--partition_best,
          -pw/--partition_worst,
          -pr/--partition_reasonable,
          """,
    )
    sys.exit(1)

partition_str = args.partition

# Grid parameters
d = 34
n = 2 * np.ones(d, dtype=int)
binsize = np.ones(d, dtype=int)
liml = np.zeros(d)
grid = GridParms(n, binsize, liml)

# Set up the partition tree
tree = Tree(partition_str, grid)

r_out = np.ones(tree.n_internal_nodes, dtype="int") * args.rank
n_basisfunctions = np.ones(r_out.size, dtype="int")
tree.initialize(reaction_system, r_out)


def eval_x(x: np.ndarray, grid: GridParms):
    result = 1.0 / grid.dx()
    # pos0 = np.argwhere(grid.species==0) # HMGB
    # pos4 = np.argwhere(grid.species==4) # RAS
    # pos25 = np.argwhere(grid.species==25) # P54
    # if pos0.size > 0:
    #     result *= (1.0 if x[pos0] == 1 else 0.0)
    # if pos4.size > 0:
    #     result *= (1.0 if x[pos4] == 1 else 0.0)
    # if pos25.size > 0:
    #     result *= (1.0 if x[pos25] == 0 else 0.0)
    return result


# Low-rank initial conditions
initial_conditions = InitialCondition(tree, n_basisfunctions)

for Q in initial_conditions.Q:
    Q[0, 0, 0] = 1.0

for n_node in range(tree.n_external_nodes):
    vec_index = np.zeros(initial_conditions.external_nodes[n_node].grid.d())
    for i in range(initial_conditions.external_nodes[n_node].grid.dx()):
        initial_conditions.X[n_node][i, :] = eval_x(
            vec_index, initial_conditions.external_nodes[n_node].grid
        )
        incrVecIndex(
            vec_index,
            initial_conditions.external_nodes[n_node].grid.n,
            initial_conditions.external_nodes[n_node].grid.d(),
        )

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
