"""Script for setting the initial conditions for the cascade model."""

import argparse
import sys

import numpy as np

import examples.models.kinetic.cascade as model
from src.grid import GridParms
from src.index_functions import incrVecIndex
from src.initial_condition import InitialCondition
from src.tree import Tree

partition = [
    "((((0 1)(2 3))(4 5))((6 7)(8 9)))(((10 11)(12 13))((14 15)((16 17)(18 19))))",
    "(0 1)((2 3)((4 5)((6 7)((8 9)((10 11)((12 13)((14 15)((16 17)(18 19)))))))))",
    (
        "(((((0)(1))((2)(3)))((4)(5)))(((6)(7))((8)(9))))"
        "((((10)(11))((12)(13)))(((14)(15))(((16)(17))((18)(19)))))"
    ),
    (
        "((0)(1))(((2)(3))(((4)(5))(((6)(7))(((8)(9))"
        "(((10)(11))(((12)(13))(((14)(15))(((16)(17))((18)(19))))))))))"
    ),
]

parser = argparse.ArgumentParser(
    prog="set_cascade",
    usage="python3 examples/kinetic/set_cascade.py --rank 5",
    description="This script sets the initial conditions for the cascade model.",
)

for i, p in enumerate(partition):
    parser.add_argument(
        "-p" + str(i),
        "--partition" + str(i),
        action="store_const",
        const=p,
        required=False,
        help="Set the partition string to " + p,
        dest="partition",
    )
parser.add_argument(
    "-r",
    "--rank",
    nargs="+",
    type=int,
    required=True,
    help="Specify the ranks of the internal nodes",
)
args = parser.parse_args()

if args.partition is None:
    print("usage:", parser.usage)
    print(
        parser.prog + ":",
        "error: the following arguments are required: -p[n]/--partition[n], n=0,...,"
        + str(len(partition) - 1),
    )
    sys.exit(1)

partition_str = args.partition

# Grid parameters
d = 20
n = 64 * np.ones(d, dtype="int")
binsize = np.ones(d, dtype="int")
liml = np.zeros(d)
grid = GridParms(n, binsize, liml)

# Set up the partition tree
tree = Tree(partition_str, grid)

r_out = np.ones(tree.n_internal_nodes, dtype="int") * args.rank
n_basisfunctions = np.ones(r_out.size, dtype="int")
tree.initialize(model.reaction_system, r_out)


def eval_x(x: np.ndarray):
    return 1.0 if np.all(x == 0) else 0.0


# Low-rank initial conditions
initial_conditions = InitialCondition(tree, n_basisfunctions)

for Q in initial_conditions.Q:
    Q[0, 0, 0] = 1.0

idx = 0
for node in range(tree.n_external_nodes):
    vec_index = np.zeros(initial_conditions.external_nodes[node].grid.d())
    for i in range(initial_conditions.external_nodes[node].grid.dx()):
        initial_conditions.X[node][i, :] = eval_x(vec_index)
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
