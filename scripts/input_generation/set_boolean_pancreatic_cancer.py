"""Script for setting the initial conditions for the Boolean pancreatic cancer model."""
import argparse
import numpy as np
import sys

from scripts.grid_class import GridParms
from scripts.initial_condition_class import InitialCondition
from scripts.tree_class import Tree
from scripts.index_functions import incrVecIndex

import scripts.models.boolean_pancreatic_cancer as model

partition = ['(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)(17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33)',
             '((0 1 2 3 4 5 6 7 8)(9 10 11 12 13 14 15 16))((17 18 19 20 21 22 23 24 25)(26 27 28 29 30 31 32 33))',
             '(((0 1 2 3 4)(5 6 7 8))(((9 10 11 12)(13 14 15 16)))(((17 18 19 20 21)(22 23 24 25))((I(26 27 28 29)(30 31 32 33)))']

parser = argparse.ArgumentParser(
                    prog='set_bax',
                    usage='python3 scripts/input_generation/set_boolean_pancreatic_cancer.py --rank 5',
                    description='This script sets the initial conditions for the Boolean pancreatic cancer model.')

for i, p in enumerate(partition):
    parser.add_argument('-p'+str(i), 
                        '--partition'+str(i), 
                        action='store_const', 
                        const=p,
                        required=False, 
                        help='Set the partition string to '+p,
                        dest='partition', 
                        )
parser.add_argument('-p', 
                    '--partition', 
                    type=str, 
                    required=False, 
                    help='Specify a general partition string',
                    dest='partition', 
                    )

parser.add_argument('-r', 
                    '--rank', 
                    type=int, 
                    required=True, 
                    help="Specify the ranks of the internal nodes",
                    )
args = parser.parse_args()

if args.partition == None:
    print("usage:", parser.usage)
    print(parser.prog+":", "error: the following arguments are required: -p/--partition` or -p[n]/--partition[n], n=0,...,"+str(len(partition)-1))
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
tree.initialize(model.reaction_system, r_out)

def eval_x(x: np.ndarray, grid: GridParms):
    result = 1.0 / grid.dx()
    pos0 = np.argwhere(grid.species==0) # HMGB
    pos4 = np.argwhere(grid.species==4) # RAS
    pos25 = np.argwhere(grid.species==25) # P54
    if pos0.size > 0:
        result *= (1.0 if x[pos0] == 1 else 0.0)
    if pos4.size > 0:
        result *= (1.0 if x[pos4] == 1 else 0.0)
    if pos25.size > 0:
        result *= (1.0 if x[pos25] == 0 else 0.0)
    return result

# Low-rank initial conditions
initial_conditions = InitialCondition(tree, n_basisfunctions)

for Q in initial_conditions.Q:
    Q[0, 0, 0] = 1.0

for n_node in range(tree.n_external_nodes):
    vec_index = np.zeros(initial_conditions.external_nodes[n_node].grid.d())
    for i in range(initial_conditions.external_nodes[n_node].grid.dx()):
        initial_conditions.X[n_node][i, :] = eval_x(vec_index, initial_conditions.external_nodes[n_node].grid)
        incrVecIndex(vec_index, initial_conditions.external_nodes[n_node].grid.n, initial_conditions.external_nodes[n_node].grid.d())

# Calculate norm
_, marginal_distribution = tree.calculateObservables(np.zeros(tree.root.grid.d(), dtype="int"))
norm = np.sum(marginal_distribution[0])
print("norm:", norm)
tree.root.Q[0, 0, 0] /= norm

# Print tree and write it to a netCDF file
print(tree)
tree.write()