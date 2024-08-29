"""Script for setting the initial conditions for the Boolean apoptosis model."""
import argparse
import numpy as np
import sys

import scripts.boolean_helper
from scripts.grid_class import GridParms
from scripts.initial_condition_class import InitialCondition
from scripts.tree_class import Tree
from scripts.index_functions import incrVecIndex

reaction_system = scripts.boolean_helper.convertRulesToReactions("scripts/models/boolean_rulefiles/apoptosis.hpp")

p_best = '((0 2 3 4 5 6 7 8 9 12 20)(21 25 26 27 28 29 30 31 39 40))((1 10 11 14 15 22 23 32 33 38)(13 16 17 18 19 24 34 35 36 37))'
p_worst = '((0 1 4 7 11 12 14 16 23 33 34)(2 3 5 8 13 19 20 24 32 35))((6 9 15 17 21 22 25 26 31 37)(10 18 27 28 29 30 36 38 39 40))'
p_reasonable = '((0 2 3 4 5 6 7 8 12 20)(1 9 10 11 14 15 16 17 22 23))((13 18 19 24 32 33 34 35 36 37 38)(21 25 26 27 28 29 30 31 39 40))'

parser = argparse.ArgumentParser(
                    prog='set_bax',
                    usage='python3 scripts/input_generation/set_boolean_pancreatic_cancer.py --parition_best --rank 5',
                    description='This script sets the initial conditions for the Boolean apoptosis model.')

parser.add_argument('-p_best', 
                    '--partition_best', 
                    action='store_const', 
                    const=p_best,
                    required=False, 
                    help='Set the partition string to the best partition in terms of entropy',
                    dest='partition', 
                    )

parser.add_argument('-p_worst',
                    '--partition_worst', 
                    action='store_const', 
                    const=p_worst,
                    required=False, 
                    help='Set the partition string to the worst partition in terms of entropy',
                    dest='partition', 
                    )

parser.add_argument('-p_reasonable',
                    '--partition_reasonable', 
                    action='store_const', 
                    const=p_reasonable,
                    required=False, 
                    help='Set the partition string to the best partition in terms of Kerninghan-Lin counts',
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
    print(parser.prog+":",
          """
          error: one of the following arguments is required:
          -p/--partition`,
          -p_best/--partition_best,
          -p_worst/--partition_worst,
          -p_reasonable/--partition_reasonable,
          """)
    sys.exit(1)

partition_str = args.partition

# Grid parameters
d = 41
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
    pos0 = np.argwhere(grid.species==0) # TNF
    pos1 = np.argwhere(grid.species==1) # GF
    pos40 = np.argwhere(grid.species==40) # DNAdam
    if pos0.size > 0:
        result *= (1.0 if x[pos0] == 1 else 0.0)
    if pos1.size > 0:
        result *= (1.0 if x[pos1] == 1 else 0.0)
    if pos40.size > 0:
        result *= (1.0 if x[pos40] == 0 else 0.0)
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
norm = np.sum(marginal_distribution[tree.species_names[0]])
print("norm:", norm)
tree.root.Q[0, 0, 0] /= norm

# Print tree and write it to a netCDF file
print(tree)
tree.write()