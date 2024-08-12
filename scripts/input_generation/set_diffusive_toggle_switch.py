"""Script for setting the initial conditions for the diffusive toggle switch model."""
import argparse
import numpy as np
import sys

from scripts.grid_class import GridParms
from scripts.initial_condition_class import InitialCondition
from scripts.tree_class import Tree
from scripts.index_functions import incrVecIndex

import scripts.models.diffusive_toggle_switch as model

partition = ['(((0)(1))((2)(3)))(((4)(5))((6)(7)))', 
             '(0)((1)((2)((3)((4)((5)((6)(7)))))))']

parser = argparse.ArgumentParser(
                    prog='set_diffusive_toggle_switch',
                    usage='python3 scripts/input_generation/set_diffusive_toggle_switch.py --partition0 --rank 5',
                    description='This script sets the initial conditions for the diffusive toggle switch model.')

for i, p in enumerate(partition):
    parser.add_argument('-p'+str(i), 
                        '--partition'+str(i), 
                        action='store_const', 
                        const=p,
                        required=False, 
                        help='Set the partition string to '+p,
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
    print(parser.prog+":", "error: the following arguments are required: -p[n]/--partition[n], n=0,...,"+str(len(partition)-1))
    sys.exit(1)

partition_str = args.partition
r_out = np.ones(7, dtype="int") * args.rank
n_basisfunctions = np.ones(r_out.size, dtype="int")

# Grid parameters
n = np.array([51, 51, 51, 51, 51, 51, 51, 51])
d = n.size
binsize = np.ones(d, dtype=int)
liml = np.zeros(d)
grid = GridParms(n, binsize, liml)

# Set up the partition tree
tree = Tree(partition_str, grid)
tree.initialize(model.reaction_system, r_out)

C = 0.2
Cinv = 1 / C
mu = np.array([30, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0])

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
        initial_conditions.X[node][i, :] = eval_x(vec_index, mu_perm[idx : idx+len(vec_index)])
        incrVecIndex(vec_index, initial_conditions.external_nodes[node].grid.n, initial_conditions.external_nodes[node].grid.d())
    idx += len(vec_index)

# Calculate norm
_, marginal_distribution = tree.calculateObservables(np.zeros(tree.root.grid.d(), dtype="int"))
norm = np.sum(marginal_distribution[0])
print("norm:", norm)
tree.root.Q[0, 0, 0] /= norm

# Print tree and write it to a netCDF file
print(tree)
tree.write()