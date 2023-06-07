"""
Script for setting parameters (`parameters.hpp`) and initial conditions (`input.nc`) 
for the lambda phage model with initial conditions of Kronecker delta form.
Call in project root with: `python3 scripts/input/example_setups/set_lp_kd.py --flags`
"""

from scripts.input.input_helper import *
from scripts.input.parameters_helper import ParametersParser

parser = ParametersParser(
    description="Configure parameters."
)

# Add arguments for setting parameters with `ParametersParser``
parser.add_tstar()
parser.add_tau()
parser.add_snapshot()
parser.add_secondorder()
parser.add_substeps()
parser.add_filename()
params = parser.parse_args()

# Set in-script parameters
params.kModel = "reactions_lp"
r = params.kR = 4
m1 = params.kM1 = 2
m2 = params.kM2 = 3
n1 = params.kN1 = (16, 41)
n2 = params.kN2 = (6, 11, 21)
bin1 = params.kBinsize1 = (1, 1)
bin2 = params.kBinsize2 = (1, 1, 1)
liml1 = params.kLiml1 = (0, 70)
liml2 = params.kLiml2 = (0, 0, 0)

x10 = (0, 85)
x20 = (0, 5, 10)

# Change `parameters.hpp`
params.configure()

# Set up `grid_info`
grid = grid_info(r, m1, m2, n1, n2, bin1, bin2, liml1, liml2)

# Set input/initial conditions
SetInputKD(x10, x20, grid)
