"""
Script for setting parameters (`parameters.hpp`) and initial conditions (`input.nc`) 
for the BAX pore aggregation model.
Call in project root with: `python3 scripts/input/example_setups/set_bax.py --flags`
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
params.kModel = "reactions_bax_seq"
r = params.kR = 4
m1 = params.kM1 = 7
m2 = params.kM2 = 4
n1 = params.kN1 = (6, 6, 6, 6, 6, 21, 6)
n2 = params.kN2 = (91, 36, 21, 11)
bin1 = params.kBinsize1 = (1, 1, 1, 1, 1, 1, 1)
bin2 = params.kBinsize2 = (1, 1, 1, 1)
liml1 = params.kLiml1 = (0, 0, 0, 0, 0, 90, 0)
liml2 = params.kLiml2 = (0, 0, 0, 0)

x10 = (0, 0, 0, 0, 0, 100, 0)
x20 = (80, 0, 0, 0)

# Change `parameters.hpp`
params.configure()

# Set up `GridInfo`
grid = GridInfo(r, m1, m2, n1, n2, bin1, bin2, liml1, liml2)

# Set input/initial conditions
SetInputKD(x10, x20, grid)
