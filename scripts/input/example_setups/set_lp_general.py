"""
Script for setting parameters (`parameters.hpp`) and initial conditions (`input.nc`) 
for the lambda phage model with initial conditions according to Jahnke & Huisinga 2008.
Call in project root with: `python3 scripts/input/example_setups/set_lp_general.py --flags`
"""

import numpy as np
from scipy.special import factorial

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
n2 = params.kN2 = (11, 11, 11)
bin1 = params.kBinsize1 = (1, 1)
bin2 = params.kBinsize2 = (1, 1, 1)
liml1 = params.kLiml1 = (0, 0)
liml2 = params.kLiml2 = (0, 0, 0)

def eval_p0(x):
    abs_x = np.sum(x)
    if (abs_x <= 3):
        p0 = factorial(3) * (0.05 ** abs_x) * \
            ((1.0 - 5 * 0.05) ** (3 - abs_x)) / (np.prod(factorial(x)) * factorial(3 - abs_x))
    else:
        p0 = 0.0
    return p0

# Change `parameters.hpp`
params.configure()

# Set up `GridInfo`
grid = GridInfo(r, m1, m2, n1, n2, bin1, bin2, liml1, liml2)

# Set input/initial conditions
SetInputGeneral(eval_p0, grid)
