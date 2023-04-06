"""
Script for setting parameters (`parameters.hpp`) and initial conditions (`input.nc`) 
for the Tyson oscillator model.
Call in `input` folder with: `python3 example_setups/set_tyson.py --flags`
"""

import numpy as np
import sys

sys.path.append('../input')

from input_helper import *
from parameters_helper import ParametersParser

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
params.kModel = "reactions_tyson"
r = params.kR = 6
m1 = params.kM1 = 3
m2 = params.kM2 = 2
n1 = params.kN1 = (21, 61, 11)
n2 = params.kN2 = (71, 61)
bin1 = params.kBinsize1 = (1, 1, 1)
bin2 = params.kBinsize2 = (1, 1)
liml1 = params.kLiml1 = (0, 0, 0)
liml2 = params.kLiml2 = (0, 0)

x10 = (15, 0, 0)
x20 = (60, 0)

# Change `parameters.hpp`
params.configure()

# Set up `grid_info`
grid = grid_info(r, m1, m2, n1, n2, bin1, bin2, liml1, liml2)

# Set input/initial conditions
SetInputKD(x10, x20, grid)
