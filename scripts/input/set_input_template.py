"""
Template for setting parameters (`parameters.hpp`) and initial conditions (`input.nc`).
Modify code marked with TODO according to your needs.
Call in project root with: `python3 scripts/input/set_input_template.py --flags` (use --h for help)
"""

import numpy as np
from scipy.special import factorial

from scripts.input.input_helper import *
from scripts.input.parameters_helper import ParametersParser

parser = ParametersParser(
    description="Configure parameters."
)

# Add arguments for setting parameters with `ParametersParser``
# TODO: add or remove flags according to your setup
parser.add_tstar()
parser.add_tau()
parser.add_snapshot()
parser.add_secondorder()
parser.add_substeps()
parser.add_filename()
params = parser.parse_args()

# Set in-script parameters
# TODO: adjust parameters according to your model
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

# TODO: set/modify `x10` and `x20` when using `SetInputKD(x10, x20, grid)`
x10 = (0, 85)
x20 = (0, 5, 10)

#TODO: set/modify `eval_p0(x)` when using `SetInputGeneral(eval_p0, grid)`
def eval_p0(x):
    abs_x = np.sum(x)
    if (abs_x <= 3):
        p0 = factorial(3) * (0.05 ** abs_x) * ((1.0 - 5 * 0.05) ** (3 - abs_x)) / (np.prod(factorial(x)) * factorial(3 - abs_x))
    else:
        p0 = 0.0
    return p0

# TODO: set/modify `eval_x10(x)`, `eval_x20(x)` and `s`
# when using `SetInputGeneral(eval_x10, eval_x20, s, grid)`
C = 0.2
Cinv = 1 / C
mu = np.array([0, 0, 0, 0, 0, 100, 0, 80, 0, 0, 0])

def eval_x10(x: np.ndarray) -> float:
    return np.exp(-0.5 * np.dot(np.transpose(x - mu[:m1]), np.dot(Cinv, (x - mu[:m1]))))

def eval_x20(x: np.ndarray) -> float:
    return np.exp(-0.5 * np.dot(np.transpose(x - mu[m1:]), np.dot(Cinv, (x - mu[m1:]))))

s = [1]

# Change `parameters.hpp`
params.configure()

# Set up `GridInfo`
grid = GridInfo(r, m1, m2, n1, n2, bin1, bin2, liml1, liml2)

# Set input/initial conditions
# TODO: choose either `SetInputKD`, `SetInputFunc` or `SetInputGeneral`
setInputKD(x10, x20, grid)
# setInputFunc([eval_x10], [eval_x20], s, grid)
# setInputGeneral(eval_p0, grid)
