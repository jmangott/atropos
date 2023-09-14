"""
Script for setting parameters (`parameters.hpp`) and initial conditions (`input.nc`) 
for the toggle switch model with initial conditions according to Jahnke & Huisinga 2008.
Call in project root with: `python3 scripts/input/example_setups/set_ts.py --flags`
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
params.kModel = "reactions_ts"
r = params.kR = 5
m1 = params.kM1 = 1
m2 = params.kM2 = 1
n1 = params.kN1 = (51,)
n2 = params.kN2 = (51,)
bin1 = params.kBinsize1 = (1,)
bin2 = params.kBinsize2 = (1,)
liml1 = params.kLiml1 = (0,)
liml2 = params.kLiml2 = (0,)

def eval_p0(x):
    C = 0.5 * np.array([[75, -15], [-15, 75]])
    Cinv = np.linalg.inv(C)
    mu = np.array([30, 5])
    p0 = np.exp(-0.5 * np.dot(np.transpose(x - mu), np.dot(Cinv, (x - mu))))
    return p0

# Change `parameters.hpp`
params.configure()

# Set up `GridInfo`
grid = GridInfo(r, m1, m2, n1, n2, bin1, bin2, liml1, liml2)

# Set input/initial conditions
setInputGeneral(eval_p0, grid)
