"""
Script for setting parameters (`parameters.hpp`) and initial conditions (`input.nc`) 
for the BAX pore assembly model.
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
m1 = params.kM1 = 6
m2 = params.kM2 = 5
n1 = params.kN1 = (16, 6, 6, 6, 106, 106)
n2 = params.kN2 = (86, 21, 21, 21, 21)
bin1 = params.kBinsize1 = (1, 1, 1, 1, 1, 1)
bin2 = params.kBinsize2 = (1, 1, 1, 1, 1)
liml1 = params.kLiml1 = (0, 0, 0, 0, 0, 0)
liml2 = params.kLiml2 = (0, 0, 0, 0, 0)

C = 0.2
Cinv = 1 / C
mu = np.array([0, 0, 0, 0, 100, 0, 80, 0, 0, 0, 0])

def eval_x(x: np.ndarray, mu: np.ndarray) -> float:
    return np.exp(-0.5 * Cinv * np.dot(np.transpose(x - mu), (x - mu)))

# Change `parameters.hpp`
params.configure()

# Set up `GridInfo`
grid = GridInfo(r, m1, m2, n1, n2, bin1, bin2, liml1, liml2)

# Set input/initial conditions
setInputFunc([lambda x: eval_x(x, mu[:m1])], [lambda x: eval_x(x, mu[m1:])], np.array([1]), grid)
