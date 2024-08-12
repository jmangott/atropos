"""Script for generating the SSA solution for a given model."""
import argparse
import numpy as np
from scipy.special import factorial
import sys
import time

from ssa_helper import *

parser = argparse.ArgumentParser(
                    prog='ssa',
                    usage='python3 scripts/reference_solutions/ssa.py --model lambda_phage --nruns 10000',
                    description='This script generates the SSA solution for a given model.')

parser.add_argument('-m', 
                    '--model', 
                    type=str, 
                    required=True, 
                    help='Specify a model ("toggle_switch", "lambda_phage", "diffusive_toggle_switch", "enzymatic_futile_cycle", "cascade")',
                    )
parser.add_argument('-n', 
                    '--n_runs', 
                    type=int, 
                    required=True, 
                    help="Specify the number of runs/trajectories",
                    )
args = parser.parse_args()

model = args.model
sweeps = args.n_runs

if model == "toggle_switch":
    import scripts.models.toggle_switch_gillespy as gillespy_model

    liml = np.array([0, 0])
    n = np.array([51, 51])

    mu = np.array([30.0, 5.0])
    C = 0.2
    Cinv = 1 / C

    sum_vec = np.zeros(n.size)
    for i, n_el in enumerate(n):
        for x in range(n_el):
            x_vec = np.array([x])
            sum_vec[i] += np.exp(-0.5 * Cinv * np.linalg.norm(x + liml[i] - mu[i]) ** 2)
    gamma = 1 / np.prod(sum_vec)

    def eval_P0(x: np.ndarray) -> float:
        return gamma * np.exp(-0.5 * Cinv * np.linalg.norm(x - mu) ** 2)

elif model == "lambda_phage":
    import scripts.models.lambda_phage_gillespy as gillespy_model

    liml = np.array([0, 0, 0, 0, 0])
    n = np.array([16, 41, 11, 11, 11])

    def eval_P0(x: np.ndarray) -> float:
        abs_x = np.sum(x)
        if (abs_x <= 3):
            P0 = factorial(3) * (0.05 ** abs_x) * \
                ((1.0 - 5 * 0.05) ** (3 - abs_x)) / \
                (np.prod(factorial(x)) * factorial(3 - abs_x))
        else:
            P0 = 0.0
        return P0

elif model == "diffusive_toggle_switch":
    import scripts.models.diffusive_toggle_switch_gillespy as gillespy_model

    liml = np.array([28, 3, 0, 0, 0, 0, 0, 0])
    n = np.array([5, 5, 4, 4, 4, 4, 4, 4])

    mu = np.array([30.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    C = 0.2
    Cinv = 1 / C

    sum_vec = np.zeros(n.size)
    for i, n_el in enumerate(n):
        for x in range(n_el):
            x_vec = np.array([x])
            sum_vec[i] += np.exp(-0.5 * Cinv * np.linalg.norm(x + liml[i] - mu[i]) ** 2)
    gamma = 1 / np.prod(sum_vec)

    def eval_P0(x: np.ndarray) -> float:
        return gamma * np.exp(-0.5 * Cinv * np.linalg.norm(x - mu) ** 2)

elif model == "enzymatic_futile_cycle":
    import scripts.models.enzymatic_futile_cycle_gillespy as gillespy_model

    liml = np.array([30, 90, 2, 2, 0, 0])
    n = np.ones(liml.size, dtype="int")

    def eval_P0(x: np.ndarray) -> float:
        return 1.0

elif model == "cascade":
    import scripts.models.cascade_gillespy as gillespy_model
    liml = np.zeros(20)
    n = np.ones(liml.size, dtype="int")

    def eval_P0(x: np.ndarray) -> float:
        return 1.0

else:
    print(parser.prog+":", 'error: the following arguments for `model` are allowed: "toggle_switch", "lambda_phage", "diffusive_toggle_switch", "enzymatic_futile_cycle" or "cascade"')
    sys.exit(1)

m = n.size # len(pysb_model.model.observables)
observables = [obs for obs in gillespy_model.model.get_all_species().keys()]

n_runs, n_runs_tot = calculateNRuns(eval_P0, sweeps, n, liml)

t_start = time.time_ns()
result = runGillespySSA(n_runs, n_runs_tot, n, liml, gillespy_model.model, observables)
t_stop = time.time_ns()
wall_time = (t_stop - t_start) * 1e-9
print("Time elapsed:", str(wall_time)+"s")
with open("scripts/reference_solutions/" + model + "_ssa_{:.0e}".format(sweeps) + ".npz", "wb") as f:
    np.savez(f, result=result, wall_time=wall_time)