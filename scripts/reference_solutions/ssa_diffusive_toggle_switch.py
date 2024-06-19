import numpy as np
import matplotlib.pyplot as plt
import scripts.models.diffusive_toggle_switch_pysb as pysb_model
from pysb.simulator import StochKitSimulator

from ssa_helper import *

n_time = 1 + 1
tspan = np.linspace(0.0, 500, n_time)

# This defines the small sample space
liml = np.array([28, 3, 0, 0, 0, 0, 0, 0])
n = np.array([5, 5, 4, 4, 4, 4, 4, 4])

sweeps = 1e5
meta_runs = 1

m = n.size # len(pysb_model.model.observables)
observables = [obs.name for obs in pysb_model.model.observables]

mu = np.array([30.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
C = 0.2
Cinv = 1 / C

# Calculate the normalization
sum_vec = np.zeros(n.size)
for i, n_el in enumerate(n):
    for x in range(n_el):
        x_vec = np.array([x])
        sum_vec[i] += np.exp(-0.5 * Cinv * np.linalg.norm(x + liml[i] - mu[i]) ** 2)
gamma = 1 / np.prod(sum_vec)

@njit
def eval_P0(x: np.ndarray) -> float:
    return gamma * np.exp(-0.5 * Cinv * np.linalg.norm(x - mu) ** 2)

n_runs, n_runs_tot = calculateNRuns(eval_P0, sweeps, n, liml)
result_tot = np.empty((meta_runs * n_runs_tot, n_time, m), dtype="int64")

for i in range(meta_runs):
    result_tot[i * n_runs_tot : (i + 1) * n_runs_tot] = runSSA(n_runs, n_runs_tot, tspan, n, liml, pysb_model.model, observables)
    print(i)

np.save("scripts/reference_solutions/dts_ssa_1e5.npy", result_tot)