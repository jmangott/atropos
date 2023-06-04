import numpy as np
# import pysb_model_ts as pysb_model
import pysb_model_lp as pysb_model
# import pysb_model_bax_seq as pysb_model
from pysb.simulator import StochKitSimulator
import time
from scripts.index_functions import VecIndexToCombIndex

n_time = 100 + 1
n_runs = 100000
tspan = np.linspace(0.0, 100.0, n_time)
fname = "bax_stochkit_result"

d = len(pysb_model.model.observables)
slice_vec = [int(pysb_model.model.initials[i].value.value) for i in range(d)]

# Run StochKit simulation
start_time = time.perf_counter()
sim = StochKitSimulator(pysb_model.model, tspan=tspan)
simulation_result = sim.run(n_runs=n_runs)
end_time = time.perf_counter()
print("CPU time:", end_time - start_time, "seconds")

# Convert the result into a numpy array
result = np.zeros((n_runs, n_time, d), dtype="int64")
for i_runs in range(n_runs):
    for i_obs, obs in enumerate(pysb_model.model.observables):
        result[i_runs, :, i_obs] = simulation_result.observables[i_runs][obs.name]

# Calculate lower and upper population bounds
n_max = np.zeros(d, dtype="int64")
n_min = np.zeros(d, dtype="int64")
n = np.zeros(d, dtype="int64")
for i in range(d):
    n_max[i] = np.amax(result[:, :, i])
    n_min[i] = np.amin(result[:, :, i])
n = n_max - n_min + 1
dx_tot = np.prod(n[0:2])

P_marginal = [np.zeros((n_el, n_time), dtype="float64") for n_el in n]
P_marginal2D = [np.zeros(dx_tot, dtype="float64") for t in range(n_time)]
P_marginal2D_mat = [np.zeros((n[0], n[1]), dtype="float64") for t in range(n_time)]
P_sliced = [np.zeros((n_el, n_time), dtype="float64") for n_el in n]

# Calculate marginal probability distributions
for k in range(d):
    for j in range(n_time):
        lin_dset = np.zeros((n_runs), dtype="int64")
        for i in range(n_runs):
            lin_dset[i] = result[i, j, k] - n_min[k]
        P_marginal[k][:, j] = np.bincount(lin_dset, minlength=n[k])

# Normalize P_marginal for all sampling times
for i in range(d):
    for j in range(n_time):
        P_marginal[i][:, j] /= np.sum(P_marginal[i][:, j])

for j in range(n_time):
    # Linearize the population numbers (i.e. a unique number is assigned to every configuration),
    # then count the linearized numbers with np.bincount for all sampling times
    lin_dset = np.zeros((n_runs), dtype="int64")
    for i in range(n_runs):
        lin_dset[i] = VecIndexToCombIndex(result[i, j, 0:2] - n_min[0:2], n[0:2])
    P_marginal2D[j] = np.bincount(lin_dset, minlength=dx_tot).astype(np.float64)

# Normalize P for all sampling times
for i in range(n_time):
    P_marginal2D[i] /= np.sum(P_marginal2D[i])
    P_marginal2D_mat[i] = np.reshape(P_marginal2D[i], n[0:2], order='F').T

# Calculate sliced probability distributions
for j in range(n_time):
    lin_dset = np.zeros((n_runs), dtype="int64")
    for i in range(n_runs):
        lin_dset[i] = VecIndexToCombIndex(result[i, j, :] - n_min, n)
    for k in range(d):
        comp_vec = slice_vec.copy() - n_min
        for l in range(n[k]):
            comp_vec[k] = l
            comp_index = VecIndexToCombIndex(comp_vec, n)
            P_sliced[k][l, j] = np.sum(np.equal(lin_dset, comp_index))

# Normalize P_sliced
for i in range(d):
    P_sliced[i] /= n_runs

# Save result
with open("scripts/reference_solutions/" + fname + ".npy", "wb") as f:
    np.save(f, P_marginal)
    np.save(f, P_sliced)
    # np.save(f, P_marginal2D_mat)
    np.save(f, n)
    np.save(f, n_min)
    np.save(f, n_max)