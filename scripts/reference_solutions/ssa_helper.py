"""Helper module for solving the CME with various SSA implementations."""
import gillespy2
from numba import njit, int64, float64
# from numba.experimental import jitclass
import numpy as np
import pysb.core
from pysb.simulator import StochKitSimulator

from scripts.index_functions import vecIndexToCombIndex, incrVecIndex


@njit
def _calculateObservables(ssa_result: np.ndarray, n_time: int, n_runs: int, m: int, n: np.ndarray, n_min: np.ndarray, slice_vec: np.ndarray, idx_2D: np.ndarray):
    dx_tot2 = np.prod(n[idx_2D])
    P_marginal = [[np.zeros(n_el, dtype="float64") for n_el in n] for _ in range(n_time)]
    P_marginal2D = np.zeros((n_time, dx_tot2), dtype="float64")
    P_sliced = [[np.zeros(n_el, dtype="float64") for n_el in n] for _ in range(n_time)]
    P_sliced2D = np.zeros((n_time, dx_tot2), dtype="float64")

    for j in range(n_time):
        lin_dset_marginal2D = np.zeros(n_runs, dtype="int64")
        lin_dset_sliced = np.zeros(n_runs, dtype="int64")

        for i in range(n_runs):
            # Linearize the population numbers for all times and runs
            lin_dset_marginal2D[i] = vecIndexToCombIndex(ssa_result[i, j, idx_2D] - n_min[idx_2D], n[idx_2D])
            lin_dset_sliced[i] = vecIndexToCombIndex(ssa_result[i, j, :] - n_min, n)

        for k in range(m):
            P_marginal[j][k] = np.bincount(
                ssa_result[:, j, k] - n_min[k], minlength=n[k], weights=np.ones(n_runs, dtype="float64"))
            P_marginal[j][k] /= np.sum(P_marginal[j][k])

            # Leave `slice_vec` untouched, store the normalization in `comp_vec`
            comp_vec = np.copy(slice_vec) - n_min

            for l in range(n[k]):
                # Calculate `P_sliced` for species `l`,
                # modify `comp_vec` for obtaining the correct slice and calculate a combined index
                comp_vec[k] = l
                comp_index = vecIndexToCombIndex(comp_vec, n)
                P_sliced[j][k][l] = np.sum(
                    np.equal(lin_dset_sliced, comp_index)) / n_runs
                # P_sliced has to be divided by the total number of runs (it is not normalized to 1)

        comp_vec2D = np.copy(slice_vec) - n_min
        vec_index2D = np.zeros(2)
        for k in range(dx_tot2):
            comp_vec2D[idx_2D] = vec_index2D
            comp_index2D = vecIndexToCombIndex(comp_vec2D, n)
            P_sliced2D[j][k] = np.sum(
                np.equal(lin_dset_sliced, comp_index2D)) / n_runs
            incrVecIndex(vec_index2D, n[idx_2D], 2)

        # For a given time, `P_marginal2D` is stored as a vector (in the row-major convention)
        P_marginal2D[j] = np.bincount(
            lin_dset_marginal2D, minlength=dx_tot2).astype(np.float64)
        P_marginal2D[j] /= np.sum(P_marginal2D[j])

    return P_marginal, P_marginal2D, P_sliced, P_sliced2D


@njit
def _calculateFullDistribution(ssa_result: np.ndarray, n_time: int, n_runs: int, n: np.ndarray, n_min: np.ndarray):
    dx_tot = np.prod(n)
    P_full = np.zeros((n_time, dx_tot), dtype="float64")

    for j in range(n_time):
        lin_dset_sliced = np.zeros(n_runs, dtype="int64")

        for i in range(n_runs):
            # Linearize the population numbers for all times and runs
            lin_dset_sliced[i] = vecIndexToCombIndex(ssa_result[i, j, :] - n_min, n)

        P_full[j] = np.bincount(lin_dset_sliced, minlength=dx_tot, weights=np.ones(n_runs, dtype="float64"))
        P_full[j] /= n_runs

    return P_full


class SSASol:
    """Helper class for calculating marginal and sliced distributions and the full probability distribution."""
    def __init__(self, _ssa_result: np.ndarray):
        self.ssa_result = _ssa_result
        self.n_runs = self.ssa_result.shape[0]
        self.n_time = self.ssa_result.shape[1]
        self.m = self.ssa_result.shape[2]

        self.n_max = np.zeros(self.m, dtype="int64")
        self.n_min = np.zeros(self.m, dtype="int64")
        for i in range(self.m):
            self.n_max[i] = np.amax(self.ssa_result[:, :, i])
            self.n_min[i] = np.amin(self.ssa_result[:, :, i])
        self.n = self.n_max - self.n_min + 1

    def calculateObservables(self, slice_vec: np.ndarray, idx_2D: np.ndarray):
        """Calculate marginal and sliced distributions."""
        return _calculateObservables(self.ssa_result, self.n_time, self.n_runs, self.m, self.n, self.n_min, slice_vec, idx_2D)

    def calculateFullDistribution(self):
        """
        Calculates the full probability distribution.
        NOTE: This function should be used only for small systems.
        """
        return _calculateFullDistribution(self.ssa_result, self.n_time, self.n_runs, self.n, self.n_min)


def calculateNRuns(eval_P0: callable, sweeps: int, interval: np.ndarray, liml: np.ndarray) -> tuple[np.ndarray, int]:
    """
    Helper function for determining how many runs with a certain initial state have to be performed (`n_runs`) to sample the given initial distribution `eval_P0`.
    """
    dx = np.prod(interval)
    m = interval.size
    vec_index = np.zeros(m)

    n_runs = np.empty(dx, dtype="int64")
    n_runs_tot = 0
    P_lin = np.empty(dx)
    # Evaluate the probability function for all x in the sample space
    for i in range(dx):
        P_lin[i] = eval_P0(vec_index + liml)
        # Multiply the probability by `sweeps` and round it to the nearest integer
        n_runs[i] = np.around(P_lin[i] * sweeps)
        n_runs_tot += n_runs[i]
        incrVecIndex(vec_index, interval, m)

    err = np.linalg.norm(n_runs / n_runs_tot - P_lin)
    print("n_runs_tot:", n_runs_tot)
    print("error:", err)
    return n_runs, n_runs_tot


def runStochkitSSA(n_runs: np.ndarray, n_runs_tot: int, tspan: np.ndarray, interval: np.ndarray, liml: np.ndarray, model: pysb.core.Model, observables: list) -> np.ndarray:
    """
    Run the Stochkit2 SSA implementation with starting values in accordance with the initial distribution described by `eval_P0`. `sweeps` is both an approximation and a control parameter for the total number of SSA sweeps. The actual number of total sweeps is given by `n_runs_tot`. `interval` and `liml` determine the sample space, which should cover most of the initial distribution.
    
    Output: numpy.ndarray with shape `(n_runs_tot, n_time, interval.size)`, which collects all trajectories.
    """
    dx = np.prod(interval)
    m = interval.size
    i_runs = 0
    n_time = tspan.size

    result = np.empty((n_runs_tot, n_time, m), dtype="int64")
    vec_index = np.zeros(m)

    # Loop through the sample space and perform `n_runs[i]` SSA sweeps for state `i`
    for i in range(dx):
        if n_runs[i] > 0.0:
            # Modify the initial conditions
            state = vec_index + liml
            for j, ic in enumerate(model.initials):
                parameter_name = ic.value.name
                model.parameters[parameter_name].value = state[j]

            sim = StochKitSimulator(model, tspan=tspan)
            simulation_result = sim.run(n_runs=n_runs[i])

            # Convert the result into a numpy array
            for j_runs in range(n_runs[i]):
                for i_obs, obs in enumerate(observables):
                    result[i_runs + j_runs, :, i_obs] = simulation_result.observables[j_runs][obs]
            
            i_runs += n_runs[i]
        incrVecIndex(vec_index, interval, m)
    return result


def runGillespySSA(n_runs: np.ndarray, n_runs_tot: int, interval: np.ndarray, liml: np.ndarray, model: gillespy2.core.model.Model, observables: list) -> np.ndarray:
    """
    Run the Gillespy2 SSA implementation with starting values in accordance with the initial distribution described by `eval_P0`. `sweeps` is both an approximation and a control parameter for the total number of SSA sweeps. The actual number of total sweeps is given by `n_runs_tot`. `interval` and `liml` determine the sample space, which should cover most of the initial distribution.
    
    Output: numpy.ndarray with shape `(n_runs_tot, n_time, interval.size)`, which collects all trajectories.
    """
    dx = np.prod(interval)
    m = interval.size
    i_runs = 0
    n_time = len(model.tspan)

    result = np.empty((n_runs_tot, n_time, m), dtype="int64")
    vec_index = np.zeros(m)

    # Loop through the sample space and perform `n_runs[i]` SSA sweeps for state `i`
    for i in range(dx):
        if n_runs[i] > 0.0:
            # Modify the initial conditions
            state = vec_index + liml
            for j, species in enumerate(model.get_all_species().keys()):
                model.get_species(species).initial_value=state[j]

            simulation_result = model.run(number_of_trajectories=n_runs[i], algorithm="SSA")

            # Convert the result into a numpy array
            for j_runs, trajectory in enumerate(simulation_result):
                for i_obs, obs in enumerate(observables):
                    result[i_runs + j_runs, :, i_obs] = trajectory[obs]
            
            i_runs += n_runs[i]
        incrVecIndex(vec_index, interval, m)
    return result