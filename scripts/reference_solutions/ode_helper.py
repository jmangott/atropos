"""Helper module for solving the full CME."""
from numba import njit
import numpy as np
import matplotlib.pyplot as plt

from scripts.index_functions import IncrVecIndex

@njit
def CalculateShift(nu: np.ndarray, interval: np.ndarray) -> float:
    """Calculate linearized shift for a given stoichiometric vector `nu`."""
    shift = 0
    stride = 1
    m = interval.size
    for i in range(m):
        shift += nu[i] * stride
        stride *= interval[i]
    return shift


@njit
def ShiftArray(input_array: np.ndarray, nu: np.ndarray, interval: np.ndarray) -> np.ndarray:
    """Calculate the shifted probability distribution in the CME for a given stoichiometric vector `nu`."""
    output_array = np.copy(input_array)
    shift = CalculateShift(nu, interval)
    n_rows = output_array.size
    m = interval.size
    vec_index = np.zeros(m, dtype="int64")
    for i in range(n_rows):
        if (shift < 0 and i - shift < n_rows) or (shift >= 0 and i - shift >= 0):
            output_array[i] = input_array[i - shift]
        else:
            output_array[i] = 0.0
            IncrVecIndex(vec_index, interval, m)
            continue
        for k in range(m):
            if (nu[k] > 0 and vec_index[k] - nu[k] < 0) or (nu[k] < 0 and vec_index[k] - nu[k] >= interval[k]):
                output_array[i] = 0.0
                break
        IncrVecIndex(vec_index, interval, m)
    return output_array


@njit
def EvaluateProp(prop_fun: callable, nu: np.ndarray, interval: np.ndarray) -> np.ndarray:
    """Evaluate a given propensity function `prop_fun`."""
    dx = np.prod(interval)
    output_array = np.zeros(dx)
    m = interval.size
    vec_index = np.zeros(m, dtype="int64")
    for i in range(dx):
        output_array[i] = prop_fun(vec_index - nu)
        IncrVecIndex(vec_index, interval, m)
    return output_array


def ConstructP0(eval_P0: callable, interval: np.ndarray) -> np.ndarray:
    """Set up the initial probability distribution according to a given function `eval_P0`."""
    dx = np.prod(interval)
    m = interval.size
    P0 = np.zeros(dx)
    vec_index = np.zeros(m)
    for i in range(dx):
        P0[i] = eval_P0(vec_index)
        IncrVecIndex(vec_index, interval, m)
    return P0 / np.sum(P0)


# @njit
def CalculateObservables(y: np.ndarray, interval: np.ndarray, r: int, i2D: np.ndarray, m1: int, slice_vec: np.ndarray):
    """Calculate marginal and sliced distributions and the best approximation."""
    dx = np.prod(interval)
    m = interval.size

    P_full = y
    P_marginal = [[np.zeros(n_el) for n_el in interval] for _ in range(y.shape[0])]
    P_marginal2D = [np.zeros((interval[0], interval[1])) for _ in range(y.shape[0])]
    P_sliced = [[np.zeros(n_el) for n_el in interval] for _ in range(y.shape[0])]
    P_sliced2D = [np.zeros((interval[0], interval[1])) for _ in range(y.shape[0])]
    P_best_approximation = np.zeros(y.shape)

    vec_index_c = np.zeros(m - 2, dtype="int64")
    vec_index_k = np.zeros(m - 1, dtype="int64")
    slice_vec_c = np.zeros(m - 2, dtype="int64")
    slice_vec_k = np.zeros(m - 1, dtype="int64")

    slice_vec_c[:i2D[0]] = slice_vec[:i2D[0]]
    slice_vec_c[i2D[0]:i2D[1]-1] = slice_vec[i2D[0]+1:i2D[1]]
    slice_vec_c[i2D[1]-1:] = slice_vec[i2D[1]+1:]

    for i in range(y.shape[0]):
        vec_index = np.zeros(m, dtype="int64")
        for j in range(dx):
            P_marginal2D[i][vec_index[i2D[0]], vec_index[i2D[1]]] += y[i, j]
            vec_index_c[:i2D[0]] = vec_index[:i2D[0]]
            vec_index_c[i2D[0]:i2D[1]-1] = vec_index[i2D[0]+1:i2D[1]]
            vec_index_c[i2D[1]-1:] = vec_index[i2D[1]+1:]

            if np.all(vec_index_c == slice_vec_c):
                P_sliced2D[i][vec_index[i2D[0]], vec_index[i2D[1]]] = y[i, j]
            for k in range(m):
                P_marginal[i][k][vec_index[k]] += y[i, j]
                vec_index_k[:k] = vec_index[:k]
                vec_index_k[k:] = vec_index[k+1:]
                slice_vec_k[:k] = slice_vec[:k]
                slice_vec_k[k:] = slice_vec[k+1:]
                if np.all(vec_index_k == slice_vec_k):
                    P_sliced[i][k][vec_index[k]] = y[i, j]

            IncrVecIndex(vec_index, interval, m)

        P = y[i, :].reshape((np.prod(interval[m1:]), np.prod(interval[:m1])))
        u, s, vh = np.linalg.svd(P, full_matrices=False)
        # Use only the first `r` singular values
        X1 = u[:, :r]
        S = s[:r]
        X2h = vh[:r, :]
        P_best_approximation[i, :] = ((X1 * S) @ X2h).flatten()

    return P_full, P_marginal, P_marginal2D, P_sliced, P_sliced2D, P_best_approximation