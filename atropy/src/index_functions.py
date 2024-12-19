"""Helper module for index calculations."""

import numpy as np
from numba import njit


@njit
def combIndexToVecIndex(comb_index: int, interval: np.ndarray) -> np.ndarray:
    """Calculates a vector index from a given combined index `comb_index`."""
    vec_index = np.zeros(interval.size, dtype="int64")
    for i, int_ele in enumerate(interval):
        if i == interval.size - 1:
            vec_index[i] = comb_index
        else:
            vec_index[i] = comb_index % int_ele
            comb_index = comb_index // int_ele
    return vec_index


@njit
def vecIndexToCombIndex(vec_index: np.ndarray, interval: np.ndarray) -> int:
    """Calculates a combined index from a given vector index `vec_index`."""
    stride = 1
    comb_index = 0
    for i, int_ele in enumerate(interval):
        comb_index += stride * vec_index[i]
        stride *= int_ele
    return comb_index


@njit
def vecIndexToState(
    vec_index: np.ndarray, liml: np.ndarray, binsize: np.ndarray
) -> np.ndarray:
    """Calculates a state from a given vector index `vec_index`."""
    state = liml + binsize * vec_index
    return state


@njit
def incrVecIndex(vec_index: np.ndarray, interval: np.ndarray, dim: int) -> None:
    """
    Increases a given vector index `vec_index`
    (i.e. for a given `vec_index`, calculate `res = vecIndexToCombIndex(vec_index) + 1`
    and subsequently `combIndexToVecIndex(res)`).
    """
    for k in range(dim):
        vec_index[k] += 1
        if vec_index[k] < interval[k]:
            return
        vec_index[k] = 0
    if dim > 0:
        vec_index[dim - 1] += 1
    return


def tensorUnfold(tensor, mode):
    """
    Cf. https://stackoverflow.com/questions/49970141/using-numpy-reshape-to-perform-3rd-rank-tensor-unfold-operation
    """
    return np.reshape(np.moveaxis(tensor, mode, 0), (tensor.shape[mode], -1), order="F")
