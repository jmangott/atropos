from numba import njit
import numpy as np

@njit
def CombIndexToVecIndex(comb_index, interval):
    vec_index = np.zeros(interval.size, dtype="int64")
    for i, int_ele in enumerate(interval):
        if i == interval.size - 1:
            vec_index[i] = comb_index
        else:
            vec_index[i] = comb_index % int_ele
            comb_index = comb_index // int_ele
    return vec_index

@njit
def VecIndexToCombIndex(vec_index, interval):
    stride = 1
    comb_index = 0
    for i, int_ele in enumerate(interval):
        comb_index += stride * vec_index[i]
        stride *= int_ele
    return comb_index