"""Helper moodule for setting the initial conditions."""
import netCDF4 as nc
import numpy as np
import os
import xarray as xr

from scripts.index_functions import incrVecIndex, vecIndexToState

class GridInfo:
    """Class for storing grid parameters."""
    def __init__(self, _r: int, _m1: int, _m2: int, _n1: tuple, _n2: tuple, _binsize1: tuple, _binsize2: tuple, _liml1: tuple, _liml2: tuple):
        self.r = _r
        self.m1 = _m1
        self.m2 = _m2
        self.d = self.m1 + self.m2
        self.n1 = np.asarray(_n1, dtype="int64")
        self.n2 = np.asarray(_n2, dtype="int64")
        self.n = np.concatenate((self.n1, self.n2))
        self.binsize1 = np.asarray(_binsize1, dtype="int64")
        self.binsize2 = np.asarray(_binsize2, dtype="int64")
        self.binsize = np.concatenate((self.binsize1, self.binsize2))
        self.liml1 = np.asarray(_liml1, dtype="int64")
        self.liml2 = np.asarray(_liml2, dtype="int64")
        self.liml = np.concatenate((self.liml1, self.liml2))
        self.dx1 = np.prod(self.n1)
        self.dx2 = np.prod(self.n2)
        self.dx = np.prod(self.n)


def writeNC(x1: np.ndarray, x2: np.ndarray, s: np.ndarray) -> None:
    """
    Write X1, X2 and S to a netCDF file `input.nc`.
    NOTE: `x1` and `x2` have to be transposed (i.e. shape=(r, dx1/2)),
    as Ensign works with column-major order arrays.
    """
    if not os.path.exists("input"):
        os.makedirs("input")

    ds = xr.Dataset(
        {
            "xx1": (["n_basisfunctions", "dx1"], x1),
            "xx2": (["n_basisfunctions", "dx2"], x2),
            "ss": (["r", "r"], s)
        }
    )
    ds.to_netcdf("input/input.nc")


def setInputKD(x10: tuple, x20: tuple, grid: GridInfo) -> None:
    """
    Creates a netCDF file containing X1, X2 and S
    for an initial distribution of Kronecker delta form.
    """
    i1 = ((np.asarray(x10, dtype="int64") - grid.liml1) / grid.binsize1).astype(int)
    i2 = ((np.asarray(x20, dtype="int64") - grid.liml2) / grid.binsize2).astype(int)

    x1 = np.zeros((1, grid.dx1), dtype="float64")
    x2 = np.zeros((1, grid.dx2), dtype="float64")

    stride = 1
    idx1 = 0
    for i, el in enumerate(i1):
        idx1 += el * stride
        stride *= grid.n1[i]

    stride = 1
    idx2 = 0
    for i, el in enumerate(i2):
        idx2 += el * stride
        stride *= grid.n2[i]

    # Transpose x1 and x2, as Ensign works with column-major order arrays
    x1[0, idx1] = 1.0 / np.prod(grid.binsize1)
    x2[0, idx2] = 1.0 / np.prod(grid.binsize2)
    s = np.eye(grid.r)

    writeNC(x1, x2, s)


def setInputFunc(eval_X10: tuple[callable], eval_X20: tuple[callable], sdiag: np.ndarray, grid: GridInfo) -> None:
    """
    Creates a .netCDF file containing X1, X2 and S. The initial distribution has to be of the form 
    as P(x_1, x_2) = \sum_{i=0}^r X1_i(x_1) * S_i * X2_i(x_2). 
    The elements (callables) of the tuples `eval_X10` and `eval_X20` describe the functional dependency of X1_i(x_1) and X2_i(x_2) on the truncated state space. `sdiag` contains the elements on the diagonal of S.
    """
    r = sdiag.size

    x1 = np.zeros((r, grid.dx1), dtype="float64")
    x2 = np.zeros((r, grid.dx2), dtype="float64")
    vec_index1 = np.zeros(grid.m1)
    vec_index2 = np.zeros(grid.m2)

    for j in range(grid.dx1):
        state_x1 = vecIndexToState(vec_index1, grid.liml1, grid.binsize1)
        for i in range(r):
            x1[i, j] = eval_X10[i](state_x1 + (grid.binsize1 - 1.0) * 0.5)
        incrVecIndex(vec_index1, grid.n1, grid.m1)

    for j in range(grid.dx2):
        state_x2 = vecIndexToState(vec_index2, grid.liml2, grid.binsize2)
        for i in range(r):
            x2[i, j] = eval_X20[i](state_x2 + (grid.binsize2 - 1.0) * 0.5)
        incrVecIndex(vec_index2, grid.n2, grid.m2)

    # Normalize x1 and x2
    x1 /= np.sum(x1, axis=1)
    x2 /= np.sum(x2, axis=1)

    s = np.ones(grid.r)
    s[:r] = sdiag
    s = np.diag(s)

    writeNC(x1, x2, s)


def setInputGeneral(eval_p0: callable, grid: GridInfo) -> None:
    """
    Creates a .netCDF file containing X1, X2 and S for a general initial distribution.
    NOTE: This function creates the full P matrix,
    therefore it only should be used for small systems.
    """
    p0 = np.zeros(grid.dx)
    vec_index = np.zeros(grid.d)
    for i in range(grid.dx):
        state = vecIndexToState(vec_index, grid.liml, grid.binsize)
        p0[i] = eval_p0(state + (grid.binsize - 1.0) * 0.5)
        incrVecIndex(vec_index, grid.n, grid.d)

    # Normalize p0
    p0 = p0 / np.sum(p0)
    p0_mat = p0.reshape((grid.dx1, grid.dx2), order="F")

    # SVD of p0
    u, s, vh = np.linalg.svd(p0_mat, full_matrices=False)

    # Use only the first `r` singular values
    # Transpose x1 and x2, as Ensign works with column-major order arrays
    x1 = u[:, :grid.r].T
    s = np.diag(s[:grid.r])
    x2 = vh[:grid.r, :]

    writeNC(x1, x2, s)
