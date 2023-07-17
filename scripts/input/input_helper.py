"""Helper moodule for setting the initial conditions."""
import netCDF4 as nc
import numpy as np
import os

from scripts.index_functions import CombIndexToState

class grid_info:
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


def WriteNC(x1: np.ndarray, x2: np.ndarray, s: np.ndarray) -> None:
    """
    Write X1, X2 and S to a netCDF file `input.nc`.
    NOTE: `x1` and `x2` have to be transposed (i.e. shape=(r, dx1/2)),
    as Ensign works with column-major order arrays.
    """

    if not os.path.exists("input"):
        os.makedirs("input")

    ds = nc.Dataset("input/input.nc", mode="w", format="NETCDF4")
    ds.createDimension("r", np.shape(s)[0])
    ds.createDimension("n_basisfunctions", np.shape(x1)[0])
    ds.createDimension("dx1", np.shape(x1)[1])
    ds.createDimension("dx2", np.shape(x2)[1])
    xx1 = ds.createVariable("xx1", "f8", ("n_basisfunctions", "dx1"))
    xx2 = ds.createVariable("xx2", "f8", ("n_basisfunctions", "dx2"))
    ss = ds.createVariable("ss", "f8", ("r", "r"))
    xx1[:] = x1
    xx2[:] = x2
    ss[:] = s
    ds.close()


def SetInputKD(x10: tuple, x20: tuple, grid: grid_info) -> None:
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

    WriteNC(x1, x2, s)


def SetInputGeneral(eval_p0: callable, grid: grid_info) -> None:
    """
    Creates a .netCDF file containing X1, X2 and S for a general initial distribution.
    NOTE: This function creates the full P matrix,
    therefore it only should be used for small systems.
    """

    p0 = np.zeros((grid.dx1, grid.dx2))
    for i in range(grid.dx1):
        for j in range(grid.dx2):
            comb_index = i + grid.dx1 * j
            state = CombIndexToState(comb_index, grid.n, grid.liml, grid.binsize)
            p0[i, j] = eval_p0(state + (grid.binsize - 1.0) * 0.5)

    # Normalize p0
    p0 = p0 / np.sum(p0)

    # SVD of p0
    u, s, vh = np.linalg.svd(p0, full_matrices=False)

    # Use only the first `r` singular values
    # Transpose x1 and x2, as Ensign works with column-major order arrays
    x1 = u[:, :grid.r].T
    s = np.diag(s[:grid.r])
    x2 = vh[:grid.r, :]

    WriteNC(x1, x2, s)
