"""Helper module for postprocessing the DLR results."""
import matplotlib as mpl
mpl.use('pgf')
import matplotlib.pyplot as plt
# import netCDF4 as nc
from numba import njit
import numpy as np
import os
import xarray as xr

if not os.path.exists("plots"):
    os.makedirs("plots")

from scripts.index_functions import incrVecIndex, vecIndexToCombIndex

# Update Matplotlib settings
plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    "figure.figsize": (6, 4),
    "lines.linewidth": "1",
    "lines.markersize": "4.5",
    "lines.linestyle": "none",
    "axes.titlesize": "medium",
    "axes.prop_cycle": mpl.cycler(color=["r", "b", "k"], marker=["o", "x", "."])
})


class GridInfo:
    """Class for storing DLR parameters."""
    def __init__(self, _ds: xr.core.dataset.Dataset):
        self.n1 = _ds["n1"].values.astype(int)
        self.n2 = _ds["n2"].values.astype(int)
        self.dx1 = np.prod(self.n1)
        self.dx2 = np.prod(self.n2)
        self.n = np.concatenate((self.n1, self.n2))
        self.m1 = _ds.dims["m1"]
        self.m2 = _ds.dims["m2"]
        self.r = _ds.dims["r"]
        self.d = self.m1 + self.m2
        self.bin = _ds["binsize"].values
        self.liml = _ds["liml"].values.astype(int)
        self.t = _ds["t"].values
        self.dt = _ds["dt"].values


@njit
def calculateXSum(input_array: np.ndarray, n: np.ndarray) -> list[np.ndarray]:
    """Integrates `input_array` (X1 or X2) for each species over all the remaining species."""
    r = input_array.shape[0]
    m = n.size
    dx = np.prod(n)
    vec_index = np.zeros(m, dtype=np.int64)
    output_array = [np.zeros((n_el, r), dtype="float64") for n_el in n]
    for i in range(dx):
        for k_vec, vec_el in enumerate(vec_index):
            output_array[k_vec][vec_el, :] += input_array[:, i]
        incrVecIndex(vec_index, n, m)
    return output_array


# @njit
def calculateXSum2D(input_array: np.ndarray, n: np.ndarray, idx_2D: np.ndarray):
    """Integrates `input_array` (X1 or X2) for two species (given by `idx_2D`) over all the remaining species."""
    r = input_array.shape[0]
    m = n.size
    dx = np.prod(n)
    dx2 = np.prod(n[idx_2D])
    vec_index = np.zeros(m, dtype="int64")
    output_array = np.zeros((dx2, r))
    for i in range(dx):
        comb_index = vecIndexToCombIndex(vec_index[idx_2D], n[idx_2D])
        output_array[comb_index, :] += input_array[:, i]
        incrVecIndex(vec_index, n, m)
    return output_array


@njit
def calculateXSlice(input_array: np.ndarray, n: np.ndarray, slice_vec: np.ndarray) -> list[np.ndarray]:
    """Evaluates `input_array` (X1 or X2) for each species at `slice_vec` for all the remaining species."""
    r = input_array.shape[0]
    output_array = [np.zeros((n_el, r), dtype="float64") for n_el in n]
    for k, n_el in enumerate(n):
        vec_index = np.copy(slice_vec)
        for i in range(n_el):
            vec_index[k] = i
            comb_index = vecIndexToCombIndex(vec_index, n)
            output_array[k][i, :] = input_array[:, comb_index]
    return output_array


# @njit
def calculateXSlice2D(input_array: np.ndarray, n: np.ndarray, slice_vec: np.ndarray, idx_2D: np.ndarray):
    """Evaluates `input_array` (X1 or X2) for two species (given by `idx_2D`) at `slice_vec` for all the remaining species."""
    r = input_array.shape[0]
    dx2 = np.prod(n[idx_2D])
    output_array = np.zeros((dx2, r))
    comp_vec = np.copy(slice_vec)
    vec_index = np.zeros(2, dtype="int64")
    for i in range(dx2):
        comp_vec[idx_2D] = vec_index
        comb_index = vecIndexToCombIndex(comp_vec, n)
        output_array[i, :] += input_array[:, comb_index]
        incrVecIndex(vec_index, n[idx_2D], 2)
    return output_array


class LRSol:
    def __init__(self, _ds: xr.core.dataset.Dataset, _grid: GridInfo):
        """
        Helper class for calculating marginal and sliced distributions and the full probability distribution.
        """
        self.grid = _grid
        self.X1 = _ds['X'].values
        self.S = _ds['S'].values
        self.X2 = _ds['V'].values


    def fullDistribution(self) -> np.ndarray:
        """
        Calculates the full probability distribution.
        NOTE: This function should be used only for small systems.
        """
        P_full = np.matmul(np.transpose(self.X1), np.matmul(self.S, self.X2))
        
        return P_full.flatten(order='F')


    def marginalDistributions(self) -> list[np.ndarray]:
        """Calculates marginal distributions."""
        P_marginal = [np.zeros(n_el) for n_el in self.grid.n]

        X1_sum_tot = np.sum(self.X1, axis=1)
        X2_sum_tot = np.sum(self.X2, axis=1)

        X1_sum = calculateXSum(self.X1, self.grid.n1)
        X2_sum = calculateXSum(self.X2, self.grid.n2)

        for i in range(self.grid.m1):
            P_marginal[i] = np.matmul(X1_sum[i], np.matmul(self.S, np.transpose(X2_sum_tot)))

        for i in range(self.grid.m2):
            P_marginal[self.grid.m1 + i] = np.matmul(X1_sum_tot, np.matmul(self.S, np.transpose(X2_sum[i])))
        
        return P_marginal


    def marginalDistribution2D(self, idx_2D: np.ndarray) -> np.ndarray:
        """Calculates a 2D marginal distribution."""
        if  (idx_2D[0] < self.grid.m1) and (idx_2D[1] < self.grid.m1):
            X1_sum = calculateXSum2D(self.X1, self.grid.n1, idx_2D)
            X2_sum= np.sum(self.X2, axis=1)

        elif (idx_2D[0] >= self.grid.m1) and (idx_2D[1] >= self.grid.m1):
            X1_sum = np.sum(self.X1, axis=1)
            X2_sum = calculateXSum2D(self.X2, self.grid.n2, idx_2D - self.grid.m1)

        else:
            X1_sum = calculateXSum(self.X1, self.grid.n1)[idx_2D[0]]
            X2_sum = calculateXSum(self.X2, self.grid.n2)[idx_2D[1] - self.grid.m1]

        P_marginal2D_vec = np.matmul(X2_sum, np.matmul(self.S, X1_sum))
        P_marginal2D = np.reshape(P_marginal2D_vec, self.grid.n[idx_2D], order="F")

        return P_marginal2D


    def slicedDistributions(self, slice_vec: np.ndarray) -> list[np.ndarray]:
        """Calculates sliced distributions for a given `slice_vec`."""
        P_sliced = [np.zeros(n_el) for n_el in self.grid.n]

        slice_vec_index = ((np.asarray(slice_vec, dtype="int64") - self.grid.liml) / self.grid.bin).astype(int)
        X1_slice = calculateXSlice(self.X1, self.grid.n1, slice_vec_index[:self.grid.m1])
        X2_slice = calculateXSlice(self.X2, self.grid.n2, slice_vec_index[self.grid.m1:])

        X1_slice_tot = X1_slice[0][slice_vec_index[0], :]
        X2_slice_tot = X2_slice[0][slice_vec_index[self.grid.m1], :]

        for i in range(self.grid.m1):
            P_sliced[i] = np.matmul(X1_slice[i], np.matmul(self.S, np.transpose(X2_slice_tot)))

        for i in range(self.grid.m2):
            P_sliced[self.grid.m1 + i] = np.matmul(X1_slice_tot, np.matmul(self.S, np.transpose(X2_slice[i])))

        return P_sliced


    def slicedDistribution2D(self, slice_vec: np.ndarray, idx_2D: np.ndarray) -> np.ndarray:
        """
        Calculates a 2D sliced distribution for a given `slice_vec`.
        NOTE: This is 'hard-coded' for the first partition and assumes that only two species lie in that partition.
        """
        slice_vec_index = ((np.asarray(slice_vec, dtype="int64") - self.grid.liml) / self.grid.bin).astype(int)

        if (idx_2D[0] < self.grid.m1) and (idx_2D[1] < self.grid.m1):
            X1_slice = calculateXSlice2D(self.X1, self.grid.n1, slice_vec_index[:self.grid.m1], idx_2D)
            comb_index = vecIndexToCombIndex(slice_vec_index[self.grid.m1:], self.grid.n2)
            X2_slice = self.X2[:, comb_index]

        elif (idx_2D[0] >= self.grid.m1) and (idx_2D[1] >= self.grid.m1):
            comb_index = vecIndexToCombIndex(slice_vec_index[:self.grid.m1], self.grid.n1)
            X1_slice = self.X1[:, comb_index]
            X2_slice = calculateXSlice2D(self.X2, self.grid.n2, slice_vec_index[self.grid.m1:], idx_2D - self.grid.m1)

        else:
            X1_slice = calculateXSlice(self.X1, self.grid.n1, slice_vec_index[:self.grid.m1])[idx_2D[0]]
            X2_slice = calculateXSlice(self.X2, self.grid.n2, slice_vec_index[self.grid.m1:])[idx_2D[1] - self.grid.m1]

        P_sliced2D_vec = np.matmul(X2_slice, np.matmul(self.S, X1_slice))
        P_sliced2D = np.reshape(P_sliced2D_vec, self.grid.n[idx_2D], order="F")
        return P_sliced2D


def plotP2D(axs, P: list[np.ndarray], mesh: tuple, title: list[str]):
    """Plots multiple datasets stored in `P` on a common `mesh` as 2D contour plots."""
    levels = np.linspace(np.amin([np.amin(P_el) for P_el in P]), np.amax([np.amax(P_el) for P_el in P]), 9)
    for i, ax in enumerate(axs.flatten()):
        if i < len(P):
            ax.contour(mesh[0], mesh[1], P[i], levels=levels)
            ax.set_title(title[i])


def plotP1D(ax, P: list, mesh: tuple, label: list, **kwargs) -> tuple[plt.Figure, np.ndarray]:
    """Helper function for `plotP1Dmult`."""
    plt.setp(ax, **kwargs)
    for i, P_el in enumerate(P):
        ax.plot(mesh, P_el, label=label[i])


def plotP1Dmult(axs, P: list[np.ndarray], mesh: list, label: list, idx: list, Plabel: str, *, stats: bool = False):
    """
    Plots a subset given by `idx` of all datasets stored in `P` on a common `mesh`. If a dataset is defined on a larger mesh, then the outlying points will be truncated. If the dataset is defined on a smaller mesh, all additional points will be set to 0. When `stats` is set to `True`, then the last dataset `P[-1]` is treated as a reference solution and the maximal error for every other dataset will be calculated with respect to `P[-1]`.
    """
    for i, ax in enumerate(axs.flatten()):
        if i < len(idx):
            j = idx[i]
            xlabel = "$x_{{" + str(j + 1) + "}}$"
            ylabel = "$P_{{\mathrm{{" + Plabel + "}}}}(x_{{" + str(j + 1) + "}})$"

            # Extend mesh[k] so that mesh[0] and mesh[k] cover the same interval (k > 0)
            for k in range(1, len(mesh)):
                if mesh[0][j][0] < mesh[k][j][0]:
                    np.insert(mesh[k][j], 0, mesh[0][j][0])
                if mesh[k][j][-1] < mesh[0][j][-1]:
                    np.append(mesh[k][j], mesh[0][j][-1])

                # Extend P_ssa so that P and P_ssa are evaluated on the same mesh
            P_new = [np.interp(mesh[0][j], mesh[k][j], P[k][j]) for k in range(len(mesh))]
            plotP1D(ax, P_new, mesh[0][j], label, xlabel=xlabel, ylabel=ylabel)

            # Calculate maximal difference and print it in title
            if stats == True:
                title = ""
                for k in range(len(mesh) - 1):
                    max_error = np.max(np.abs(P[k][j] - P[-1][j]))
                    if np.isnan(max_error):
                        title += "$\mathrm{{max}}(|\\textrm{{{}}}-\\textrm{{{}}}|)$ = NaN\n".format(label[k], label[-1])
                    else:
                        error_str_split = "{:.2e}".format(max_error).split('e')
                        error_str = "{} \\cdot 10^{{{}}}".format(error_str_split[0], int(error_str_split[1]))
                        title += "$\mathrm{{max}}(|\\textrm{{{}}}-\\textrm{{{}}}|)$ = ${}$\n".format(label[k], label[-1], error_str)
                ax.set_title(title)
