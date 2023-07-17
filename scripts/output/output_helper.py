"""Helper module for postprocessing the DLR results."""
# TODO: for large netCDF files, this is VERY slow, 
# maybe precompute the 1D marginal and sliced distributions?
import matplotlib as mpl
mpl.use('pgf')
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import os

if not os.path.exists("plots"):
    os.makedirs("plots")

from scripts.index_functions import IncrVecIndex, VecIndexToCombIndex

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
    def __init__(self, ds: nc._netCDF4.Dataset):
        self.n1 = ds["n1"][:]
        self.n2 = ds["n2"][:]
        self.dx1 = np.prod(self.n1)
        self.dx2 = np.prod(self.n2)
        self.n = np.concatenate((self.n1, self.n2))
        self.m1 = ds.dimensions["m1"].size
        self.m2 = ds.dimensions["m2"].size
        self.r = ds.dimensions["r"].size
        self.d = self.m1 + self.m2
        self.bin = ds["binsize"][:]
        self.liml = np.asarray(ds["liml"][:], dtype="int64")
        self.t = ds["t"][0]
        self.dt = ds["dt"][0]


class Observables:
    def __init__(self, _ds: nc._netCDF4.Dataset, _grid: GridInfo):
        """
        Helper class for calculating marginal and sliced distributions and the full probability distribution.
        """
        self.X1 = _ds['X']
        self.S = _ds['S']
        self.X2 = _ds['V']
        self.grid = _grid


    def fullDistribution(self) -> np.ndarray:
        """
        Calculates full probability distribution.
        NOTE: This function should be used only for small systems.
        """
        P_full = np.matmul(np.transpose(self.X1), np.matmul(self.S, self.X2))
        
        return P_full.flatten(order='F')


    def marginalDistributions(self) -> list[np.ndarray]:
        """Calculates marginal distributions."""
        X1_sum = [np.zeros((n_el, self.grid.r), dtype="float64") for n_el in self.grid.n1]
        X2_sum = [np.zeros((n_el, self.grid.r), dtype="float64") for n_el in self.grid.n2]
        P_marginal = [np.zeros(n_el) for n_el in self.grid.n]

        vec_index = np.zeros(self.grid.m1, dtype="int64")
        for i in range(self.grid.dx1):
            for k_vec, vec_el in enumerate(vec_index):
                X1_sum[k_vec][vec_el, :] += self.X1[:, i]
            IncrVecIndex(vec_index, self.grid.n1, self.grid.m1)

        vec_index = np.zeros(self.grid.m2, dtype="int64")
        for i in range(self.grid.dx2):
            for k_vec, vec_el in enumerate(vec_index):
                X2_sum[k_vec][vec_el, :] += self.X2[:, i]
            IncrVecIndex(vec_index, self.grid.n2, self.grid.m2)

        X1_sum_tot = np.sum(X1_sum[0], axis=0)
        X2_sum_tot = np.sum(X2_sum[0], axis=0)

        for i in range(self.grid.m1):
            P_marginal[i] = np.matmul(X1_sum[i], np.matmul(self.S, np.transpose(X2_sum_tot)))

        for i in range(self.grid.m2):
            P_marginal[self.grid.m1 + i] = np.matmul(X1_sum_tot, np.matmul(self.S, np.transpose(X2_sum[i])))
        
        return P_marginal


    def marginalDistribution2D(self) -> np.ndarray:
        """
        Calculates a 2D marginal distribution. 
        NOTE: This is 'hard-coded' for the first partition and assumes that only two species lie in that partition.
        """
        X2_sum_tot = np.sum(self.X2, axis=1)
        P_marginal2D_vec = np.matmul(np.transpose(self.X1), np.matmul(self.S, X2_sum_tot))
        P_marginal2D = np.reshape(P_marginal2D_vec, self.grid.n1, order="F")

        return P_marginal2D


    def slicedDistributions(self, slice_vec: np.ndarray) -> list[np.ndarray]:
        """Calculates sliced distributions for a given `slice_vec`."""
        X1_slice = [np.zeros((n_el, self.grid.r), dtype="float64") for n_el in self.grid.n1]
        X2_slice = [np.zeros((n_el, self.grid.r), dtype="float64") for n_el in self.grid.n2]
        P_sliced = [np.zeros(n_el) for n_el in self.grid.n]

        slice_vec_index = ((np.asarray(slice_vec, dtype="int64") - self.grid.liml) / self.grid.bin).astype(int)

        for k, n_el in enumerate(self.grid.n1):
            vec_index = slice_vec_index[:self.grid.m1].copy()
            for i in range(n_el):
                vec_index[k] = i
                comb_index = VecIndexToCombIndex(vec_index, self.grid.n1)
                X1_slice[k][i, :] = self.X1[:, comb_index]

        for k, n_el in enumerate(self.grid.n2):
            vec_index = slice_vec_index[self.grid.m1:].copy()
            for i in range(n_el):
                vec_index[k] = i
                comb_index = VecIndexToCombIndex(vec_index, self.grid.n2)
                X2_slice[k][i, :] = self.X2[:, comb_index]

        X1_slice_tot = X1_slice[0][slice_vec_index[0], :]
        X2_slice_tot = X2_slice[0][slice_vec_index[self.grid.m1], :]

        for i in range(self.grid.m1):
            P_sliced[i] = np.matmul(X1_slice[i], np.matmul(self.S, np.transpose(X2_slice_tot)))

        for i in range(self.grid.m2):
            P_sliced[self.grid.m1 + i] = np.matmul(X1_slice_tot, np.matmul(self.S, np.transpose(X2_slice[i])))

        return P_sliced


    def slicedDistribution2D(self, slice_vec: np.ndarray) -> np.ndarray:
        """
        Calculates a 2D sliced distribution for a given `slice_vec`.
        NOTE: This is 'hard-coded' for the first partition and assumes that only two species lie in that partition.
        """
        slice_vec_index = ((np.asarray(slice_vec, dtype="int64") - self.grid.liml) / self.grid.bin).astype(int)
        comb_index = VecIndexToCombIndex(slice_vec_index[self.grid.m1:], self.grid.n2)
        X2_slice_tot = self.X2[:, comb_index]

        P_sliced2D_vec = np.matmul(np.transpose(self.X1), np.matmul(self.S, X2_slice_tot))
        P_sliced2D = np.reshape(P_sliced2D_vec, self.grid.n1, order="F")
        return P_sliced2D


def plotP2D(axs, P: list[np.ndarray], mesh: tuple, title: list[str]) -> tuple[plt.Figure, np.ndarray]:
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


def plotP1Dmult(axs, P: list[np.ndarray], mesh: list, label: list, idx: list, Plabel: str, *, stats: bool = False) -> tuple[plt.Figure, np.ndarray]:
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
                P[k][j] = np.interp(mesh[0][j], mesh[k][j], P[k][j])
            P_new = [P[k][j] for k in range(len(mesh))]
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
