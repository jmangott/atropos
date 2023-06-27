# TODO: add descriptions, generalize `marginalDistribution2D`
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
    "scatter.marker": "2",
    "lines.marker": "2",
    "lines.linestyle": "none",
    "axes.titlesize": "medium" 
})


class grid_info:
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
        self.liml = ds["liml"][:]
        self.t = ds["t"][0]
        self.dt = ds["dt"][0]


def marginalDistributions(ds: nc._netCDF4.Dataset, grid: grid_info) -> list[np.ndarray]:
    """Calculates marginal distributions."""
    X1_sum = [np.zeros((n_el, grid.r), dtype="float64") for n_el in grid.n1]
    X2_sum = [np.zeros((n_el, grid.r), dtype="float64") for n_el in grid.n2]
    P = [np.zeros(n_el) for n_el in grid.n]

    S = ds['S']

    vec_index = np.zeros(grid.m1, dtype="int64")
    for i in range(grid.dx1):
        line = ds['X'][:, i]
        for k_vec, vec_el in enumerate(vec_index):
            X1_sum[k_vec][vec_el, :] += line
        IncrVecIndex(vec_index, grid.n1, grid.m1)

    vec_index = np.zeros(grid.m2, dtype="int64")
    for i in range(grid.dx2):
        line = ds['V'][:, i]
        for k_vec, vec_el in enumerate(vec_index):
            X2_sum[k_vec][vec_el, :] += line
        IncrVecIndex(vec_index, grid.n2, grid.m2)

    X1_sum_tot = np.sum(X1_sum[0], axis=0)
    X2_sum_tot = np.sum(X2_sum[0], axis=0)

    for i in range(grid.m1):
        P[i] = np.matmul(X1_sum[i], np.matmul(S, np.transpose(X2_sum_tot)))

    for i in range(grid.m2):
        P[grid.m1 + i] = np.matmul(X1_sum_tot, np.matmul(S, np.transpose(X2_sum[i])))
    
    return P


def fullDistribution(ds: nc._netCDF4.Dataset, grid: grid_info) -> list[np.ndarray]:
    """
    Calculates full probability distribution.
    NOTE: This function should be used only for small systems.
    """
    S = ds['S']
    X1 = ds['X']
    X2 = ds['V']
    P = np.matmul(np.transpose(X1), np.matmul(S, X2))
    
    return P


def slicedDistributions(ds: nc._netCDF4.Dataset, grid: grid_info, slice_vec: np.ndarray) -> list[np.ndarray]:
    """Calculate sliced distributions for a given `slice_vec`."""
    X1_slice = [np.zeros((n_el, grid.r), dtype="float64") for n_el in grid.n1]
    X2_slice = [np.zeros((n_el, grid.r), dtype="float64") for n_el in grid.n2]
    P = [np.zeros(n_el) for n_el in grid.n]

    S = ds['S']

    slice_vec_index = ((np.asarray(slice_vec, dtype="int64") - grid.liml) / grid.bin).astype(int)

    for k, n_el in enumerate(grid.n1):
        vec_index = slice_vec_index[:grid.m1].copy()
        for i in range(n_el):
            vec_index[k] = i
            comb_index = VecIndexToCombIndex(vec_index, grid.n1)
            line = ds['X'][:, comb_index]
            X1_slice[k][i, :] = line

    for k, n_el in enumerate(grid.n2):
        vec_index = slice_vec_index[grid.m1:].copy()
        for i in range(n_el):
            vec_index[k] = i
            comb_index = VecIndexToCombIndex(vec_index, grid.n2)
            line = ds['V'][:, comb_index]
            X2_slice[k][i, :] = line

    X1_slice_tot = X1_slice[0][slice_vec_index[0], :]
    X2_slice_tot = X2_slice[0][slice_vec_index[grid.m1], :]

    for i in range(grid.m1):
        P[i] = np.matmul(X1_slice[i], np.matmul(S, np.transpose(X2_slice_tot)))

    for i in range(grid.m2):
        P[grid.m1 + i] = np.matmul(X1_slice_tot, np.matmul(S, np.transpose(X2_slice[i])))

    return P


def marginalDistribution2D(ds: nc._netCDF4.Dataset, grid: grid_info) -> np.ndarray:
    """Calculates a 2D marginal distribution."""

    P = np.zeros(grid.dx1)

    X1_sum = np.zeros((grid.dx1, grid.r), dtype="float64")
    X2_sum_tot = np.zeros((grid.r), dtype="float64")

    S = ds['S']

    for i in range(grid.dx1):
        line = ds['X'][:, i]
        X1_sum[i, :] += line

    for i in range(grid.dx2):
        line = ds['V'][:, i]
        X2_sum_tot += line

    for i in range(grid.dx1):
        P[i] = np.matmul(X1_sum[i], np.matmul(S, np.transpose(X2_sum_tot)))

    Pmat = np.reshape(P, grid.n1, order='F')

    return Pmat


def plotP2D(P: np.ndarray, P_ref: np.ndarray, mesh: tuple, title: list[str]) -> tuple[plt.Figure, np.ndarray]:
    fig, ax = plt.subplots(1, 2)
    plt.setp(ax.flat, xlabel='$x_1$', ylabel='$x_2$')
    levels = np.linspace(np.amin([np.amin(P), np.amin(P_ref)]), np.amax(
        [np.amax(P), np.amax(P_ref)]), 9)
    c1 = ax[0].contour(mesh[0], mesh[1], P, levels=levels)
    ax[1].contour(mesh[0], mesh[1], P_ref, levels=levels)
    ax[0].set_title(title[0])
    ax[1].set_title(title[1])
    return fig, ax


def plotP1D(ax, P: np.ndarray, P_ref: np.ndarray, mesh: tuple, label: list, **kwargs) -> tuple[plt.Figure, np.ndarray]:
    plt.setp(ax, **kwargs)
    ax.plot(mesh, P, 'ro', label=label[0])
    ax.plot(mesh, P_ref, 'bx', label=label[1])
    return ax


def plotP1Dmult(axs, P: np.ndarray, P_ref: np.ndarray, grid: any, mesh_ref: any, label: list, idx: list, Plabel: str) -> tuple[plt.Figure, np.ndarray]:
    for i, ax in enumerate(axs.flatten()):
        if i < len(idx):
            j = idx[i]
            xlabel = "$x_{{" + str(j + 1) + "}}$"
            ylabel = "$P_{{\mathrm{{" + Plabel + "}}}}(x_{{" + str(j + 1) + "}})$"
            mesh = grid.bin[j] * range(grid.n[j]) + grid.liml[j]

            # extend mesh_ref so that mesh and mesh_ref cover the same interval
            if mesh[0] < mesh_ref[j][0]:
                np.insert(mesh_ref[j], 0, mesh[0])
            if mesh_ref[j][-1] < mesh[-1]:
                np.append(mesh_ref[j], mesh[-1])

            # extend P_ref so that P and P_ref are evaluated on the same mesh
            P_ref_inter = np.interp(mesh, mesh_ref[j], P_ref[j])

            plotP1D(ax, P[j], P_ref_inter, mesh, label, xlabel=xlabel, ylabel=ylabel)

            # calculate max. difference and print it in title
            max_error = np.max(np.abs(P[j] - P_ref_inter))
            error_str_split = "{:.2e}".format(max_error).split('e')
            error_str = "{} \\times 10^{{{}}}".format(error_str_split[0], int(error_str_split[1]))
            ax.set_title("max. difference = ${}$".format(error_str))
    return axs
