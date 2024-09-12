"""Helper module for postprocessing the DLR results."""

import networkx as nx
import numpy as np
import glob
import os
import xarray as xr

if not os.path.exists("plots"):
    os.makedirs("plots")

from src.index_functions import vecIndexToState
from src.grid import GridParms
from src.tree import Tree, Node, ExternalNode, InternalNode


def groupPath(id) -> str:
    path = ""
    for i, id_element in enumerate(id):
        path = path + "/" + id[:i] + id_element
    return path[1:]


def __readTree(node: Node, filename: str):
    gp = groupPath(node.id)
    if isinstance(node, ExternalNode):
        with xr.open_dataset(filename, group=gp) as ds:
            node.X = ds["X"].values.T

    elif isinstance(node, InternalNode):
        with xr.open_dataset(filename, group=gp) as ds:
            node.Q = ds["Q"].values.T
            __readTree(node.child[0], filename)
            __readTree(node.child[1], filename)
    return


def readTree(filename: str) -> Tree:
    with xr.open_dataset(filename) as ds:
        partition_str = ds.attrs["partition_str"]
        species = ds["species"].values
        sorted = np.argsort(species)
        n = ds["n"].values[sorted]
        binsize = ds["binsize"].values[sorted]
        liml = ds["liml"].values[sorted]
        grid = GridParms(n, binsize, liml)
        tree = Tree(partition_str, grid)

        tree.species_names = ds["species_names"].values
        tree.root.Q = ds["Q"].values.T
        __readTree(tree.root.child[0], filename)
        __readTree(tree.root.child[1], filename)

        return tree


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


def convertToSeconds(time_string):
    factor = [3600.0, 60.0, 1.0, 0.001]
    unit = ["h", "mins", "s", "ms"]
    seconds = 0.0
    for ts, f, u in zip(time_string, factor, unit):
        seconds += float(ts[: -len(u)]) * f
    return seconds


class TimeSeries:
    def __init__(self, _foldername):
        if os.path.exists(_foldername):
            self.foldername = _foldername
            self.time = []
            self.__list_of_files = sorted(
                glob.glob(_foldername + "/*.nc"), key=self.__getT
            )
            self.time.sort()
            self.__number_of_files = len(self.__list_of_files)
        else:
            raise Exception("`_foldername` does not exist")

    def __getT(self, filename):
        with xr.open_dataset(filename) as ds:
            t = float(ds["t"].values)
            self.time.append(t)
            return t

    def getMaxMassErr(self):
        with open(self.foldername + "/diagnostics.txt") as file:
            for line in file:
                if line.startswith("max(norm - 1.0):"):
                    return float(line.split()[-1])

    def getWallTime(self):
        with open(self.foldername + "/diagnostics.txt") as file:
            for line in file:
                if line.startswith("Time elapsed:"):
                    return convertToSeconds(line.split()[-4:])

    def getTau(self):
        with xr.open_dataset(self.__list_of_files[0]) as ds:
            return float(ds["tau"].values)

    def getD(self):
        with xr.open_dataset(self.__list_of_files[0]) as ds:
            return ds["n"].values.size

    def getSpeciesNames(self):
        with xr.open_dataset(self.__list_of_files[0]) as ds:
            return ds["species_names"].values

    def getDx(self):
        with xr.open_dataset(self.__list_of_files[0]) as ds:
            return np.prod(ds["n"].values.astype(int))

    def getMassErr(self):
        mass_error = np.zeros(self.__number_of_files)
        for i, filename in enumerate(self.__list_of_files):
            with xr.open_dataset(filename) as ds:
                mass_error[i] = float(ds["dm"].values)
        return mass_error

    def calculateMoments(self):
        n_moments = 2
        moments = [
            {
                name: np.zeros((self.__number_of_files))
                for name in self.getSpeciesNames()
            }
            for _ in range(n_moments)
        ]

        for i, filename in enumerate(self.__list_of_files):
            tree = readTree(filename)
            slice_vec = np.zeros(tree.grid.d(), dtype="int")
            _, marginal_distribution = tree.calculateObservables(slice_vec)
            for j, n_j in enumerate(tree.grid.n):
                name = tree.species_names[j]
                state = vecIndexToState(
                    np.arange(n_j), tree.grid.liml[j], tree.grid.binsize[j]
                )
                for m in range(n_moments):
                    moments[m][name][i] = np.dot(
                        marginal_distribution[name], state ** (m + 1)
                    )
        return moments

    def calculateFullDistribution(self):
        P = np.zeros((self.__number_of_files, self.getDx()))

        for i, filename in enumerate(self.__list_of_files):
            tree = readTree(filename)
            P[i, :] = tree.calculateFullDistribution()

        return P


def calculateDistributionError(
    filename, ref_sliced_distribution, ref_marginal_distribution, slice_vec, ssa_sol
):
    tree = readTree(filename)
    sliced_err = np.zeros(tree.grid.d())
    marginal_err = np.zeros(tree.grid.d())
    sliced, marginal = tree.calculateObservables(slice_vec)
    for i in range(tree.grid.d()):
        sliced_err[i] = np.linalg.norm(
            sliced[i][ssa_sol.n_min[i] : ssa_sol.n_min[i] + ssa_sol.n[i]]
            - ref_sliced_distribution[i][: tree.grid.n[i]]
        )  # Frobenius norm
        marginal_err[i] = np.linalg.norm(
            marginal[i][ssa_sol.n_min[i] : ssa_sol.n_min[i] + ssa_sol.n[i]]
            - ref_marginal_distribution[i][: tree.grid.n[i]]
        )  # Frobenius norm
    return sliced_err, marginal_err


def printEntropyCuts(tree: Tree):
    entropy = tree.calculateEntropy(tree.root)
    entropy0 = tree.calculateEntropy(tree.root.child[0])
    entropy1 = tree.calculateEntropy(tree.root.child[1])
    total_entropy = entropy + entropy0 + entropy1

    cuts = nx.cut_size(
        tree.G,
        [tree.species_names[s] for s in tree.root.child[0].grid.species],
        [tree.species_names[s] for s in tree.root.child[1].grid.species],
    )
    G0 = nx.subgraph(
        tree.G, [tree.species_names[s] for s in tree.root.child[0].grid.species]
    )
    cuts0 = nx.cut_size(
        G0,
        [tree.species_names[s] for s in tree.root.child[0].child[0].grid.species],
        [tree.species_names[s] for s in tree.root.child[0].child[1].grid.species],
    )
    G1 = nx.subgraph(
        tree.G, [tree.species_names[s] for s in tree.root.child[1].grid.species]
    )
    cuts1 = nx.cut_size(
        G1,
        [tree.species_names[s] for s in tree.root.child[1].child[0].grid.species],
        [tree.species_names[s] for s in tree.root.child[1].child[1].grid.species],
    )
    total_cuts = cuts + cuts0 + cuts1

    print(
        "Total entropy: {}, total number of cuts: {}".format(total_entropy, total_cuts)
    )
