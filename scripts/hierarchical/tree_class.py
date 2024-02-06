"""
Contains the `Tree` class for storing the low-rank approximation od the initial probability distribution as a binary tree according to a prescribed partition and writing it to a netCDF file.
"""
import copy
from datatree import DataTree
import numpy as np
import numpy.typing as npt
import os
import regex
from typing import Union
import xarray as xr

from scripts.hierarchical.grid_class import GridParms
from scripts.hierarchical.id_class import Id
from scripts.hierarchical.reaction_class import ReactionSystem
from scripts.index_functions import incrVecIndex, vecIndexToCombIndex

class Node:
    def __init__(self, _id: Id, _grid: GridParms):
        self.child = [None] * 2
        self.id = _id
        self.grid = _grid


class InternalNode(Node):
    def __init__(self, _parent: 'InternalNode', _id: Id, _grid: GridParms):
        super().__init__(_id, _grid)
        self.parent = _parent
        self.Q = np.zeros((0, 0, 0))

    def rankIn(self):
        return self.Q.shape[-1]

    def rankOut(self):
        return self.Q.shape[0]

class ExternalNode(Node):
    def __init__(self, _parent: InternalNode, _id: Id, _grid: GridParms):
        super().__init__(_id, _grid)
        self.parent = _parent
        self.X = np.zeros((0, 0))
        self.propensity = []

    def rankIn(self):
        return self.X.shape[-1]


"""
The `tree` class stores the initial condition for the hierarchical DLR approximation of the chemical master equation. Note that the input for the grid parameters (`_grid`) must follow the same ordering convention for the species as the reaction system (`_reaction_system`), when initializing the tree via the `initialize` method, but the partition string for dividing the reaction networks allows also permutation.
"""
class Tree:
    # TODO: Move the `_grid` dependency to `initialize` method
    def __init__(self,
                 _partition_str: str,
                 _grid: GridParms
                ):

        # Test whether `_partition_str` is a valid input string
        p = self.__removeBrackets(_partition_str)
        
        self.species = np.copy(p)  # Create a deep copy of `p` before sorting it
        p.sort()

        p_diff = np.diff(p)
        if not regex.fullmatch(r"\((?:\d+(?:\s\d)*|(?R))+\)\((?:\d+(?:\s\d)*|(?R))+\)", _partition_str):
            raise Exception("Not a valid `_partition_str`")
        if p[0] != 0:
            raise Exception("The smallest value of `_partition_str` has to be 0")
        elif np.any(p_diff != 1):
            raise Exception("Not all species covered by `_partition_str`")

        # Test whether parentheses are closed
        sigma = 0
        for ele in _partition_str:
            sigma += self.__parsingHelper(ele)
            if sigma < 0:
                raise Exception("Left bracket is missing in `_partition_str`")
        if sigma != 0:
            raise Exception("Right bracket is missing is `_partition_str`")

        # Test whether the dimension of `_partition_str` and `_grid.d` match
        if p.size != _grid.d():
            raise Exception(
                "Dimensions of `_partition_str` and `_grid.d` do not match")

        self.n_internal_nodes = 1  # 1 is for the root node
        self.n_external_nodes = 0

        self.reaction_system = None
        self.partition_str = _partition_str
        self.grid = copy.deepcopy(_grid)
        self.root = InternalNode(None, Id(""), self.grid.permute(self.species))
        self.__build(self.root, self.partition_str)

    @staticmethod
    def __parsingHelper(input_str):
        if (input_str == "("):
            sigma = 1
        elif (input_str == ")"):
            sigma = -1
        else:
            sigma = 0
        return sigma

    @staticmethod
    def __removeBrackets(input_str):
        n = input_str.replace("(", " ").replace(")", " ").split()
        n = np.array([int(ele) for ele in n])
        return n

    def __build(self, node: InternalNode, partition_str):
        sigma = 0
        i = 0
        for i, ele in enumerate(partition_str):
            sigma += self.__parsingHelper(ele)
            if sigma == 0:
                break

        partition_str0 = partition_str[1:i]
        partition_str1 = partition_str[i+2:-1]
        p0_size = self.__removeBrackets(partition_str0).size

        grid0 = GridParms(
            node.grid.n[:p0_size], node.grid.binsize[:p0_size], node.grid.liml[:p0_size], node.grid.species[:p0_size])
        grid1 = GridParms(
            node.grid.n[p0_size:], node.grid.binsize[p0_size:], node.grid.liml[p0_size:], node.grid.species[p0_size:])
                
        if (partition_str0[0] == "("):
            node.child[0] = InternalNode(node, node.id + 0, grid0)
            self.__build(node.child[0], partition_str0)
            self.n_internal_nodes += 1
        else:
            node.child[0] = ExternalNode(node, node.id + 0, grid0)
            self.n_external_nodes += 1

        if (partition_str1[0] == "("):
            node.child[1] = InternalNode(node, node.id + 1, grid1)
            self.__build(node.child[1], partition_str1)
            self.n_internal_nodes += 1
        else:
            node.child[1] = ExternalNode(node, node.id + 1, grid1)
            self.n_external_nodes += 1

        return

    @staticmethod
    def __calculatePropensity(node: ExternalNode, reaction_system: ReactionSystem):
        node.propensity = [None] * reaction_system.size()
        for mu, reaction in enumerate(reaction_system.reactions):
            n_dep = node.grid.n[node.grid.dep[:, mu]]
            dx_dep = np.prod(n_dep)
            node.propensity[mu] = np.ones(dx_dep)
            vec_index = np.zeros(n_dep.size)
            reactants = [reactant for reactant in node.grid.species if reactant in reaction.propensity.keys()]
            for i in range(dx_dep):
                for j, reactant in enumerate(reactants):
                    node.propensity[mu][i] *= reaction.propensity[reactant](vec_index[j])
                incrVecIndex(vec_index, n_dep, n_dep.size)
    
    def __initialize(self, node: Node, reaction_system: ReactionSystem):
        p0_size = node.parent.child[0].grid.d()

        idx = node.parent.child.index(node)
        if idx == 0:
            sl = slice(0,p0_size)
        elif idx == 1:
            sl = slice(p0_size, node.parent.grid.d())

        node.grid = GridParms(
            node.parent.grid.n[sl], node.parent.grid.binsize[sl], node.parent.grid.liml[sl], node.parent.grid.species[sl], node.parent.grid.dep[sl, :], node.parent.grid.nu[sl, :])

        if isinstance(node, InternalNode):
            next_r_out = next(self.__r_out_iter)
            node.Q.resize((next_r_out, next_r_out, node.parent.rankOut()))
            self.__initialize(node.child[0], reaction_system)
            self.__initialize(node.child[1], reaction_system)

        if isinstance(node, ExternalNode):
            node.X.resize((node.grid.dx(), node.parent.rankOut()))
            self.__calculatePropensity(node, reaction_system)

        return

    def initialize(self, reaction_system: ReactionSystem, r_out: npt.NDArray[np.int_]):
        # Test whether the dimension of `_r` is equal to n_internal_nodes
        if r_out.size != self.n_internal_nodes:
            raise Exception(
                "`r_out.size` must be equal to the number of internal nodes")

        self.grid.initialize(reaction_system)
        self.root.grid = self.grid.permute(self.species)
        self.__r_out_iter = iter(r_out)
        next_r_out = next(self.__r_out_iter)
        self.root.Q.resize((next_r_out, next_r_out, 1))
        self.__initialize(self.root.child[0], reaction_system)
        self.__initialize(self.root.child[1], reaction_system)

    def __print(self, node: Node, os: str) -> str:
        os = " ".join([os, str(type(node)), "id:", str(node.id), "n:", str(node.grid), "species:", str(node.grid.species)])
        if isinstance(node, ExternalNode):
            os = " ".join([os, "X.shape:", str(node.X.shape), "\n"])
        elif isinstance(node, InternalNode):
            os = " ".join([os, "Q.shape:", str(node.Q.shape), "\n"])
        if node.child[0]:
            os = self.__print(node.child[0], os)
        if node.child[1]:
            os = self.__print(node.child[1], os)
        return os

    def __str__(self) -> str:
        return self.__print(self.root, "")

    @staticmethod
    def __createDataset(node: Node):
        ds = xr.Dataset(
            {
                "n": (["d"], node.grid.n),
                "binsize": (["d"], node.grid.binsize),
                "liml": (["d"], node.grid.liml),
                "species": (["d"], node.grid.species),
                "dep": (["d", "n_reactions"], node.grid.dep),
                "nu": (["d", "n_reactions"], node.grid.nu)
            }
        )
        return ds

    def __write(self, node: Node, parent_dt: DataTree):
        if isinstance(node, ExternalNode):
            ds = self.__createDataset(node)
            ds["X"] = (["n_basisfunctions", "dx"], node.X.T)
            for mu, propensity in enumerate(node.propensity):
                ds["propensity_{}".format(mu)] = (["dx_{}".format(mu)], propensity)
            DataTree(name=str(node.id),
                     parent=parent_dt, data=ds)

        elif isinstance(node, InternalNode):
            ds = self.__createDataset(node)
            ds["Q"] = (["n_basisfunctions", "r_out", "r_out"], node.Q.T)
            dt = DataTree(name=str(node.id),
                          parent=parent_dt, data=ds)
            self.__write(node.child[0], dt)
            self.__write(node.child[1], dt)

        return

    def write(self, fname: str = "input/input.nc"):
        if not os.path.exists("input"):
            os.makedirs("input")

        # Undo permutation of grid for root only
        self.grid.permute(self.species)

        ds = self.__createDataset(self.root)
        ds["Q"] = (["n_basisfunctions", "r_out", "r_out"], self.root.Q.T)
        dt = DataTree(name=str(self.root.id), data=ds)
        dt.attrs["partition_str"] = self.partition_str
        self.__write(self.root.child[0], dt)
        self.__write(self.root.child[1], dt)
        dt.to_netcdf(fname)

    def __calculateObservables(self, node: Node, slice: npt.NDArray[np.int_]):
        if isinstance(node, ExternalNode):
            node.X_sum = np.sum(node.X, axis=0)
            node.X_slice = node.X[vecIndexToCombIndex(slice[node.grid.species], node.grid.n), :]
        elif isinstance(node, InternalNode):
            self.__calculateObservables(node.child[0], slice)
            self.__calculateObservables(node.child[1], slice)
            node.X_sum = np.einsum("i,ijk,j", node.child[0].X_sum, node.Q, node.child[1].X_sum)
            node.X_slice = np.einsum("i,ijk,j", node.child[0].X_slice, node.Q, node.child[1].X_slice)

    def calculateObservables(self, slice: npt.NDArray[np.int_]):
        if slice.size != self.root.grid.d():
            raise Exception(
                "`slice.size` must be equal to the number of dimensions of the root node")
        self.__calculateObservables(self.root.child[0], slice)
        self.__calculateObservables(self.root.child[1], slice)