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
from scripts.index_functions import incrVecIndex

class Node:
    def __init__(self, _id: Id, _grid: GridParms, _species: list):
        self.child = [None] * 2
        self.id = _id
        self.grid = _grid
        self.species = _species
        self.propensity = []

class InternalNode(Node):
    def __init__(self, _parent: 'InternalNode', _id: Id, _grid: GridParms, _species_limits: list, _r_in: int, _r_out: int):
        super().__init__(_id, _grid, _species_limits)
        self.parent = _parent
        self.r_in = _r_in
        self.r_out = _r_out
        self.Q = np.zeros((self.r_out, self.r_out, self.r_in))

class ExternalNode(Node):
    def __init__(self, _parent: InternalNode, _id: Id, _grid: GridParms, _species_limits: list):
        super().__init__(_id, _grid, _species_limits)
        self.parent = _parent
        self.r_in = self.parent.r_out
        self.X = np.zeros((self.grid.dx(), self.r_in))


"""
The `tree` class stores the initial condition for the hierarchical DLR approximation of the chemical master equation. Note that the input for the grid parameters (`_grid`) must follow the same ordering convention for the species as the reaction system (`_reaction_system`), but the partition string for dividing the reaction networks allows also permutation.
"""
class Tree:
    def __init__(self,
                 _reaction_system: ReactionSystem,
                 _partition_str: str,
                 _grid: GridParms,
                 _r_out: npt.NDArray[np.int_]):

        # Test whether `_partition_str` is a valid input string
        p = self.__removeBrackets(_partition_str)
        
        species = np.copy(p)  # Create a deep copy of `p` before sorting it
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

        # Calculate the number of internal nodes
        self.n_internal_nodes = self.__countInternalNodes(_partition_str, 1)

        # Test whether the dimension of `_r` is equal to n_internal_nodes
        if _r_out.size != self.n_internal_nodes:
            raise Exception(
                "`_r_out.size` must be equal to the number of internal nodes")
        
        # Calculate the number of external nodes
        self.n_external_nodes = len(regex.findall(r"\(\d+(?:\s\d)*\)", _partition_str)) + 1

        self.reaction_system = _reaction_system
        self.partition_str = _partition_str
        self.grid = copy.deepcopy(_grid)
        self.r_out = _r_out
        self.grid.initialize(self.reaction_system)
        self.grid.permute(species)
        self.root = InternalNode(None, Id(""), self.grid, species, 1, self.r_out[0])

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

    def __countInternalNodes(self, partition_str, n):
        sigma = 0
        for i, ele in enumerate(partition_str):
            sigma += self.__parsingHelper(ele)
            if sigma == 0:
                break

        partition_str0 = partition_str[1:i]
        partition_str1 = partition_str[i+2:-1]

        if (partition_str0[0] == "("):
            n = self.__countInternalNodes(partition_str0, n + 1)
        if (partition_str1[0] == "("):
            n = self.__countInternalNodes(partition_str1, n + 1)
        return n

    def __buildTree(self, node, partition_str, r_out_iter):
        sigma = 0
        for i, ele in enumerate(partition_str):
            sigma += self.__parsingHelper(ele)
            if sigma == 0:
                break

        partition_str0 = partition_str[1:i]
        partition_str1 = partition_str[i+2:-1]
        species0 = self.__removeBrackets(partition_str0)
        species1 = self.__removeBrackets(partition_str1)
        p0_size = species0.size

        grid0 = GridParms(
            node.grid.n[:p0_size], node.grid.binsize[:p0_size], node.grid.liml[:p0_size], node.grid.dep[:p0_size, :], node.grid.nu[:p0_size, :])
        grid1 = GridParms(
            node.grid.n[p0_size:], node.grid.binsize[p0_size:], node.grid.liml[p0_size:], node.grid.dep[p0_size:, :], node.grid.nu[p0_size:, :])
        
        if (partition_str0[0] == "("):
            node.child[0] = InternalNode(node, node.id + 0, grid0, species0, node.r_out, next(r_out_iter))
            self.__buildTree(node.child[0], partition_str0, r_out_iter)
        else:
            node.child[0] = ExternalNode(node, node.id + 0, grid0, species0)

        if (partition_str1[0] == "("):
            node.child[1] = InternalNode(node, node.id + 1, grid1, species1, node.r_out, next(r_out_iter))
            self.__buildTree(node.child[1], partition_str1, r_out_iter)
        else:
            node.child[1] = ExternalNode(node, node.id + 1, grid1, species1)

        self.__calculatePropensity(node.child[0], self.reaction_system)
        self.__calculatePropensity(node.child[1], self.reaction_system)
        return

    @staticmethod
    def __calculatePropensity(node: Node, reaction_system: ReactionSystem):
        node.propensity = [None] * reaction_system.size()
        for mu, reaction in enumerate(reaction_system.reactions):
            n_dep = np.array([n for i, n in enumerate(node.grid.n) if node.grid.dep[i, mu]], dtype="int")
            dx_dep = np.prod(n_dep)
            node.propensity[mu] = np.ones(dx_dep)
            vec_index = np.zeros(n_dep.size)
            reactants = [reactant for reactant in node.species if reactant in reaction.propensity.keys()]
            for i in range(dx_dep):
                for j, reactant in enumerate(reactants):
                    node.propensity[mu][i] *= reaction.propensity[reactant](vec_index[j])
                incrVecIndex(vec_index, n_dep, n_dep.size)

    def buildTree(self):
        r_out_iter = iter(self.r_out[1:])
        self.__calculatePropensity(self.root, self.reaction_system)
        self.__buildTree(self.root, self.partition_str, r_out_iter)

    def __printTree(self, node: Node, os: str) -> str:
        os = " ".join([os, str(type(node)), "id:", str(node.id), "n:", str(node.grid), "species:", str(node.species)])
        if isinstance(node, ExternalNode):
            os = " ".join([os, "X.shape:", str(node.X.shape), "\n"])
        else:
            os = " ".join([os, "Q.shape:", str(node.Q.shape), "\n"])
        if node.child[0]:
            os = self.__printTree(node.child[0], os)
        if node.child[1]:
            os = self.__printTree(node.child[1], os)
        return os

    def __str__(self) -> str:
        return self.__printTree(self.root, "")

    @staticmethod
    def __createInternalDataset(node: InternalNode):
        ds = xr.Dataset(
            {                
                "Q": (["n_basisfunctions", "r_out", "r_out"], node.Q.T),
                "n": (["d"], node.grid.n),
                "binsize": (["d"], node.grid.binsize),
                "liml": (["d"], node.grid.liml),
                "dep": (["d", "n_reactions"], node.grid.dep),
                "nu": (["d", "n_reactions"], node.grid.nu)
            }
        )
        for mu, propensity in enumerate(node.propensity):
            ds["propensity_{}".format(mu)] = (["dx_{}".format(mu)], propensity)
        return ds

    @staticmethod
    def __createExternalDataset(node: ExternalNode):
        ds = xr.Dataset(
            {
                "X": (["n_basisfunctions", "dx"], node.X.T),
                "n": (["d"], node.grid.n),
                "binsize": (["d"], node.grid.binsize),
                "liml": (["d"], node.grid.liml),
                "dep": (["d", "n_reactions"], node.grid.dep),
                "nu": (["d", "n_reactions"], node.grid.nu)
            }
        )
        for mu, propensity in enumerate(node.propensity):
            ds["propensity_{}".format(mu)] = (["dx_{}".format(mu)], propensity)
        return ds

    def __writeTree(self, node: Node, parent_dt: DataTree):
        if isinstance(node.child[0], ExternalNode):
            ds = self.__createExternalDataset(node.child[0])
            DataTree(name=str(node.child[0].id),
                     parent=parent_dt, data=ds)
        else:
            ds = self.__createInternalDataset(node.child[0])
            dt = DataTree(name=str(node.child[0].id),
                          parent=parent_dt, data=ds)
            self.__writeTree(node.child[0], dt)

        if isinstance(node.child[1], ExternalNode):
            ds = self.__createExternalDataset(node.child[1])
            DataTree(name=str(node.child[1].id),
                     parent=parent_dt, data=ds)
        else:
            ds = self.__createInternalDataset(node.child[1])
            dt = DataTree(name=str(node.child[1].id),
                          parent=parent_dt, data=ds)
            self.__writeTree(node.child[1], dt)
        return

    def writeTree(self, fname: str = "input/input.nc"):
        if not os.path.exists("input"):
            os.makedirs("input")

        ds = self.__createInternalDataset(self.root)

        dt = DataTree(name=str(self.root.id), data=ds)
        self.__writeTree(self.root, dt)
        dt.to_netcdf(fname)