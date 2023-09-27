from datatree import DataTree
import numpy.typing as npt
import os
import regex
import xarray as xr

from tree_class import *

class Partition:
    def __init__(self, 
                 _partition_str: str, 
                 _grid: GridParms, 
                 _r: npt.NDArray[np.int_]):

        # Test whether `_partition_str` is a valid input string
        p = self.__removeBrackets(_partition_str)
        p_diff = np.diff(p)
        if not regex.fullmatch(r"\((?:\d+(?:\s\d)*|(?R))+\)\((?:\d+(?:\s\d)*|(?R))+\)", _partition_str):
            raise ValueError("Not a valid input string")
        if p[0] != 0:
            raise ValueError("The first numerical value has to be 0")
        elif np.any(p_diff != 1):
            raise ValueError("Numerical values have to be in ascending order")

        # Test whether parentheses are closed
        sigma = 0
        for ele in _partition_str:
            sigma += self.__parsingHelper(ele)
            if sigma < 0:
                raise ValueError("Left bracket is missing in `_partition_str`")
        if sigma != 0:
            raise ValueError("Right bracket is missing is `_partition_str`")
            
        # Test whether the dimension of `_partition_str` and `_grid.d` match
        if p.size != _grid.d:
            raise ValueError("Dimension of `_partition_str` and `_grid` do not match")

        # Calculate the number of internal nodes
        n_internal_nodes = self.__countInternalNodes(_partition_str, 0)

        # Test whether the dimension of `_r` is equal to n_internal_nodes + 1
        if _r.size != n_internal_nodes + 1:
            raise ValueError("`_r.size` must be equal to the number of internal nodes + 1")

        self.partition_str = _partition_str
        self.grid = _grid
        self.r = _r
        self.tree = Tree(Root(self.grid, self.r[0]))

    @staticmethod
    def __parsingHelper(input_str):
        if (input_str == "("):
            sigma = 1
        elif (input_str == ")"):
            sigma =-1
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

    def __buildTree(self, node, partition_str, r_iter):
        sigma = 0
        for i, ele in enumerate(partition_str):
            sigma += self.__parsingHelper(ele)
            if sigma == 0:
                break

        partition_str0 = partition_str[1:i]
        partition_str1 = partition_str[i+2:-1]

        p0 = self.__removeBrackets(partition_str0)

        grid0 = GridParms(
            node.grid.n[:p0.size], node.grid.binsize[:p0.size], node.grid.liml[:p0.size], node.grid.dep[:, :p0.size])
        grid1 = GridParms(
            node.grid.n[p0.size:], node.grid.binsize[p0.size:], node.grid.liml[p0.size:], node.grid.dep[:, p0.size:])

        if (partition_str0[0] == "("):
            node.left = InternalNode(node, node.id + 0, grid0, next(r_iter))
            self.__buildTree(node.left, partition_str0, r_iter)
        else:
            node.left = ExternalNode(node, node.id + 0, grid0)

        if (partition_str1[0] == "("):
            node.right = InternalNode(node, node.id + 1, grid1, next(r_iter))
            self.__buildTree(node.right, partition_str1, r_iter)
        else:
            node.right = ExternalNode(node, node.id + 1, grid1)
        return

    def buildTree(self):
        r_iter = iter(self.r[1:])
        self.__buildTree(self.tree.root, self.partition_str, r_iter)

    def initializeTree(self):
        ...

    def printTree(self):
        self.tree.printTree()

    @staticmethod
    def __createInternalDataset(node: InternalNode):
        ds = xr.Dataset(
            {
                "Q": (["parent_r", "r", "r"], node.Q),
                "S": (["r", "r"], node.S),
                "n": (["d"], node.grid.n),
                "binsize": (["d"],node.grid.binsize),
                "liml": (["d"], node.grid.liml),
                "dep": (["mu", "d"], node.grid.dep)
            }
        )
        return ds
    
    @staticmethod
    def __createExternalDataset(node: ExternalNode):
        ds = xr.Dataset(
            {
                "X": (["n_basisfunctions", "dx"], node.X),
                "n": (["d"], node.grid.n),
                "binsize": (["d"],node.grid.binsize),
                "liml": (["d"], node.grid.liml),
                "dep": (["mu", "d"], node.grid.dep)
            }
        )
        return ds

    def __writeTree(self, node, parent_dt):
        if isinstance(node.left, ExternalNode):
            ds = self.__createExternalDataset(node.left)
            DataTree(name=str(node.left.id),
                     parent=parent_dt, data=ds)
        else:
            ds = self.__createInternalDataset(node.left)
            dt = DataTree(name=str(node.left.id),
                          parent=parent_dt, data=ds)
            self.__writeTree(node.left, dt)
        
        if isinstance(node.right, ExternalNode):
            ds = self.__createExternalDataset(node.right)
            DataTree(name=str(node.right.id),
                     parent=parent_dt, data=ds)
        else:
            ds = self.__createInternalDataset(node.right)
            dt = DataTree(name=str(node.right.id),
                          parent=parent_dt, data=ds)
            self.__writeTree(node.right, dt)
        return
    
    def writeTree(self, fname: str = "input/input.nc"):
        if not os.path.exists("input"):
            os.makedirs("input")

        ds = xr.Dataset(
            {
                "S": (["r", "r"], self.tree.root.S),
                "n": (["d"], self.tree.root.grid.n),
                "binsize": (["d"], self.tree.root.grid.binsize),
                "liml": (["d"], self.tree.root.grid.liml),
                "dep": (["mu", "d"], self.tree.root.grid.dep)
            }
        )
        dt = DataTree(name=str(self.tree.root.id), data=ds)
        self.__writeTree(self.tree.root, dt)
        dt.to_netcdf(fname)

if __name__ == "__main__":
    partition_str = "((0 1)(2 3))((4 5)((6)(7 8)))"
    r = np.array([5, 4, 3, 2])
    mu = 5
    n = np.array([4, 6, 7, 8, 3, 5, 2, 9, 1])
    binsize = np.ones(n.size, dtype=int)
    liml = np.zeros(n.size)
    dep = np.ones((mu, n.size), dtype=bool)
    grid = GridParms(n, binsize, liml, dep)
    j = Partition(partition_str, grid, r)
    j.buildTree()
    j.printTree()
    j.writeTree()