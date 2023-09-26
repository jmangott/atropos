import numpy as np
from typing import Union

from grid_class import GridParms
from id_class import Id

class Node:
    left = None
    right = None
    def __init__(self, _id: Id, _grid: GridParms):
        self.id = _id
        self.grid = _grid

class Root(Node):
    def __init__(self, _grid: GridParms, _r: int):
        super().__init__(Id(""), _grid)
        self.r = _r
        self.S = np.zeros((_r, _r))

class InternalNode(Node):
    def __init__(self, _parent: Union[Root, 'InternalNode'], _id: Id, _grid: GridParms, _r: int):
        super().__init__(_id, _grid)
        self.parent = _parent
        self.r = _r
        self.Q = np.zeros((self.parent.r, self.r, self.r))
        self.S = np.zeros((self.r, self.r))

class ExternalNode(Node):
    def __init__(self, _parent: Union[Root, 'InternalNode'], _id: Id, _grid: GridParms):
        super().__init__(_id, _grid)
        self.parent = _parent
        self.X = np.zeros((self.parent.r, self.grid.dx))

class Tree:
    def __init__(self, _root: Root):
        self.root = _root

    def __printTree(self, node: Node):
        if isinstance(node, ExternalNode):
            print(type(node), "id:", node.id, "n:", node.grid)
        else:
            print(type(node), "id:", node.id, "n:", node.grid, "r:", node.r)
        if node.left:
            self.__printTree(node.left)
        if node.right:
            self.__printTree(node.right)
    
    def printTree(self):
        self.__printTree(self.root)