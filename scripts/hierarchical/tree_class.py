import numpy as np
from typing import Union

from grid_class import GridParms

class Node:
    left = None
    right = None
    def __init__(self, _data):
        self.data = _data

    def printTree(self):
        print(self.data)
        if self.left:
            self.left.printTree()
        if self.right:
            self.right.printTree()

class Root(Node):
    def __init__(self, _grid: GridParms, _r: int, _data):
        super().__init__(_data)
        self.grid = _grid
        self.r = _r
        self.S = np.zeros((_r, _r))

class InternalNode(Node):
    def __init__(self, _parent: Union[Root, 'InternalNode'], _grid: GridParms, _r: int, _data):
        super().__init__(_data)
        self.parent = _parent
        self.grid = _grid
        self.r = _r
        self.Q = np.zeros((self.parent.r, self.r, self.r))
        self.S = np.zeros((self.r, self.r))

class ExternalNode(Node):
    def __init__(self, _parent: Union[Root, 'InternalNode'], _grid: GridParms, _data):
        super().__init__(_data)
        self.parent = _parent
        self.grid = _grid
        self.X = np.zeros((self.parent.r, self.grid.dx))

class Tree:
    def __init__(self, _root: Root):
        self.root = _root
    
    def printTree(self):
        self.root.printTree()