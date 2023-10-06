"""Contains the `GridParms` class for storing grid parameters for a partition of the reaction network."""
import numpy as np
import numpy.typing as npt

class GridParms:
    def __init__(self, _n: npt.NDArray[np.int_], _binsize: npt.NDArray[np.int_], _liml: npt.NDArray[np.float_], _dep: npt.NDArray[np.bool_]):
        if np.any(np.not_equal([_binsize.size, _liml.size, _dep.shape[1]], _n.size)):
            raise ValueError("Input arrays must be of equal length")
        
        if np.any(_n == 0):
            raise ValueError("Grid size `n` must be larger than 0")

        self.n = _n
        self.binsize = _binsize
        self.liml = _liml
        self.dep = _dep
        self.dx = np.prod(self.n)
        self.h_mult = np.prod(self.binsize)
        self.d = self.n.size
    
    def __str__(self):
        return ", ".join([str(ele) for ele in self.n])