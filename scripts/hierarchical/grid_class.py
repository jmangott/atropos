import numpy as np
import numpy.typing as npt

class GridParms:
    def __init__(self, _n: npt.NDArray[np.int_], _binsize: npt.NDArray[np.int_], _liml: npt.NDArray[np.float_]):
        self.n = _n
        self.binsize = _binsize
        self.liml = _liml
        self.dx = np.prod(self.n)
        self.h_mult = np.prod(self.binsize)
        self.d = self.n.size