"""Contains the `GridParms` class for storing grid parameters for a partition of the reaction network."""
import numpy as np
import numpy.typing as npt

from reaction_class import *

class GridParms:
    def __init__(self, _n: npt.NDArray[np.int_], _binsize: npt.NDArray[np.int_], _liml: npt.NDArray[np.float_], _dep: npt.NDArray[np.bool_] = None, _nu: npt.NDArray[np.int_] = None):
        if np.any(np.not_equal([_binsize.size, _liml.size], _n.size)):
            raise ValueError("Input arrays must be of equal length")
        
        if np.any(_n == 0):
            raise ValueError("Grid size `n` must be larger than 0")

        self.n = _n
        self.binsize = _binsize
        self.liml = _liml
        self.dx = np.prod(self.n)
        self.h_mult = np.prod(self.binsize)
        self.d = self.n.size

        # These values have to be set according to a given reaction system,
        # this is done with the `Initialize` method
        self.dep = _dep
        self.nu = _nu
        if self.dep is None:
            self.n_reactions = None
        else:
            self.n_reactions = self.dep.shape[0] 

    def initialize(self, reaction_system: ReactionSystem):
        for reaction in reaction_system.reactions:
            if (self.d != len(reaction.nu)):
                raise ValueError("`d` and `len(nu)` must be equal")

        self.n_reactions = reaction_system.size()
        self.dep = np.zeros((self.d, self.n_reactions), dtype="bool")
        self.nu = np.zeros((self.d, self.n_reactions), dtype = "int")

        for mu, reaction in enumerate(reaction_system.reactions):
            for i in reaction.propensity.keys():
                self.dep[i, mu] = True
            self.nu[:, mu] = reaction.nu
    
    def __str__(self):
        return "[" + ", ".join([str(ele) for ele in self.n]) + "]"