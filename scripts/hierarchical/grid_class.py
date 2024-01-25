"""Contains the `GridParms` class for storing grid parameters for a partition of the reaction network."""
import numpy as np
import numpy.typing as npt

from scripts.hierarchical.reaction_class import ReactionSystem

class GridParms:
    def __init__(self, _n: npt.NDArray[np.int_], _binsize: npt.NDArray[np.int_], _liml: npt.NDArray[np.float_], _species: npt.NDArray[np.int_] = None, _dep: npt.NDArray[np.bool_] = None, _nu: npt.NDArray[np.int_] = None):
        if np.any(np.not_equal([_binsize.size, _liml.size], _n.size)):
            raise Exception("Input arrays must be of equal length")
        
        if np.any(_n == 0):
            raise Exception("Grid size `n` must be larger than 0")

        self.n = _n
        self.binsize = _binsize
        self.liml = _liml
        self.species = _species

        if _species is None:
            self.species = np.arange(self.d(), dtype="int")
        else:
            if _species.size != _n.size:
                raise Exception("Input arrays must be of equal length")
            self.species = _species

        # These values have to be set according to a given reaction system,
        # this is done with the `initialize` method
        self.dep = _dep
        self.nu = _nu
    
    def dx(self):
        return np.prod(self.n)

    def h_mult(self):
        return np.prod(self.binsize)
    
    def d(self):
        return self.n.size
    
    def n_reactions(self):
        return self.dep.shape[1]

    def initialize(self, reaction_system: ReactionSystem):
        for reaction in reaction_system.reactions:
            if (self.d() != len(reaction.nu)):
                raise Exception("`d` and `len(nu)` must be equal")

        self.dep = np.zeros((self.d(), reaction_system.size()), dtype="bool")
        self.nu = np.zeros((self.d(), reaction_system.size()), dtype = "int")

        for mu, reaction in enumerate(reaction_system.reactions):
            for i in reaction.propensity.keys():
                self.dep[i, mu] = True
            self.nu[:, mu] = reaction.nu

    def permute(self, permutation: npt.NDArray[np.int_]):
        self.n = self.n[permutation]
        self.binsize = self.binsize[permutation]
        self.liml = self.liml[permutation]
        self.species = self.species[permutation]
        self.dep = self.dep[permutation, :]
        self.nu = self.nu[permutation, :]

    def __str__(self):
        return str(self.n)