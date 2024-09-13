"""
Contains the `Reaction` class for a single reaction of the reaction network and the
`ReactionSystem` class for a reaction system consisting of multiple reactions.
"""

import numpy as np


class Reaction:
    def __init__(self, _propensity: dict, _nu: np.ndarray[int]):
        self.propensity = dict(sorted(_propensity.items()))
        self.nu = _nu


class ReactionSystem:
    def __init__(self, _reactions: list[Reaction], _species_names: list[str]):
        for reaction in _reactions:
            if reaction.nu.size != len(_species_names):
                raise ValueError("Dimensions of `nu` and `_species_names` do not match")

        self.reactions = _reactions
        self.species_names = _species_names

    def size(self):
        return len(self.reactions)

    def d(self):
        return len(self.species_names)
