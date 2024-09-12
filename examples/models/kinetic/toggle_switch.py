import numpy as np
from src.reaction import Reaction, ReactionSystem

b = 0.4
c = 0.05

reactions = [None] * 4

reactions[0] = Reaction({0: lambda x: c * x}, np.array([-1, 0], dtype="int"))
reactions[1] = Reaction({1: lambda x: c * x}, np.array([0, -1], dtype="int"))
reactions[2] = Reaction({1: lambda x: b / (b + x)}, np.array([1, 0], dtype="int"))
reactions[3] = Reaction({0: lambda x: b / (b + x)}, np.array([0, 1], dtype="int"))

species_names = ["A", "B"]
reaction_system = ReactionSystem(reactions, species_names)
