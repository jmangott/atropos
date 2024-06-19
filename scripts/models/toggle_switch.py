import numpy as np

from scripts.reaction_class import Reaction, ReactionSystem

b = 0.4
c = 0.05

reactions = [None] * 4

reactions[0] = Reaction({0: lambda x: c * x}, [-1, 0])
reactions[1] = Reaction({1: lambda x: c * x}, [0, -1])
reactions[2] = Reaction({1: lambda x: b / (b + x)}, [1, 0])
reactions[3] = Reaction({0: lambda x: b / (b + x)}, [0, 1])

species_names = ["A", "B"]
reaction_system = ReactionSystem(reactions, species_names)