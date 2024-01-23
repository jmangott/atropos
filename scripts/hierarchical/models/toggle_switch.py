from scripts.hierarchical.reaction_class import Reaction, ReactionSystem

kB = 0.4
kC = 0.05

reactions = [None] * 4

reactions[0] = Reaction({0: lambda x: kC * x}, [-1, 0])

reactions[1] = Reaction({1: lambda x: kC * x}, [0, -1])

reactions[2] = Reaction({1: lambda x: kB / (kB + x)}, [1, 0])

reactions[3] = Reaction({0: lambda x: kB / (kB + x)}, [0, 1])

species_names = ["S1", "S2"]

reaction_system = ReactionSystem(reactions, species_names)