import numpy as np

from scripts.reaction_class import Reaction, ReactionSystem

b = 0.4
c = 0.05
d0 = 0.01
d1 = 0.01

reactions = [None] * (8 * 4 - 4)

# Compartment 0
reactions[0] = Reaction({0: lambda x: c * x}, [-1, 0, 0, 0, 0, 0, 0, 0])
reactions[1] = Reaction({1: lambda x: c * x}, [0, -1, 0, 0, 0, 0, 0, 0])
reactions[2] = Reaction({1: lambda x: b / (b + x)}, [1, 0, 0, 0, 0, 0, 0, 0])
reactions[3] = Reaction({0: lambda x: b / (b + x)}, [0, 1, 0, 0, 0, 0, 0, 0])
reactions[4] = Reaction({0: lambda x: d0 * x}, [-1, 0, 1, 0, 0, 0, 0, 0])
reactions[5] = Reaction({1: lambda x: d1 * x}, [0, -1, 0, 1, 0, 0, 0, 0])

# Compartment 1
reactions[6] = Reaction({2: lambda x: c * x}, [0, 0, -1, 0, 0, 0, 0, 0])
reactions[7] = Reaction({3: lambda x: c * x}, [0, 0, 0, -1, 0, 0, 0, 0])
reactions[8] = Reaction({3: lambda x: b / (b + x)}, [0, 0, 1, 0, 0, 0, 0, 0])
reactions[9] = Reaction({2: lambda x: b / (b + x)}, [0, 0, 0, 1, 0, 0, 0, 0])
reactions[10] = Reaction({2: lambda x: d0 * x}, [1, 0, -1, 0, 0, 0, 0, 0])
reactions[11] = Reaction({3: lambda x: d1 * x}, [0, 1, 0, -1, 0, 0, 0, 0])
reactions[12] = Reaction({2: lambda x: d0 * x}, [0, 0, -1, 0, 1, 0, 0, 0])
reactions[13] = Reaction({3: lambda x: d1 * x}, [0, 0, 0, -1, 0, 1, 0, 0])

# Compartment 2
reactions[14] = Reaction({4: lambda x: c * x}, [0, 0, 0, 0, -1, 0, 0, 0])
reactions[15] = Reaction({5: lambda x: c * x}, [0, 0, 0, 0, 0, -1, 0, 0])
reactions[16] = Reaction({5: lambda x: b / (b + x)}, [0, 0, 0, 0, 1, 0, 0, 0])
reactions[17] = Reaction({4: lambda x: b / (b + x)}, [0, 0, 0, 0, 0, 1, 0, 0])
reactions[18] = Reaction({4: lambda x: d0 * x}, [0, 0, 1, 0, -1, 0, 0, 0])
reactions[19] = Reaction({5: lambda x: d1 * x}, [0, 0, 0, 1, 0, -1, 0, 0])
reactions[20] = Reaction({4: lambda x: d0 * x}, [0, 0, 0, 0, -1, 0, 1, 0])
reactions[21] = Reaction({5: lambda x: d1 * x}, [0, 0, 0, 0, 0, -1, 0, 1])

# Compartment 3
reactions[22] = Reaction({6: lambda x: c * x}, [0, 0, 0, 0, 0, 0, -1, 0])
reactions[23] = Reaction({7: lambda x: c * x}, [0, 0, 0, 0, 0, 0, 0, -1])
reactions[24] = Reaction({7: lambda x: b / (b + x)}, [0, 0, 0, 0, 0, 0, 1, 0])
reactions[25] = Reaction({6: lambda x: b / (b + x)}, [0, 0, 0, 0, 0, 0, 0, 1])
reactions[26] = Reaction({6: lambda x: d0 * x}, [0, 0, 0, 0, 1, 0, -1, 0])
reactions[27] = Reaction({7: lambda x: d1 * x}, [0, 0, 0, 0, 0, 1, 0, -1])

species_names = ["A0", "B0", "A1", "B1", "A2", "B2", "A3", "B3"]
reaction_system = ReactionSystem(reactions, species_names)