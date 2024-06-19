import numpy as np
from scripts.reaction_class import Reaction, ReactionSystem

kp1 = 0.4
kp2 = 100.0
kp3 = 100.0
km1 = 2.0
km2 = 1.0
km3 = 50.0

reactions = [None] * 6

reactions[0] = Reaction({0: lambda x: np.sqrt(kp1) * x, 1: lambda x: np.sqrt(kp1) * x}, [-1, -1, 1, 0, 0, 0])
reactions[1] = Reaction({2: lambda x: kp2 * x}, [1, 1, -1, 0, 0, 0])
reactions[2] = Reaction({3: lambda x: np.sqrt(km1) * x, 4: lambda x: np.sqrt(km1) * x}, [0, 0, 0, -1, -1, 1])
reactions[3] = Reaction({5: lambda x: km2 * x}, [0, 0, 0, 1, 1, -1])
reactions[4] = Reaction({2: lambda x: kp3 * x}, [0, 1, -1, 1, 0, 0])
reactions[5] = Reaction({5: lambda x: km3 * x}, [1, 0, 0, 0, 1, -1])

species_names = ["S0", "S1", "S2", "S3", "S4", "S5"]

reaction_system = ReactionSystem(reactions, species_names)