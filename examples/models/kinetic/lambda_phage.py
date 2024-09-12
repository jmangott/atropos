import numpy as np
from src.reaction import Reaction, ReactionSystem

kA0 = 0.5
kA1 = 1.0
kA2 = 0.15
kA3 = 0.3
kA4 = 0.3
kB0 = 0.12
kB1 = 0.6
kB2 = 1.0
kB3 = 1.0
kB4 = 1.0
kC0 = 0.0025
kC1 = 0.0007
kC2 = 0.0231
kC3 = 0.01
kC4 = 0.01

reactions = [None] * 10

reactions[0] = Reaction(
    {1: lambda x: kA0 * kB0 / (kB0 + x)}, np.array([1, 0, 0, 0, 0], dtype="int")
)

reactions[1] = Reaction(
    {0: lambda x: kB1 / (kB1 + x), 4: lambda x: kA1 + x},
    np.array([0, 1, 0, 0, 0], dtype="int"),
)

reactions[2] = Reaction(
    {1: lambda x: kA2 * kB2 * x / (kB2 * x + 1.0)},
    np.array([0, 0, 1, 0, 0], dtype="int"),
)

reactions[3] = Reaction(
    {2: lambda x: kA3 * kB3 * x / (kB3 * x + 1.0)},
    np.array([0, 0, 0, 1, 0], dtype="int"),
)

reactions[4] = Reaction(
    {2: lambda x: kA4 * kB4 * x / (kB4 * x + 1.0)},
    np.array([0, 0, 0, 0, 1], dtype="int"),
)

reactions[5] = Reaction({0: lambda x: kC0 * x}, np.array([-1, 0, 0, 0, 0], dtype="int"))

reactions[6] = Reaction({1: lambda x: kC1 * x}, np.array([0, -1, 0, 0, 0], dtype="int"))

reactions[7] = Reaction({2: lambda x: kC2 * x}, np.array([0, 0, -1, 0, 0], dtype="int"))

reactions[8] = Reaction({3: lambda x: kC3 * x}, np.array([0, 0, 0, -1, 0], dtype="int"))

reactions[9] = Reaction({4: lambda x: kC4 * x}, np.array([0, 0, 0, 0, -1], dtype="int"))

species_names = ["S1", "S2", "S3", "S4", "S5"]

reaction_system = ReactionSystem(reactions, species_names)
