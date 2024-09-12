import numpy as np
from src.reaction import Reaction, ReactionSystem

d = 20

k0 = 0.7
kp = 5.0
km = 0.07

reactions = [None] * 2 * d

# Creation
nu_c0 = np.zeros(d, dtype="int")
nu_c0[0] = 1
reactions[0] = Reaction({0: lambda x: k0}, nu_c0)

for i in range(1, d):
    nu_c = np.zeros(d, dtype="int")
    nu_c[i] = 1
    reactions[i] = Reaction({(i - 1): lambda x: x / (kp + x)}, nu_c)

# Annihilation
for i in range(d):
    nu_a = np.zeros(d, dtype="int")
    nu_a[i] = -1
    reactions[i + d] = Reaction({i: lambda x: km * x}, nu_a)

species_names = ["S{}".format(i) for i in range(d)]

reaction_system = ReactionSystem(reactions, species_names)
