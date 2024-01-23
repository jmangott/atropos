import numpy as np
from scripts.hierarchical.reaction_class import Reaction, ReactionSystem

kA0 = 0.0002
kA1 = 0.001
kA2 = 0.0002
kA3 = 0.001
kA4 = 0.0002
kA5 = 0.001
kA6 = 0.0002
kA7 = 0.001
kA8 = 0.0002
kA9 = 0.001

kA10 = 3.0e-5
kA11 = 0.001
kA12 = 10.0
kA13 = 3.0e-5
kA14 = 0.001
kA15 = 10.0
kA16 = 3.0e-5
kA17 = 0.001
kA18 = 10.0

reactions = [None] * 19

reactions[0] = Reaction({0: lambda x: 0.5 * kA0 * x * (x - 1.0)}, [-2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])

reactions[1] = Reaction({1: lambda x: kA1 * x}, [2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0])

reactions[2] = Reaction({0: lambda x: np.sqrt(kA2) * x, 1: lambda x: np.sqrt(kA2) * x}, [-1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0])

reactions[3] = Reaction({2: lambda x: kA3 * x}, [1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0])

reactions[4] = Reaction({0: lambda x: np.sqrt(kA4) * x, 2: lambda x: np.sqrt(kA4) * x }, [-1, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0])

reactions[5] = Reaction({3: lambda x: kA5 * x }, [1, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0])

reactions[6] = Reaction({0: lambda x: np.sqrt(kA6) * x, 3: lambda x: np.sqrt(kA6) * x}, [-1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0])

reactions[7] = Reaction({4: lambda x: kA7 * x}, [1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0])

reactions[8] = Reaction({0: lambda x: np.sqrt(kA8) * x, 4: lambda x: np.sqrt(kA8) * x}, [-1, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0])

reactions[9] = Reaction({5: lambda x: kA9 * x}, [1, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0])

reactions[10] = Reaction({3: lambda x: np.sqrt(kA10) * x, 9: lambda x: np.sqrt(kA10) * x}, [0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0])

reactions[11] = Reaction({6: lambda x: kA11 * x}, [0, 0, 0, 1, 0, 0, -1, 0, 0, 1, 0])

reactions[12] = Reaction({6: lambda x: kA12 * x}, [0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1])

reactions[13] = Reaction({4: lambda x: np.sqrt(kA13) * x, 9: lambda x: np.sqrt(kA13) * x}, [0, 0, 0, 0, -1, 0, 0, 1, 0, -1, 0])

reactions[14] = Reaction({7: lambda x: kA14 * x}, [0, 0, 0, 0, 1, 0, 0, -1, 0, 1, 0])

reactions[15] = Reaction({7: lambda x: kA15 * x}, [0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 1])

reactions[16] = Reaction({5: lambda x: np.sqrt(kA16) * x, 9: lambda x: np.sqrt(kA16) * x}, [0, 0, 0, 0, 0, -1, 0, 0, 1, -1, 0])

reactions[17] = Reaction({8: lambda x: kA17 * x}, [0, 0, 0, 0, 0, 1, 0, 0, -1, 1, 0])

reactions[18] = Reaction({8: lambda x: kA18 * x}, [0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 1])

species_names = ["S" + str(i) for i in range(11)]

reaction_system = ReactionSystem(reactions, species_names)