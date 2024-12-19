import numpy as np
import sympy as sp

from atropy.src.generator import Model, Partitioning, run, species

"""
Define variables
"""
S0, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10 = species(
    "S0, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10"
)

"""
Generate model
"""
model = Model((S0, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10))

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


### Full reaction generation:

# model.add_reaction(2*S0, S1, {S0: 0.5 * kA0 * S0 * (S0 - 1.0)})
# model.add_reaction(S1, 2*S0, {S1: kA1 * S1})
# model.add_reaction(S0 + S1, S2, {S0: np.sqrt(kA2) * S0, S1: np.sqrt(kA2) * S1})
# model.add_reaction(S2, S0 + S1, {S2: kA3 * S2})
# model.add_reaction(S0 + S2, S3, {S0: np.sqrt(kA4) * S0, S2: np.sqrt(kA4) * S2})
# model.add_reaction(S3, S0 + S2, {S3: kA5 * S3})
# model.add_reaction(S0 + S3, S4, {S0: np.sqrt(kA6) * S0, S3: np.sqrt(kA6) * S3})
# model.add_reaction(S4, S0 + S3, {S4: kA7 * S4})
# model.add_reaction(S0 + S4, S5, {S0: np.sqrt(kA8) * S0, S4: np.sqrt(kA8) * S4})
# model.add_reaction(S5, S0 + S4, {S5: kA9 * S5})
# model.add_reaction(S3 + S9, S6, {S3: np.sqrt(kA10) * S3, S9: np.sqrt(kA10) * S9})
# model.add_reaction(S6, S3 + S9, {S6: kA11 * S6})
# model.add_reaction(S6, S3 + S10, {S6: kA12 * S6})
# model.add_reaction(S4 + S9, S7, {S4: np.sqrt(kA13) * S4, S9: np.sqrt(kA13) * S9})
# model.add_reaction(S7, S4 + S9, {S7: kA14 * S7})
# model.add_reaction(S7, S4 + S10, {S7: kA15 * S7})
# model.add_reaction(S5 + S9, S8, {S5: np.sqrt(kA16) * S5, S9: np.sqrt(kA16) * S9})
# model.add_reaction(S8, S5 + S9, {S8: kA17 * S8})
# model.add_reaction(S8, S5 + S10, {S8: kA18 * S8})


### Shorter reaction generation:

model.add_reaction(2 * S0, S1, kA0)
model.add_reaction(S1, 2 * S0, kA1)
model.add_reaction(S0 + S1, S2, kA2)
model.add_reaction(S2, S0 + S1, kA3)
model.add_reaction(S0 + S2, S3, kA4)
model.add_reaction(S3, S0 + S2, kA5)
model.add_reaction(S0 + S3, S4, kA6)
model.add_reaction(S4, S0 + S3, kA7)
model.add_reaction(S0 + S4, S5, kA8)
model.add_reaction(S5, S0 + S4, kA9)
model.add_reaction(S3 + S9, S6, kA10)
model.add_reaction(S6, S3 + S9, kA11)
model.add_reaction(S6, S3 + S10, kA12)
model.add_reaction(S4 + S9, S7, kA13)
model.add_reaction(S7, S4 + S9, kA14)
model.add_reaction(S7, S4 + S10, kA15)
model.add_reaction(S5 + S9, S8, kA16)
model.add_reaction(S8, S5 + S9, kA17)
model.add_reaction(S8, S5 + S10, kA18)

model.generate_reaction_system()


"""
Generate tree and initial condition
"""
r = np.array([5, 4, 3])
p0 = "(S0 S1 S2)(((S3 S4 S6 S7)(S5 S8))(S9 S10))"
partitioning = Partitioning(p0, r, model)

n = np.array([46, 16, 16, 11, 11, 11, 4, 4, 4, 56, 56])
d = n.size
binsize = np.ones(d, dtype=int)
liml = np.zeros(d)
partitioning.add_grid_params(n, binsize, liml)

partitioning.generate_tree()

n_basisfunctions = np.ones(r.size, dtype="int")

partitioning.generate_initial_condition(n_basisfunctions)

C = 0.2
Cinv = 1 / C
mu = np.array([40, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0])
mu_perm = mu[partitioning.tree.species]

polynomials_dict = {
    S0: sp.exp(-0.5 * Cinv * (S0 - mu_perm[0]) ** 2),
    S1: sp.exp(-0.5 * Cinv * (S1 - mu_perm[1]) ** 2),
    S2: sp.exp(-0.5 * Cinv * (S2 - mu_perm[2]) ** 2),
    S3: sp.exp(-0.5 * Cinv * (S3 - mu_perm[3]) ** 2),
    S4: sp.exp(-0.5 * Cinv * (S4 - mu_perm[4]) ** 2),
    S5: sp.exp(-0.5 * Cinv * (S5 - mu_perm[5]) ** 2),
    S6: sp.exp(-0.5 * Cinv * (S6 - mu_perm[6]) ** 2),
    S7: sp.exp(-0.5 * Cinv * (S7 - mu_perm[7]) ** 2),
    S8: sp.exp(-0.5 * Cinv * (S8 - mu_perm[8]) ** 2),
    S9: sp.exp(-0.5 * Cinv * (S9 - mu_perm[9]) ** 2),
    S10: sp.exp(-0.5 * Cinv * (S10 - mu_perm[10]) ** 2),
}
partitioning.set_initial_condition(polynomials_dict)

_, marginal_distribution = partitioning.tree.calculateObservables(
    np.zeros(partitioning.tree.root.grid.d(), dtype="int")
)
norm = np.sum(marginal_distribution[partitioning.tree.species_names[0]])
print("norm:", norm)
partitioning.tree.root.Q[0, 0, 0] /= norm

print(partitioning.tree)


"""
write input file and run
"""
run(partitioning, "bax_Stefan", 1e-3, 1, snapshot=10, method="implicit_Euler")
