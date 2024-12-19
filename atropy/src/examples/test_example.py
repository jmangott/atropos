import numpy as np

from atropy.src.generator import Model, Partitioning, run, species

"""
Define some symbols, to then try out with those symbols
"""
NF, GR, Oo, H = species("NF, GR, Oo, H")


"""
Try out model
"""
model = Model((NF, GR, Oo, H))
print(model.species, model.reactions)

model.add_reaction(3 * NF + 7 * GR, 2 * H + Oo, {NF: NF**2, H: 1 / (1 + H**2)})
print(model.reactions)

model.add_reactions([7 * H + GR, 2 * H], [NF, 3 * Oo + 2 * NF], [{H: H}, {Oo: 3 * Oo}])
print(model.reactions)

model.generate_reaction_system()
print(model.reaction_system)


"""
Try out partitioning
"""
r = np.array([5, 5])
partitioning = Partitioning("((NF GR)(Oo))(H)", r, model)
print(partitioning.partition)

n = np.array([16, 11, 11, 11])
d = n.size
binsize = np.ones(d, dtype=int)
liml = np.zeros(d)
partitioning.add_grid_params(n, binsize, liml)
print(partitioning.grid)

partitioning.generate_tree()
print(partitioning.tree)

n_basisfunctions = np.ones(r.size, dtype="int")
partitioning.generate_initial_condition(n_basisfunctions)
print(partitioning.initial_conditions)

partitioning.set_initial_condition({NF: 3 * NF, GR: 2 * GR, Oo: Oo, H: 3 * H})


"""
Try out run
"""
run(partitioning, "test_example", 1e-3, 1, method="implicit_Euler")
