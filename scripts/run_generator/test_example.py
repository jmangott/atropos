from scripts.generator_class import Model, Partitioning, run, species
import numpy as np
import sympy as sp



"""
Define some symbols, to then try out with those symbols
"""
NF, GR, O, H = species("NF, GR, O, H")



"""
Try out model
"""
model = Model((NF, GR, O, H))
print(model.species, model.reactions)

model.add_reaction(3*NF + 7*GR, 2*H + O, {NF: NF**2, H: 1/(1 + H**2)})
print(model.reactions)

model.add_reactions([7*H + GR, 2*H],[NF, 3*O + 2*NF],[{H: H}, {O: 3*O}])
print(model.reactions)

model.generate_reaction_system()
print(model.reaction_system)



"""
Try out partitioning
"""
r = np.array([5,5])
partitioning = Partitioning('((NF GR)(O))(H)', r, model)
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

partitioning.set_initial_condition({NF: 3*NF, GR: 2*GR, O: O, H: 3*H})


"""
Try out run
"""
run(partitioning, 'test_example', 1e-3, 1, method = "implicit_Euler")