"""Gillespy2 model for the cascade example."""
import gillespy2

model = gillespy2.Model(name="cascade")

d = 20

k0 = gillespy2.Parameter(name='k0', expression=0.7)
kp = gillespy2.Parameter(name='kp', expression=5.0)
km = gillespy2.Parameter(name='km', expression=0.07)

model.add_parameter([k0, kp, km])

species = [gillespy2.Species(name='S{:02d}'.format(i), initial_value=0) for i in range(d)]
model.add_species(species)

reactions = []

# Creation
reactions.append(gillespy2.Reaction(name="r00", reactants={}, products={species[0]: 1},
                               propensity_function="k0"))

for i in range(1, d):
    reactions.append(gillespy2.Reaction(name="r{:02d}".format(i), reactants={}, products={species[i]: 1}, propensity_function="S{idx:02d} / (kp + S{idx:02d})".format(idx=i-1)))

# Annihilation
for i in range(d):
    reactions.append(gillespy2.Reaction(name="r{:02d}".format(i + d), reactants={species[i]: 1}, products={}, propensity_function="km * S{:02d}".format(i)))

model.add_reaction(reactions)

tspan = gillespy2.TimeSpan.linspace(350.0, 2)
model.timespan(tspan)