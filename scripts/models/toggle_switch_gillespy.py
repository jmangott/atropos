"""Gillespy2 model for the lambda phage example."""
import gillespy2

model = gillespy2.Model(name="toggle_switch")

b = gillespy2.Parameter(name='b', expression=0.4)
c = gillespy2.Parameter(name='c', expression=0.05)

model.add_parameter([b, c])

A = gillespy2.Species(name='A', initial_value=0)
B = gillespy2.Species(name='B', initial_value=0)

model.add_species([A, B])

r0 = gillespy2.Reaction(name="r0", reactants={A:1}, products={},
                               propensity_function="c * A")
r1 = gillespy2.Reaction(name="r1", reactants={B:1}, products={},
                               propensity_function="c * B")
r2 = gillespy2.Reaction(name="r2", reactants={}, products={A:1},
                               propensity_function="b / (b + B)")
r3 = gillespy2.Reaction(name="r3", reactants={}, products={B:1},
                               propensity_function="b / (b + A)")

model.add_reaction([r0, r1, r2, r3])

tspan = gillespy2.TimeSpan.linspace(500, 2)
model.timespan(tspan)