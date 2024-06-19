"""Gillespy2 model for the lambda phage example."""
import gillespy2

model = gillespy2.Model(name="lambda_phage")

a1 = gillespy2.Parameter(name='a1', expression=0.5)
a2 = gillespy2.Parameter(name='a2', expression=1.0)
a3 = gillespy2.Parameter(name='a3', expression=0.15)
a4 = gillespy2.Parameter(name='a4', expression=0.3)
a5 = gillespy2.Parameter(name='a5', expression=0.3)

b1 = gillespy2.Parameter(name='b1', expression=0.12)
b2 = gillespy2.Parameter(name='b2', expression=0.6)
b3 = gillespy2.Parameter(name='b3', expression=1.0)
b4 = gillespy2.Parameter(name='b4', expression=1.0)
b5 = gillespy2.Parameter(name='b5', expression=1.0)

c1 = gillespy2.Parameter(name='c1', expression=0.0025)
c2 = gillespy2.Parameter(name='c2', expression=0.0007)
c3 = gillespy2.Parameter(name='c3', expression=0.0231)
c4 = gillespy2.Parameter(name='c4', expression=0.01)
c5 = gillespy2.Parameter(name='c5', expression=0.01)

model.add_parameter([a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, c1, c2, c3, c4, c5])

S1 = gillespy2.Species(name='S1', initial_value=0)
S2 = gillespy2.Species(name='S2', initial_value=0)
S3 = gillespy2.Species(name='S3', initial_value=0)
S4 = gillespy2.Species(name='S4', initial_value=0)
S5 = gillespy2.Species(name='S5', initial_value=0)

model.add_species([S1, S2, S3, S4, S5])

reaction1 = gillespy2.Reaction(name="r1", reactants={}, products={S1:1},
                               propensity_function="a1 * b1 / (b1 + S2)")
reaction2 = gillespy2.Reaction(name="r2", reactants={}, products={S2:1},
                               propensity_function="(a2 + S5) * b2 / (b2 + S1)")
reaction3 = gillespy2.Reaction(name="r3", reactants={}, products={S3:1},
                               propensity_function="a3 * b3 * S2 / (b3 * S2 + 1.0)")
reaction4 = gillespy2.Reaction(name="r4", reactants={}, products={S4:1},
                               propensity_function="a4 * b4 * S3 / (b4 * S3 + 1.0)")
reaction5 = gillespy2.Reaction(name="r5", reactants={}, products={S5:1},
                               propensity_function="a5 * b5 * S3 / (b5 * S3 + 1.0)")
reaction6 = gillespy2.Reaction(name="r6", reactants={S1:1}, products={},
                               propensity_function="c1 * S1")
reaction7 = gillespy2.Reaction(name="r7", reactants={S2:1}, products={},
                               propensity_function="c2 * S2")
reaction8 = gillespy2.Reaction(name="r8", reactants={S3:1}, products={},
                               propensity_function="c3 * S3")
reaction9 = gillespy2.Reaction(name="r9", reactants={S4:1}, products={},
                               propensity_function="c4 * S4")
reaction10 = gillespy2.Reaction(name="r10", reactants={S5:1}, products={},
                               propensity_function="c5 * S5")

model.add_reaction([reaction1, reaction2, reaction3, reaction4, reaction5, reaction6, reaction7, reaction8, reaction9, reaction10])

tspan = gillespy2.TimeSpan.linspace(10, 11)
model.timespan(tspan)