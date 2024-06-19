"""Gillespy2 model for the enzymatic futile cycle example."""
import gillespy2

model = gillespy2.Model(name="enzymatic_futile_cycle")

kp1 = gillespy2.Parameter(name='kp1', expression=0.4)
kp2 = gillespy2.Parameter(name='kp2', expression=100.0)
kp3 = gillespy2.Parameter(name='kp3', expression=100.0)
km1 = gillespy2.Parameter(name='km1', expression=2.0)
km2 = gillespy2.Parameter(name='km2', expression=1.0)
km3 = gillespy2.Parameter(name='km3', expression=50.0)

model.add_parameter([kp1, kp2, kp3, km1, km2, km3])

S1 = gillespy2.Species(name='S1', initial_value=0)
S2 = gillespy2.Species(name='S2', initial_value=0)
S3 = gillespy2.Species(name='S3', initial_value=0)
S4 = gillespy2.Species(name='S4', initial_value=0)
S5 = gillespy2.Species(name='S5', initial_value=0)
S6 = gillespy2.Species(name='S6', initial_value=0)

model.add_species([S1, S2, S3, S4, S5, S6])

reaction1 = gillespy2.Reaction(name="r1", reactants={S1:1, S3:1}, products={S5:1},
                               propensity_function="kp1 * S1 * S3")
reaction2 = gillespy2.Reaction(name="r2", reactants={S5:1}, products={S1:1, S3:1},
                               propensity_function="kp2 * S5")
reaction3 = gillespy2.Reaction(name="r3", reactants={S2:1, S4:1}, products={S6:1},
                               propensity_function="km1 * S2 * S4")
reaction4 = gillespy2.Reaction(name="r4", reactants={S6:1}, products={S2:1, S4:1},
                               propensity_function="km2 * S6")
reaction5 = gillespy2.Reaction(name="r5", reactants={S5:1}, products={S2:1, S3:1},
                               propensity_function="kp3 * S5")
reaction6 = gillespy2.Reaction(name="r6", reactants={S6:1}, products={S1:1, S4:1},
                               propensity_function="km3 * S6")

model.add_reaction([reaction1, reaction2, reaction3, reaction4, reaction5, reaction6])

tspan = gillespy2.TimeSpan.linspace(10, 11)
model.timespan(tspan)