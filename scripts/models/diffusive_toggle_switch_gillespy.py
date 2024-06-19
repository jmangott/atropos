"""Gillespy2 model for the lambd0 phage example."""
import gillespy2

model = gillespy2.Model(name="diffusive_toggle_switch")

b = gillespy2.Parameter(name='b', expression=0.4)
c = gillespy2.Parameter(name='c', expression=0.05)
d0 = gillespy2.Parameter(name='d0', expression=0.01)
d1 = gillespy2.Parameter(name='d1', expression=0.01)

model.add_parameter([b, c, d0, d1])

S0 = gillespy2.Species(name='S0', initial_value=0)
S1 = gillespy2.Species(name='S1', initial_value=0)
S2 = gillespy2.Species(name='S2', initial_value=0)
S3 = gillespy2.Species(name='S3', initial_value=0)
S4 = gillespy2.Species(name='S4', initial_value=0)
S5 = gillespy2.Species(name='S5', initial_value=0)
S6 = gillespy2.Species(name='S6', initial_value=0)
S7 = gillespy2.Species(name='S7', initial_value=0)

model.add_species([S0, S1, S2, S3, S4, S5, S6, S7])

# Compartment 0
r0 = gillespy2.Reaction(name="r0", reactants={S0:1}, products={},
                               propensity_function="c * S0")
r1 = gillespy2.Reaction(name="r1", reactants={S1:1}, products={},
                               propensity_function="c * S1")
r2 = gillespy2.Reaction(name="r2", reactants={}, products={S0:1},
                               propensity_function="b / (b + S1)")
r3 = gillespy2.Reaction(name="r3", reactants={}, products={S1:1},
                               propensity_function="b / (b + S0)")
r4 = gillespy2.Reaction(name="r4", reactants={S0:1}, products={S2:1},
                               propensity_function="d0 * S0")
r5 = gillespy2.Reaction(name="r5", reactants={S1:1}, products={S3:1},
                               propensity_function="d1 * S1")

# Compartment 1
r6 = gillespy2.Reaction(name="r6", reactants={S2:1}, products={},
                               propensity_function="c * S2")
r7 = gillespy2.Reaction(name="r7", reactants={S3:1}, products={},
                               propensity_function="c * S3")
r8 = gillespy2.Reaction(name="r8", reactants={}, products={S2:1},
                               propensity_function="b / (b + S3)")
r9 = gillespy2.Reaction(name="r9", reactants={}, products={S3:1},
                               propensity_function="b / (b + S2)")
r10 = gillespy2.Reaction(name="r10", reactants={S2:1}, products={S0:1},
                               propensity_function="d0 * S2")
r11 = gillespy2.Reaction(name="r11", reactants={S3:1}, products={S1:1},
                               propensity_function="d1 * S3")
r12 = gillespy2.Reaction(name="r12", reactants={S2:1}, products={S4:1},
                               propensity_function="d0 * S2")
r13 = gillespy2.Reaction(name="r13", reactants={S3:1}, products={S5:1},
                               propensity_function="d1 * S3")

# Compartment 2
r14 = gillespy2.Reaction(name="r14", reactants={S4:1}, products={},
                               propensity_function="c * S4")
r15 = gillespy2.Reaction(name="r15", reactants={S5:1}, products={},
                               propensity_function="c * S5")
r16 = gillespy2.Reaction(name="r16", reactants={}, products={S4:1},
                               propensity_function="b / (b + S5)")
r17 = gillespy2.Reaction(name="r17", reactants={}, products={S5:1},
                               propensity_function="b / (b + S4)")
r18 = gillespy2.Reaction(name="r18", reactants={S4:1}, products={S2:1},
                               propensity_function="d0 * S4")
r19 = gillespy2.Reaction(name="r19", reactants={S5:1}, products={S3:1},
                               propensity_function="d1 * S5")
r20 = gillespy2.Reaction(name="r20", reactants={S4:1}, products={S6:1},
                               propensity_function="d0 * S4")
r21 = gillespy2.Reaction(name="r21", reactants={S5:1}, products={S7:1},
                               propensity_function="d1 * S5")

# Compartment 3
r22 = gillespy2.Reaction(name="r22", reactants={S6:1}, products={},
                               propensity_function="c * S6")
r23 = gillespy2.Reaction(name="r23", reactants={S7:1}, products={},
                               propensity_function="c * S7")
r24 = gillespy2.Reaction(name="r24", reactants={}, products={S6:1},
                               propensity_function="b / (b + S7)")
r25 = gillespy2.Reaction(name="r25", reactants={}, products={S7:1},
                               propensity_function="b / (b + S6)")
r26 = gillespy2.Reaction(name="r26", reactants={S6:1}, products={S4:1},
                               propensity_function="d0 * S6")
r27 = gillespy2.Reaction(name="r27", reactants={S7:1}, products={S5:1},
                               propensity_function="d1 * S7")

model.add_reaction([r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18, r19, r20, r21, r22, r23, r24, r25, r26, r27])

tspan = gillespy2.TimeSpan.linspace(500, 2)
model.timespan(tspan)