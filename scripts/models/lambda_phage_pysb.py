"""PySB model for the lambda phage example."""
import matplotlib.pyplot as plt
import numpy as np
from pysb import *
from pysb.macros import *
from pysb.integrate import odesolve

Model()

Monomer('S1')
Monomer('S2')
Monomer('S3')
Monomer('S4')
Monomer('S5')

Parameter('a1', 0.5)
Parameter('a2', 1.0)
Parameter('a3', 0.15)
Parameter('a4', 0.3)
Parameter('a5', 0.3)

Parameter('b1', 0.12)
Parameter('b2', 0.6)
Parameter('b3', 1.0)
Parameter('b4', 1.0)
Parameter('b5', 1.0)

Parameter('c1', 0.0025)
Parameter('c2', 0.0007)
Parameter('c3', 0.0231)
Parameter('c4', 0.01)
Parameter('c5', 0.01)

Observable('x1', S1)
Observable('x2', S2)
Observable('x3', S3)
Observable('x4', S4)
Observable('x5', S5)

Expression('alpha1', a1 * b1 / (b1 + x2))
Expression('alpha2', (a2 + x5) * b2 / (b2 + x1))
Expression('alpha3', a3 * b3 * x2 / (b3 * x2 + 1.0))
Expression('alpha4', a4 * b4 * x3 / (b4 * x3 + 1.0))
Expression('alpha5', a5 * b5 * x3 / (b5 * x3 + 1.0))
Expression('alpha6', c1)
Expression('alpha7', c2)
Expression('alpha8', c3)
Expression('alpha9', c4)
Expression('alpha10', c5)

Rule('rule1', None >> S1(), alpha1)
Rule('rule2', None >> S2(), alpha2)
Rule('rule3', None >> S3(), alpha3)
Rule('rule4', None >> S4(), alpha4)
Rule('rule5', None >> S5(), alpha5)
Rule('rule6', S1() >> None, alpha6)
Rule('rule7', S2() >> None, alpha7)
Rule('rule8', S3() >> None, alpha8)
Rule('rule9', S4() >> None, alpha9)
Rule('rule10', S5() >> None, alpha10)

Parameter('x1_0', 0.15)
Parameter('x2_0', 0.15)
Parameter('x3_0', 0.15)
Parameter('x4_0', 0.15)
Parameter('x5_0', 0.15)

Initial(S1(), x1_0)
Initial(S2(), x2_0)
Initial(S3(), x3_0)
Initial(S4(), x4_0)
Initial(S5(), x5_0)

if __name__ == '__main__':
    tspan = np.linspace(0, 10, 100)
    y = odesolve(model, tspan)

    plt.plot(tspan, y['x1'], label="$S_1$")
    plt.plot(tspan, y['x2'], label="$S_2$")
    plt.plot(tspan, y['x3'], label="$S_3$")
    plt.plot(tspan, y['x4'], label="$S_4$")
    plt.plot(tspan, y['x5'], label="$S_5$")
    plt.xlabel('time')
    plt.ylabel('population')
    plt.legend()
    plt.show()