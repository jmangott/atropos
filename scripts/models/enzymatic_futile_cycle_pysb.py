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
Monomer('S6')

Parameter('kp1', 0.4)
Parameter('kp2', 100.0)
Parameter('kp3', 100.0)
Parameter('km1', 2.0)
Parameter('km2', 1.0)
Parameter('km3', 50.0)

Observable('x1', S1)
Observable('x2', S2)
Observable('x3', S3)
Observable('x4', S4)
Observable('x5', S5)
Observable('x6', S6)

Rule('rule1', S1() + S3() >> S5(), kp1)
Rule('rule2', S5() >> S1() + S3(), kp2)
Rule('rule3', S2() + S4() >> S6(), km1)
Rule('rule4', S6() >> S2() + S4(), km2)
Rule('rule5', S5() >> S2() + S3(), kp3)
Rule('rule6', S6() >> S1() + S4(), km3)

Parameter('x1_0', 30.0)
Parameter('x2_0', 90.0)
Parameter('x3_0', 2.0)
Parameter('x4_0', 2.0)
Parameter('x5_0', 0.0)
Parameter('x6_0', 0.0)

Initial(S1(), x1_0)
Initial(S2(), x2_0)
Initial(S3(), x3_0)
Initial(S4(), x4_0)
Initial(S5(), x5_0)
Initial(S6(), x6_0)

if __name__ == '__main__':
    tspan = np.linspace(0, 10.0, 100)
    y = odesolve(model, tspan)

    plt.plot(tspan, y['x1'], label="$S_1$")
    plt.plot(tspan, y['x2'], label="$S_2$")
    plt.plot(tspan, y['x3'], label="$S_3$")
    plt.plot(tspan, y['x4'], label="$S_4$")
    plt.plot(tspan, y['x5'], label="$S_5$")
    plt.plot(tspan, y['x6'], label="$S_5$")
    plt.xlabel('time')
    plt.ylabel('population')
    plt.legend()
    plt.show()