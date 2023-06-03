import matplotlib.pyplot as plt
import numpy as np
from pysb import *
from pysb.macros import *
from pysb.integrate import odesolve

Model()

Monomer('S1')
Monomer('S2')

Parameter('b', 0.4)
Parameter('c', 0.05)

Observable('x1', S1)
Observable('x2', S2)

Expression('alpha1', c)
Expression('alpha2', c)
Expression('alpha3', b / (b + x2))
Expression('alpha4', b / (b + x1))

Rule('rule1', S1() >> None, alpha1)
Rule('rule2', S2() >> None, alpha2)
Rule('rule3', None >> S1(), alpha3)
Rule('rule4', None >> S2(), alpha4)

Parameter('x1_0', 1.0)
Parameter('x2_0', 0.0)

Initial(S1(), x1_0)
Initial(S2(), x2_0)

if __name__ == '__main__':
    tspan = np.linspace(0, 500, 100)
    y = odesolve(model, tspan)

    plt.plot(tspan, y['x1'], label="$S_1$")
    plt.plot(tspan, y['x2'], label="$S_2$")
    plt.xlabel('time')
    plt.ylabel('population')
    plt.legend()
    plt.show()