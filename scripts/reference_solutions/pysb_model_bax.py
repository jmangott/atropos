"""PySB model for the BAX pore assembly example."""
import matplotlib.pyplot as plt
import numpy as np
from pysb import *
from pysb.macros import *
from pysb.integrate import odesolve

max_size = 6

Model()

Monomer('Bax1')
Monomer('Bax2')
Monomer('Bax3')
Monomer('Bax4')
Monomer('Bax5')
Monomer('Bax6')
Monomer('Bax4mSmac')
Monomer('Bax5mSmac')
Monomer('Bax6mSmac')
Monomer('mSmac')
Monomer('cSmac')

Parameter('kA0', 0.0002)
Parameter('kA1', 0.001)
Parameter('kA2', 0.0002)
Parameter('kA3', 0.001)
Parameter('kA4', 0.0002)
Parameter('kA5', 0.001)
Parameter('kA6', 0.0002)
Parameter('kA7', 0.001)
Parameter('kA8', 0.0002)
Parameter('kA9', 0.001)

Parameter('kA10', 3.0e-5)
Parameter('kA11', 0.001)
Parameter('kA12', 10.0)
Parameter('kA13', 3.0e-5)
Parameter('kA14', 0.001)
Parameter('kA15', 10.0)
Parameter('kA16', 3.0e-5)
Parameter('kA17', 0.001)
Parameter('kA18', 10.0)

Observable('oBax1', Bax1)
Observable('oBax2', Bax2)
Observable('oBax3', Bax3)
Observable('oBax4', Bax4)
Observable('oBax5', Bax5)
Observable('oBax6', Bax6)
Observable('oBax4mSmac', Bax4mSmac)
Observable('oBax5mSmac', Bax5mSmac)
Observable('oBax6mSmac', Bax6mSmac)
Observable('omSmac', mSmac)
Observable('ocSmac', cSmac)

Rule('rule0', Bax1() + Bax1() >> Bax2(), kA0)
Rule('rule1', Bax2() >> Bax1() + Bax1(), kA1)

Rule('rule2', Bax1() + Bax2() >> Bax3(), kA2)
Rule('rule3', Bax3() >> Bax1() + Bax2(), kA3)

Rule('rule4', Bax1() + Bax3() >> Bax4(), kA4)
Rule('rule5', Bax4() >> Bax1() + Bax3(), kA5)

Rule('rule6', Bax1() + Bax4() >> Bax5(), kA6)
Rule('rule7', Bax5() >> Bax1() + Bax4(), kA7)

Rule('rule8', Bax1() + Bax5() >> Bax6(), kA8)
Rule('rule9', Bax6() >> Bax1() + Bax5(), kA9)

Rule('rule10', Bax4() + mSmac() >> Bax4mSmac(), kA10)
Rule('rule11', Bax4mSmac() >> Bax4() + mSmac(), kA11)
Rule('rule12', Bax4mSmac() >> Bax4() + cSmac(), kA12)

Rule('rule13', Bax5() + mSmac() >> Bax5mSmac(), kA13)
Rule('rule14', Bax5mSmac() >> Bax5() + mSmac(), kA14)
Rule('rule15', Bax5mSmac() >> Bax5() + cSmac(), kA15)

Rule('rule16', Bax6() + mSmac() >> Bax6mSmac(), kA16)
Rule('rule17', Bax6mSmac() >> Bax6() + mSmac(), kA17)
Rule('rule18', Bax6mSmac() >> Bax6() + cSmac(), kA18)

Parameter('Bax1_0', 40.0)
Parameter('Bax2_0', 0.0)
Parameter('Bax3_0', 0.0)
Parameter('Bax4_0', 0.0)
Parameter('Bax5_0', 0.0)
Parameter('Bax6_0', 0.0)
Parameter('Bax4mSmac_0', 0.0)
Parameter('Bax5mSmac_0', 0.0)
Parameter('Bax6mSmac_0', 0.0)
Parameter('mSmac_0', 50.0)
Parameter('cSmac_0', 0.0)

Initial(Bax1(), Bax1_0)
Initial(Bax2(), Bax2_0)
Initial(Bax3(), Bax3_0)
Initial(Bax4(), Bax4_0)
Initial(Bax5(), Bax5_0)
Initial(Bax6(), Bax6_0)
Initial(Bax4mSmac(), Bax4mSmac_0)
Initial(Bax5mSmac(), Bax5mSmac_0)
Initial(Bax6mSmac(), Bax6mSmac_0)
Initial(mSmac(), mSmac_0)
Initial(cSmac(), cSmac_0)

if __name__ == '__main__':
    tspan = np.arange(20000)
    x = odesolve(model, tspan)

    # Plot trajectory of each pore
    for i in range(1, max_size + 1):
        observable = 'oBax%d' % i
        # Map pore size to the central 50% of the YlOrBr color map
        color = plt.cm.YlOrBr(float(i) / max_size / 2 + 0.25)
        plt.plot(tspan, x[observable], c=color, label=observable[1:])

    # Plot BaxmSmac compounds
    plt.plot(tspan, x['oBax4mSmac'], 'k-', label='Bax4mSmac')
    plt.plot(tspan, x['oBax5mSmac'], 'k--', label='Bax5mSmac')
    plt.plot(tspan, x['oBax6mSmac'], 'k.', label='Bax6mSmac')

    # Plot Smac species
    plt.plot(tspan, x['omSmac'], c='magenta', label='mSmac')
    plt.plot(tspan, x['ocSmac'], c='cyan', label='cSmac')

    # plt.xscale('log')
    plt.legend(loc='upper left')
    plt.show()
