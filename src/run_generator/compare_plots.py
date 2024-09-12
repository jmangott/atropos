import numpy as np
from src.output_helper import readTree

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("TkAgg")

mpl.rcParams["pgf.texsystem"] = "pdflatex"

"""
Generate plots to compare lambda_phage (see if my code is correct)
"""

"""
Compare lambda phage example
"""

lambda_phage_stefan = "/home/stefan/Chemical_master_equation/kinetic-cme/output/lambda_phage_Stefan/output_t800.nc"
lambda_phage_julian = "/home/stefan/Chemical_master_equation/kinetic-cme/output/lambda_phage_Julian/output_t800.nc"

lambda_phage_stefan_tree = readTree(lambda_phage_stefan)
_, sum_lambda_phage_stefan = lambda_phage_stefan_tree.calculateObservables(
    np.zeros(lambda_phage_stefan_tree.grid.d(), dtype=int)
)

lambda_phage_julian_tree = readTree(lambda_phage_julian)
_, sum_lambda_phage_julian = lambda_phage_julian_tree.calculateObservables(
    np.zeros(lambda_phage_julian_tree.grid.d(), dtype=int)
)

plt.plot(
    np.arange(lambda_phage_stefan_tree.grid.n[0]),
    sum_lambda_phage_stefan[0],
    label="lambda_phage_stefan",
)
plt.plot(
    np.arange(lambda_phage_julian_tree.grid.n[0]),
    sum_lambda_phage_julian[0],
    label="lambda_phage_julian",
)
plt.title("Lambda phage comparison")
plt.legend()
plt.show()


"""
Compare bax example
"""

bax_stefan = (
    "/home/stefan/Chemical_master_equation/kinetic-cme/output/bax_Stefan/output_t600.nc"
)
bax_julian = (
    "/home/stefan/Chemical_master_equation/kinetic-cme/output/bax_Julian/output_t600.nc"
)

bax_stefan_tree = readTree(bax_stefan)
_, sum_bax_stefan = bax_stefan_tree.calculateObservables(
    np.zeros(bax_stefan_tree.grid.d(), dtype=int)
)

bax_julian_tree = readTree(bax_julian)
_, sum_bax_julian = bax_julian_tree.calculateObservables(
    np.zeros(bax_julian_tree.grid.d(), dtype=int)
)

plt.plot(np.arange(bax_stefan_tree.grid.n[0]), sum_bax_stefan[0], label="bax_stefan")
plt.plot(np.arange(bax_julian_tree.grid.n[0]), sum_bax_julian[0], label="bax_julian")
plt.title("Bax comparison")
plt.legend()
plt.show()
