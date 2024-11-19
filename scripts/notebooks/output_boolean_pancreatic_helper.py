import numpy as np

from scripts.tree_class import Tree
from scripts.notebooks.output_helper import *

def saveMoments(filename):
    time_series = TimeSeries(filename)
    time = time_series.time
    moments = time_series.calculateMoments()
    np.savez("{}.npz".format(filename), time=time, moments=moments)

partitions = ["b", "r", "w", "l"]
ranks = [5, 10, 20]
for p in partitions:
    for r in ranks:
        print("partition:", p, "rank:", r)
        saveMoments("output/pancreatic_matrix_p{}_r{}_e_tau1e-2".format(p, r))

# Reference solutions
ranks = [50, 60]
for r in ranks:
    print("reference solution, rank:", r)
    saveMoments("output/pancreatic_matrix_pl_r{}_e_tau1e-2".format(r))