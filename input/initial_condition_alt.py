import numpy as np
from index_functions import *

# Parameters
binsize1 = np.array([1, 1, 1])
binsize2 = np.array([1, 1])
nx1 = np.array([21, 61, 11])
nx2 = np.array([71, 61])
x1_0 = np.array([15, 0, 0], dtype="int64")
x2_0 = np.array([60, 0], dtype="int64")

# Code
dx1 = np.prod(nx1)
dx2 = np.prod(nx2)

m1 = nx1.size
m2 = nx2.size

x1 = np.zeros(dx1)
x2 = np.zeros(dx2)

stride = 1
idx1 = 0
for i, el in enumerate(x1_0):
    idx1 += el * stride
    stride *= nx1[i]

stride = 1
idx2 = 0
for i, el in enumerate(x2_0):
    idx2 += el * stride
    stride *= nx2[i]

x1[idx1] = 1.0 / np.prod(binsize1)
x2[idx2] = 1.0 / np.prod(binsize2)

fmt = '%1.16f'
np.savetxt("input/x1_input.csv", x1, delimiter = ",", fmt = fmt)
np.savetxt("input/x2_input.csv", x2, delimiter = ",", fmt = fmt)
