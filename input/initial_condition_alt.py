import numpy as np
from index_functions import *
from scipy.special import factorial

# Lambda phage
# kk1 = np.array([1, 1])
# kk2 = np.array([1, 1, 1])
# nx1 = np.array([16, 41])
# nx2 = np.array([11, 11, 11])

# x1_0 = np.array([1, 0], dtype="int64")
# x2_0 = np.array([0, 0, 0], dtype="int64")

# # TGFb8
kk1 = np.array([1, 1, 1, 1, 1, 1])
kk2 = np.array([1, 1])
nx1 = np.array([5, 5, 21, 21, 26, 21])
nx2 = np.array([151, 151])

x1_0 = np.array([2, 0, 6, 12, 16, 3], dtype="int64")
x2_0 = np.array([20, 133], dtype="int64")

dx1 = np.prod(nx1)
dx2 = np.prod(nx2)

m1 = nx1.size
m2 = nx2.size

x1 = np.zeros(dx1)
x2 = np.zeros(dx2)

for i in range(dx1):
    vec_index = CombIndexToVecIndex(i, nx1)
    if (vec_index == x1_0).all():
        x1[i] = 1.0

for i in range(dx2):
    vec_index = CombIndexToVecIndex(i, nx2)
    if (vec_index == x2_0).all():
        x2[i] = 1.0

fmt = '%1.16f'
np.savetxt("input/x1_input.csv", x1, delimiter = ",", fmt = fmt)
np.savetxt("input/x2_input.csv", x2, delimiter = ",", fmt = fmt)
