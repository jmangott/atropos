import numpy as np
from index_functions import *

# # Toggle switch
# kk1 = np.array([1])
# kk2 = np.array([1])
# nx1 = np.array([301])
# nx2 = np.array([301])

# x1_0 = np.array([0], dtype="int64")
# x2_0 = np.array([10], dtype="int64")

# # Lambda phage
# kk1 = np.array([1, 1])
# kk2 = np.array([1, 1, 1])
# nx1 = np.array([6, 41])
# nx2 = np.array([6, 11, 21])
# # nx1 = np.array([6, 151])
# # nx2 = np.array([6, 11, 21])

# # x1_0 = np.array([1, 0], dtype="int64")
# # x2_0 = np.array([0, 0, 0], dtype="int64")
# x1_0 = np.array([0, 15], dtype="int64")
# x2_0 = np.array([0, 5, 10], dtype="int64")

# # TGFb6
# kk1 = np.array([1, 1, 1, 1])
# kk2 = np.array([1, 1, 1, 1])
# nx1 = np.array([5, 5, 151, 151])
# nx2 = np.array([26, 21, 21, 21])

# x1_0 = np.array([2, 0, 20, 133], dtype="int64")
# x2_0 = np.array([16, 3, 6, 12], dtype="int64")

# Tyson
kk1 = np.array([1, 1])
kk2 = np.array([1, 1, 1])
nx1 = np.array([1521, 1901])
nx2 = np.array([16, 2041, 1251])

x1_0 = np.array([1505, 0], dtype="int64")
x2_0 = np.array([0, 2022, 0], dtype="int64")

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

x1[idx1] = 1.0
x2[idx2] = 1.0

# for i in range(dx1):
#     vec_index = CombIndexToVecIndex(i, nx1)
#     if (vec_index == x1_0).all():
#         x1[i] = 1.0

# for i in range(dx2):
#     vec_index = CombIndexToVecIndex(i, nx2)
#     if (vec_index == x2_0).all():
#         x2[i] = 1.0

fmt = '%1.16f'
np.savetxt("input/x1_input.csv", x1, delimiter = ",", fmt = fmt)
np.savetxt("input/x2_input.csv", x2, delimiter = ",", fmt = fmt)
