import numpy as np
from index_functions import *
from scipy.special import factorial

# # Toggle switch
# r = 5
# kk1 = np.array([1])
# kk2 = np.array([1])
# nx1 = np.array([51])
# nx2 = np.array([51])

# Lambda phage
r = 4
kk1 = np.array([1, 1])
kk2 = np.array([1, 1, 1])
nx1 = np.array([16, 41])
nx2 = np.array([11, 11, 11])

x_max1 = np.rint((nx1 - np.ones(nx1.size)) / kk1)
x_max2 = np.rint((nx2 - np.ones(nx2.size)) / kk2)
dim1 = np.prod(nx1)
dim2 = np.prod(nx2)

# Set up p0
# This function returns the initial condition of the problem; p0(x) = eval_p0(x)
def eval_p0(x):
    # # Toggle switch
    # C = 0.5 * np.array([[75, -15], [-15, 75]])
    # Cinv = np.linalg.inv(C)
    # mu = np.array([30, 5])
    # p0 = np.exp(-0.5 * np.dot(np.transpose(x - mu), np.dot(Cinv, (x - mu))))

    # # Lambda phage
    # abs_x = np.sum(x) # np.linalg.norm(x)
    # if (abs_x <= 3):
    #     p0 = factorial(3) * (0.05 ** abs_x) * ((1.0 - 5 * 0.05) ** (3 - abs_x)) / (np.prod(factorial(x)) * factorial(3 - abs_x))
    # else:
    #     p0 = 0.0

    # Lambda phage x1 = 1
    if ((x == [1, 0, 0, 0, 0]).all()):
        p0 = 1.0
    else:
        p0 = 0.0

    return p0

# Fill p0 according to the initial condition
p0 = np.zeros((dim1, dim2))
n1_vec = np.zeros(nx1.size, dtype = int)
n2_vec = np.zeros(nx2.size, dtype = int)
for i in range(dim1):
    n1_vec = CombIndexToVecIndex(i, nx1)
    x1 = n1_vec * x_max1 / (nx1 - 1.0)
    for j in range(dim2):
        num_j = j
        n2_vec = CombIndexToVecIndex(j, nx2)
        x2 = n2_vec * x_max2 / (nx2 - 1.0)
        x = np.concatenate([x1, x2])
        p0[i, j] = eval_p0(x)

# Normalize p0
p0 = p0 / np.sum(p0)

# SVD of p0
u, s, vh = np.linalg.svd(p0, full_matrices = False)

# Use only the first `r` singular values
X1 = u[:, :r]
S = np.diag(s[:r])
X2 = vh[:r, :].T

fmt = '%1.16f'
np.savetxt("input/x1_input.csv", X1, delimiter = ",", fmt = fmt)
np.savetxt("input/s_input.csv", S, delimiter = ",", fmt = fmt)
np.savetxt("input/x2_input.csv", X2, delimiter = ",", fmt = fmt)
