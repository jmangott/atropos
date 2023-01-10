import numpy as np
import matplotlib.pyplot as plt
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

    # Lambda phage
    abs_x = np.linalg.norm(x)
    if (abs_x <= 3):
        p0 = factorial(3) * (0.05 ** abs_x) * ((1.0 - 5 * 0.05) ** (3 - abs_x)) / (np.prod(factorial(x)) * factorial(3 - abs_x))
    else:
        p0 = 0.0

    return p0

# Fill p0 according to the initial condition
p0 = np.zeros((dim1, dim2))
n1_vec = np.zeros(nx1.size, dtype = int)
n2_vec = np.zeros(nx2.size, dtype = int)
for i in range(dim1):
    num_i = i
    for k, el_dim1 in enumerate(nx1):
        if k == nx1.size - 1:
            n1_vec[k] = num_i
        else:
            n1_vec[k] = num_i % el_dim1
            num_i = num_i // el_dim1
    x1 = n1_vec * x_max1 / (nx1 - 1.0)
    for j in range(dim2):
        num_j = j
        for l, el_dim2 in enumerate(nx2):
            if l == nx2.size - 1:
                n2_vec[l] = num_j
            else:
                n2_vec[l] = num_j % el_dim2
                num_j = num_j // el_dim2
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

# P = X1 @ S @ X2.T
# Pmat = np.zeros((nx1[0], nx1[1]))
# for i in range(nx1[0]):
#     for j in range(nx1[1]):
#         Pmat[i, j] = np.sum(P[i + nx1[0] * j, :])
# print(np.sum(P))

# xx1, xx2 = np.meshgrid(range(nx1[0]), range(nx1[1]), indexing = 'ij')
# print(np.shape(xx1))
# print(np.shape(xx2))

# plt.contourf(xx1, xx2, Pmat)
# plt.axis("scaled")
# plt.colorbar()
# plt.show()
