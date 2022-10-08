import numpy as np

r = 5
kk1 = np.array([1])
kk2 = np.array([1])
x_max1 = np.array([50])
x_max2 = np.array([50])
dim1 = kk1 * x_max1 + 1
dim2 = kk2 * x_max2 + 1

# Set up p0
# This function returns the initial condition of the problem; p0(x) = eval_p0(x)
def eval_p0(x):
    C = 0.5 * np.array([[75, -15], [-15, 75]])
    Cinv = np.linalg.inv(C)
    mu = np.array([30, 5])
    return np.exp(-0.5 * np.transpose(x - mu) @ Cinv @ (x - mu))

# Fill p0 according to the initial condition
p0 = np.zeros((np.prod(dim1), np.prod(dim2)))
n1_vec = np.zeros(dim1.size)
n2_vec = np.zeros(dim2.size)
for i in range(np.prod(dim1)):
    stride1 = 1
    for k, el_dim1 in enumerate(dim1):
        if k == dim1.size - 1:
            n1_vec[k] = i // stride1
        else:
            n1_vec[k] = (i % el_dim1) // stride1
            stride1 = stride1 * el_dim1
    x1 = n1_vec * x_max1 / (dim1 - 1.0)
    for j in range(np.prod(dim2)):
        stride2 = 1
        for l, el_dim2 in enumerate(dim2):
            if l == dim2.size - 1:
                n2_vec[l] = j // stride2
            else:
                n2_vec[l] = (j % el_dim2) // stride2
                stride2 = stride2 * el_dim2
        x2 = n2_vec * x_max2 / (dim2 - 1.0)
        x = np.concatenate([x1, x2])
        p0[i, j] = eval_p0(x)

# Normalize p0
p0 = p0 / np.sum(p0)

# SVD of p0
u, s, vh = np.linalg.svd(p0, full_matrices = False)

# Use only the first `r` singular values
fmt = '%1.8f'
np.savetxt("input/u.csv", u[:, :r], delimiter=",", fmt = fmt)
np.savetxt("input/s.csv", np.diag(s[:r]), delimiter=",", fmt = fmt)
np.savetxt("input/vh.csv", vh[:, :r], delimiter = ",", fmt = fmt)
