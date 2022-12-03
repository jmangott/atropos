import numpy as np

r = 5
kk1 = np.array([1])
kk2 = np.array([1])
nx1 = np.array([51])
nx2 = np.array([51])
x_max1 = int((nx1 - 1) / kk1)
x_max2 = int((nx2 - 1) / kk2)
dim1 = np.prod(nx1)
dim2 = np.prod(nx2)

# Set up p0
# This function returns the initial condition of the problem; p0(x) = eval_p0(x)
def eval_p0(x):
    C = 0.5 * np.array([[75, -15], [-15, 75]])
    Cinv = np.linalg.inv(C)
    mu = np.array([30, 5])
    return np.exp(-0.5 * np.dot(np.transpose(x - mu), np.dot(Cinv, (x - mu))))

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
        stride2 = 1
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
fmt = '%1.8f'
np.savetxt("x1_input.csv", u[:, :r], delimiter = ",", fmt = fmt)
np.savetxt("s_input.csv", np.diag(s[:r]), delimiter = ",", fmt = fmt)
np.savetxt("x2_input.csv", vh[:, :r], delimiter = ",", fmt = fmt)
