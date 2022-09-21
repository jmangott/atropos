import numpy as np

r = 5
x_max = 10
n_x = 50

# set up p0
def eval_p0(x):
    C = 0.5 * np.array([[75, -15], [-15, 75]])
    Cinv = np.linalg.inv(C)
    mu = np.array([30, 5])
    return np.exp(-0.5 * np.transpose(x - mu) @ Cinv @ (x - mu))

p0 = np.zeros((n_x, n_x))
for i in range(n_x):
    x1 = i * x_max / (n_x - 1.0)
    for j in range(n_x):
        x2 = j * x_max / (n_x - 1.0)
        x = np.array([x1, x2])
        p0[i, j] = eval_p0(x)

# normalize p0
p0 = p0 / np.sum(p0)

# SVD of p0
u, s, vh = np.linalg.svd(p0, full_matrices = False, hermitian = True)

# use only the first r singular values
fmt = '%1.8f'
np.savetxt("input/u.csv", u[:, :r], delimiter=",", fmt = fmt)
np.savetxt("input/s.csv", s[:r], delimiter=",", fmt = fmt)
np.savetxt("input/vh.csv", vh[:, :r], delimiter = ",", fmt = fmt)
