import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# Parameters
t = 500
n1 = 51
n2 = 51
r = 5

# Shift function
def shiftArray(arr, shift, dim):
    result = np.roll(arr, -shift)
    if abs(shift) == dim[1]:
        if shift > 0:
            result[(dim[1] - 1) * dim[0]:dim[1] * dim[0]] = 0
        elif shift < 0:
            result[0:dim[1]] = 0
    elif abs(shift) == 1:
        if shift > 0:
            result[(dim[1] - 1):dim[1] * dim[0]:dim[1]] = 0
        elif shift < 0:
            result[0:dim[1] * dim[1]:dim[0]] = 0
    return result

# RHS of the CME
def cme(t, P, dim):
    b = 0.4
    c = 0.05
    result = np.empty_like(P)
    result = c * np.kron(np.arange(dim[0]) + 1, np.ones(dim[1])) * shiftArray(P, dim[1], dim) \
        + c * np.kron(np.ones(dim[0]), np.arange(dim[1]) + 1) * shiftArray(P, 1, dim) \
            + b * np.kron(np.ones(dim[0]), np.reciprocal(b * np.ones(dim[1]) + np.arange(dim[1]))) * (shiftArray(P, -dim[1], dim) - P) \
                + b * np.kron(np.reciprocal(b * np.ones(dim[0]) + np.arange(dim[0])), np.ones(dim[1])) * (shiftArray(P, -1, dim) - P)
    result -= c * (np.kron(np.arange(dim[0]), np.ones(dim[1])) + np.kron(np.ones(dim[0]), np.arange(dim[1]))) * P
    return result

# Initial value
P0 = np.zeros(n1 * n2)
C = 0.5 * np.array([[75, -15], [-15, 75]])
Cinv = np.linalg.inv(C)
mu = np.array([30, 5])
for i in range(n1):
    for j in range(n2):
        x = np.array([i, j])
        P0[j + n2 * i] = np.exp(-0.5 *
                                np.dot(np.transpose(x - mu), np.dot(Cinv, (x - mu))))
P0 = P0 / np.sum(P0)

# Solve the system
t_eval = np.arange(0, t + 10, 10)
sol = solve_ivp(lambda t, P: cme(t, P, (n1, n2)), [0, t + 1], P0, method='RK45', t_eval=t_eval)

# Save reference solution
np.save('reference_solution/reference_solution_ts', sol.y)

# Calculate best approximation
best_approximation = np.zeros(sol.y.shape)
for i in range(t_eval.size):
    # SVD of p0
    P = sol.y[:, i].reshape((n1, n2))
    u, s, vh = np.linalg.svd(P, full_matrices=False)

    # Use only the first `r` singular values
    X1 = u[:, :r]
    S = s[:r]
    X2h = vh[:r, :]
    Pba = (X1 * S) @ X2h
    best_approximation[:, i] = Pba.flatten()

# Save best approximation
np.save('reference_solution/best_approximation_ts', best_approximation)
