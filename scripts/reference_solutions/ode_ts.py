import numpy as np
from scipy.integrate import solve_ivp

# Parameters
t = 500
n1 = 51
n2 = 51
r = 5

# Shift function
def shiftArray(arr: np.ndarray, shift: int, dim: list) -> np.ndarray:
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
            result[0:(dim[1] - 1) * dim[0] + 1:dim[1]] = 0
    return result

# RHS of the CME
def cme(t: float, P: np.ndarray, dim: list) -> np.ndarray:
    b = 0.4
    c = 0.05
    result = c * np.kron(np.arange(dim[0]) + 1, np.ones(dim[1])) * shiftArray(P, dim[1], dim) \
        + c * np.kron(np.ones(dim[0]), np.arange(dim[1]) + 1) * shiftArray(P, 1, dim) \
            + b * np.kron(np.ones(dim[0]), np.ones(dim[1]) / (b * np.ones(dim[1]) + np.arange(dim[1], dtype="f8"))) * (shiftArray(P, -dim[1], dim) - P) \
                + b * np.kron(np.ones(dim[0]) / (b * np.ones(dim[0]) + np.arange(dim[0], dtype="f8")), np.ones(dim[1])) * (shiftArray(P, -1, dim) - P)
    result -= c * (np.kron(np.arange(dim[0]), np.ones(dim[1])) + np.kron(np.ones(dim[0]), np.arange(dim[1]))) * P
    return result

# Initial value
P0 = np.zeros(n1 * n2, dtype="f8")
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
t_step = 10
t_eval = np.arange(0, t + t_step, t_step)
sol = solve_ivp(lambda t, P: cme(t, P, (n1, n2)), [0, t + 1], P0, method='RK45', t_eval=t_eval)

P_full = np.zeros((t_eval.size, n1, n2))
for i in range(t_eval.size):
    P_full[i, :, :] = sol.y[:, i].reshape((n1, n2))

# Calculate marginal distributions
P_marginal = [np.zeros((n_el, t_eval.size), dtype="float64") for n_el in (n1, n2)]
for i in range(sol.y.shape[1]):
    Pref = sol.y[:, i].reshape((n1, n2))
    P_marginal[0][:, i] = np.sum(Pref, axis=1)
    P_marginal[1][:, i] = np.sum(Pref, axis=0)

# Calculate sliced distributions
P_sliced = [np.zeros((n_el, t_eval.size), dtype="float64") for n_el in (n1, n2)]
for i in range(sol.y.shape[1]):
    Pref = sol.y[:, i].reshape((n1, n2))
    P_sliced[0][:, i] = Pref[:, 0]
    P_sliced[1][:, i] = Pref[0, :]

# Calculate best approximation
P_best_approximation = np.zeros((t_eval.size, n1, n2))
for i in range(t_eval.size):
    # SVD of p0
    P = sol.y[:, i].reshape((n1, n2))
    u, s, vh = np.linalg.svd(P, full_matrices=False)

    # Use only the first `r` singular values
    X1 = u[:, :r]
    S = s[:r]
    X2h = vh[:r, :]
    P_best_approximation[i, :, :] = (X1 * S) @ X2h

with open("scripts/reference_solutions/ts_ode_result.npy", "wb") as f:
    np.save(f, P_full)
    np.save(f, P_marginal)
    np.save(f, P_sliced)
    np.save(f, P_best_approximation)
    np.save(f, (n1, n2))