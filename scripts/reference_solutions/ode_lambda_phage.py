from scipy.integrate import solve_ivp
from scipy.special import factorial
import time

from ode_helper import *

# Simulation parameters
t = 10
n = np.array([16, 41, 11, 11, 11], dtype="int64")
slice_vec = np.array([0, 9, 1, 1, 1])
r = 4
m = n.size
m1 = 2
dx = np.prod(n)

# Reaction system
# Reaction parameters
a = np.array([0.5, 1.0, 0.15, 0.3, 0.3])
b = np.array([0.12, 0.6, 1.0, 1.0, 1.0])
c = np.array([0.0025, 0.0007, 0.0231, 0.01, 0.01])

# Stoichiometric vectors
nu0 = np.array([1, 0, 0, 0, 0])
nu1 = np.array([0, 1, 0, 0, 0])
nu2 = np.array([0, 0, 1, 0, 0])
nu3 = np.array([0, 0, 0, 1, 0])
nu4 = np.array([0, 0, 0, 0, 1])
nu5 = np.array([-1, 0, 0, 0, 0])
nu6 = np.array([0, -1, 0, 0, 0])
nu7 = np.array([0, 0, -1, 0, 0])
nu8 = np.array([0, 0, 0, -1, 0])
nu9 = np.array([0, 0, 0, 0, -1])

# Propensity functions
@njit
def prop0(x: np.ndarray) -> float:
    return a[0] * b[0] / (b[0] + x[1])
@njit
def prop1(x: np.ndarray) -> float:
    return (a[1] + x[4]) * b[1] / (b[1] + x[0])
@njit
def prop2(x: np.ndarray) -> float:
    return a[2] * b[2] * x[1] / (b[2] * x[1] + 1.0)
@njit
def prop3(x: np.ndarray) -> float:
    return a[3] * b[3] * x[2] / (b[3] * x[2] + 1.0)
@njit
def prop4(x: np.ndarray) -> float:
    return a[4] * b[4] * x[2] / (b[4] * x[2] + 1.0)
@njit
def prop5(x: np.ndarray) -> float:
    return c[0] * x[0]
@njit
def prop6(x: np.ndarray) -> float:
    return c[1] * x[1]
@njit
def prop7(x: np.ndarray) -> float:
    return c[2] * x[2]
@njit
def prop8(x: np.ndarray) -> float:
    return c[3] * x[3]
@njit
def prop9(x: np.ndarray) -> float:
    return c[4] * x[4]

# RHS of the CME
@njit
def cme(t: float, P: np.ndarray, interval: np.ndarray) -> np.ndarray:
    m = interval.size
    null = np.zeros(m)
    result = (
        evaluateProp(prop0, nu0, interval) * shiftArray(P, nu0, interval) + 
        evaluateProp(prop1, nu1, interval) * shiftArray(P, nu1, interval) +
        evaluateProp(prop2, nu2, interval) * shiftArray(P, nu2, interval) +
        evaluateProp(prop3, nu3, interval) * shiftArray(P, nu3, interval) +
        evaluateProp(prop4, nu4, interval) * shiftArray(P, nu4, interval) +
        evaluateProp(prop5, nu5, interval) * shiftArray(P, nu5, interval) +
        evaluateProp(prop6, nu6, interval) * shiftArray(P, nu6, interval) +
        evaluateProp(prop7, nu7, interval) * shiftArray(P, nu7, interval) +
        evaluateProp(prop8, nu8, interval) * shiftArray(P, nu8, interval) +
        evaluateProp(prop9, nu9, interval) * shiftArray(P, nu9, interval)
    )
    result -= (
        evaluateProp(prop0, null, interval) + 
        evaluateProp(prop1, null, interval) +
        evaluateProp(prop2, null, interval) +
        evaluateProp(prop3, null, interval) +
        evaluateProp(prop4, null, interval) +
        evaluateProp(prop5, null, interval) +
        evaluateProp(prop6, null, interval) +
        evaluateProp(prop7, null, interval) +
        evaluateProp(prop8, null, interval) +
        evaluateProp(prop9, null, interval)
    ) * P
    return result

# Set up the initial condition
def eval_P0(x: np.ndarray) -> float:
    abs_x = np.sum(x)
    if (abs_x <= 3):
        P0 = factorial(3) * (0.05 ** abs_x) * \
            ((1.0 - 5 * 0.05) ** (3 - abs_x)) / \
            (np.prod(factorial(x)) * factorial(3 - abs_x))
    else:
        P0 = 0.0
    return P0

P0 = constructP0(eval_P0, n)

# Solve the system
t_step = 1
t_eval = np.arange(0, t + t_step, t_step)

t_start = time.time_ns()
sol = solve_ivp(lambda t, P: cme(t, P, n), [0, t + 1], P0, method='RK45', t_eval=t_eval)
t_stop = time.time_ns()
wall_time = (t_stop - t_start) * 1e-9

y = np.zeros((t_eval.size, dx))
for i in range(t_eval.size):
    y[i, :] = sol.y[:, i]

# Calculate the output
P_full, P_marginal, P_marginal2D, P_sliced, P_sliced2D = calculateObservables(
    y, n, slice_vec, np.array([0, 1], dtype="int64"))
P_best_approximation = calculateBestApproximation(y, n, 5, m1)

with open("scripts/reference_solutions/lp_ode_ref_r5.npz", "wb") as f:
    np.savez(f, P_full=P_full, P_best_approximation=P_best_approximation, wall_time=wall_time)