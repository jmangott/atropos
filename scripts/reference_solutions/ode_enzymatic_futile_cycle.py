from scipy.integrate import solve_ivp
from scipy.special import factorial

from scripts.reference_solutions.ode_helper import *
from scripts.index_functions import vecIndexToCombIndex

# Simulation parameters
t = 10.0
n = np.array([128, 4, 4, 128, 4, 4], dtype="int64")
slice_vec = np.zeros(6)
m = n.size
m1 = 3
dx = np.prod(n)

# Reaction system
# Reaction parameters
kp1 = 0.4
kp2 = 100.0
kp3 = 100.0
km1 = 2.0
km2 = 1.0
km3 = 50.0

# Stoichiometric vectors
nu0 = np.array([-1, -1, 1, 0, 0, 0])
nu1 = np.array([1, 1, -1, 0, 0, 0])
nu2 = np.array([0, 0, 0, -1, -1, 1])
nu3 = np.array([0, 0, 0, 1, 1, -1])
nu4 = np.array([0, 1, -1, 1, 0, 0])
nu5 = np.array([1, 0, 0, 0, 1, -1])

# Propensity functions
@njit
def prop0(x: np.ndarray) -> float:
    return kp1 * x[0] * x[1]
@njit
def prop1(x: np.ndarray) -> float:
    return kp2 * x[2]
@njit
def prop2(x: np.ndarray) -> float:
    return km1 * x[3] * x[4]
@njit
def prop3(x: np.ndarray) -> float:
    return km2 * x[5]
@njit
def prop4(x: np.ndarray) -> float:
    return kp3 * x[2]
@njit
def prop5(x: np.ndarray) -> float:
    return km3 * x[5]

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
        evaluateProp(prop5, nu5, interval) * shiftArray(P, nu5, interval)
    )
    result -= (
        evaluateProp(prop0, null, interval) + 
        evaluateProp(prop1, null, interval) +
        evaluateProp(prop2, null, interval) +
        evaluateProp(prop3, null, interval) +
        evaluateProp(prop4, null, interval) +
        evaluateProp(prop5, null, interval)
    ) * P
    return result

# Set up the initial condition
P0 = np.zeros(dx)
idx = vecIndexToCombIndex([30, 2, 0, 90, 2, 0], n)
P0[idx] = 1.0

# Solve the system
t_step = 1
t_eval = np.arange(0, t + t_step, t_step)
sol = solve_ivp(lambda t, P: cme(t, P, n), [0, t + 1], P0, method='RK45', t_eval=t_eval)
y = sol.y.T

# Calculate output
_, P_marginal, _, _, _ = calculateObservables(y, n, slice_vec, np.array([0, 1], dtype="int64"))
with open("scripts/reference_solutions/efc_ode_ref.npz", "wb") as f:
    np.savez(f, P_marginal=P_marginal)