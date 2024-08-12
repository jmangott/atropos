from scipy.integrate import solve_ivp
import time

from ode_helper import *
plt.style.use("./scripts/output/notebooks/custom_style.mplstyle")

# Simulation parameters
t = 500
n = np.array([51, 51], dtype="int64")
slice_vec = np.array([0, 0], dtype="int64")
r = 5
m = n.size
m1 = 1
dx = np.prod(n)

# Reaction system
# Reaction parameters
b = 0.4
c = 0.05

# Stoichiometric vectors
nu0 = np.array([-1, 0])
nu1 = np.array([0, -1])
nu2 = np.array([1, 0])
nu3 = np.array([0, 1])

# Propensity functions
@njit
def prop0(x: np.ndarray) -> float:
    return c * x[0]
@njit
def prop1(x: np.ndarray) -> float:
    return c * x[1]
@njit
def prop2(x: np.ndarray) -> float:
    return b / (b + x[1])
@njit
def prop3(x: np.ndarray) -> float:
    return b / (b + x[0])

# RHS of the CME
@njit
def cme(t: float, P: np.ndarray, interval: np.ndarray) -> np.ndarray:
    m = interval.size
    null = np.zeros(m)
    result = (
        evaluateProp(prop0, nu0, interval) * shiftArray(P, nu0, interval) +
        evaluateProp(prop1, nu1, interval) * shiftArray(P, nu1, interval) +
        evaluateProp(prop2, nu2, interval) * shiftArray(P, nu2, interval) +
        evaluateProp(prop3, nu3, interval) * shiftArray(P, nu3, interval)
    )
    result -= (
        evaluateProp(prop0, null, interval) +
        evaluateProp(prop1, null, interval) +
        evaluateProp(prop2, null, interval) +
        evaluateProp(prop3, null, interval)
    ) * P
    return result

# Set up the initial condition
C = 0.5 * np.array([[75, -15], [-15, 75]])
Cinv = np.linalg.inv(C)
mu = np.array([30, 5])

def eval_P0(x: np.ndarray) -> float:
    return np.exp(-0.5 * np.dot(np.transpose(x - mu), np.dot(Cinv, (x - mu))))

P0 = constructP0(eval_P0, n)

# Solve the system
t_step = 10
t_eval = np.arange(0, t + t_step, t_step)

t_start = time.time_ns()
sol = solve_ivp(lambda t, P: cme(t, P, n), [0, t + 1], P0, method='RK45', t_eval=t_eval)
t_stop = time.time_ns()
wall_time = (t_stop - t_start) * 1e-9

y = np.zeros((t_eval.size, dx))
for i in range(t_eval.size):
    y[i, :] = sol.y[:, i]

# Calculate output
P_full, P_marginal, _, P_sliced, _, P_best_approximation = calculateObservables(
    y, n, r, m1, slice_vec, np.array([0, 1], dtype="int64"))

with open("scripts/reference_solutions/ts_ode_ref.npy", "wb") as f:
    np.save(f, P_full)
    np.save(f, P_marginal)
    np.save(f, P_sliced)
    np.save(f, P_best_approximation)
    np.save(f, n)
    np.save(f, wall_time)

P_mat = P_full[-1].reshape(n[0], n[1], order="F");
u, s, vh = np.linalg.svd(P_mat)
x0 = u[:, :5]
print(s[:5])

# Full probability distribution
xx1, xx2 = np.meshgrid(np.arange(n[0]), np.arange(n[1]))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 8))
levels = np.linspace(np.amin(P_mat), np.amax(P_mat), 9)
ax1.contour(xx1, xx2, P_mat, levels=levels)
ax1.set_title("exact")

ax2.contour(xx1, xx2, (u[:, :2] * s[:2]) @ vh[:2, :], levels=levels)
ax2.set_title("$r = 2$")
plt.setp((ax1, ax2), xlim=[0, 30], ylim=[0, 30], xlabel="$x_0$", ylabel="$x_1$", aspect="equal")

plt.subplots_adjust(wspace=0.3)
plt.savefig("plots/toggle_switch_comparison.pdf")

# First two basis functions
fig, axs = plt.subplots(2, 2, figsize=(6, 6))
axs[0, 0].plot(np.arange(n[0]), -u[:, 0])
axs[0, 0].set_xlabel("$x_0$")
axs[0, 0].set_ylabel("$P_{{0,H}}(x_0)$")
axs[0, 1].plot(np.arange(n[1]), -vh[0, :])
axs[0, 1].set_xlabel("$x_1$")
axs[0, 1].set_ylabel("$P_{{1,L}}(x_1)$")

axs[1, 0].plot(np.arange(n[0]), -u[:, 1])
axs[1, 0].set_xlabel("$x_0$")
axs[1, 0].set_ylabel("$P_{{0,L}}(x_0)$")
axs[1, 1].plot(np.arange(n[1]), -vh[1, :])
axs[1, 1].set_xlabel("$x_1$")
axs[1, 1].set_ylabel("$P_{{1,H}}(x_1)$")

fig.tight_layout()
plt.savefig("plots/toggle_switch_lr_factors.pdf")
