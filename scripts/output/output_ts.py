from output_helper import *
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

# Load datasets
with nc.Dataset("output/ts/output_t50000.nc", "r") as ds:
    grid = grid_info(ds)
    slice_vec = np.array([0, 0])
    P_full = fullDistribution(ds, grid)
    P_marginal = marginalDistributions(ds, grid)
    P_sliced = slicedDistributions(ds, grid, slice_vec)

with open("scripts/reference_solutions/ts_ode_result.npy", "rb") as f:
    P_full_ode = np.load(f, allow_pickle=True)
    P_marginal_ode = np.load(f, allow_pickle=True)
    P_sliced_ode = np.load(f, allow_pickle=True)
    P_best_approximation_ode = np.load(f, allow_pickle=True)
    n_ode = np.load(f)

# Figure 1
title = ["DLR approximation", "RK45"]
fig, axs = plt.subplots(1, 2, figsize = (8, 6))
fig.suptitle("One-dimensional marginal distributions for $t = {}$, $\Delta t = {:.2e}$".format(grid.t, grid.dt))
for i, ax in enumerate(axs.flatten()):
    if i < grid.n.size:
        ax.plot(grid.bin[i] * range(grid.n[i]) +
                grid.liml[i], P_marginal[i], "ro", label="DLRA")
    if i < P_marginal_ode.size:
        ax.plot(range(P_marginal_ode[i].shape[0]), P_marginal_ode[i][:, -1], "bx", label="RK45")
        ax.set_xlabel("$P_{{\mathrm{{MD}}}}(x_" + str(i) + ")$")
        ax.set_xlabel("$x_{{" + str(i + 1) + "}}$")

axs[0].legend()
fig.tight_layout()

# Figure 2
fig, axs = plt.subplots(1, 2, figsize=(8, 6))
fig.suptitle(
    "Sliced distributions for $x_0 = {}$, $t = {}$, $\\tau = {:.2e}$".format(np.array2string(slice_vec, separator=","), grid.t, grid.dt))
for i, ax in enumerate(axs.flatten()):
    if i < grid.n.size:
        ax.plot(grid.bin[i] * range(grid.n[i]) +
                grid.liml[i], P_sliced[i], "ro", label="DLRA")
    if i < P_sliced_ode.size:
        ax.plot(range(P_sliced_ode[i].shape[0]),
                P_sliced_ode[i][:, -1], "bx", label="RK45")
        ax.set_xlabel("$P_{{\mathrm{{MD}}}}(x_" + str(i) + ")$")
        ax.set_xlabel("$x_{{" + str(i + 1) + "}}$")

axs[0].legend()
fig.tight_layout()

# Figure 3
xx1, xx2 = np.meshgrid(np.arange(grid.n1[0]), np.arange(grid.n2[0]), indexing='ij')
fig, ax = plotP2D(P_full, P_full_ode[-1], [xx1, xx2], [xx1, xx2], title)
fig.suptitle(
    "Full probability distribution for $t = {}$, $\\tau = {:.2e}$".format(grid.t, grid.dt))
plt.setp(ax.flat, xlim=[0, 30], ylim=[0, 30], aspect="equal")
fig.tight_layout()

# Figure 4
times = np.arange(0, 500 + 10, 10)
err = np.zeros(times.size)
errba = np.zeros(times.size)

for ts, t in enumerate(times):
    filename = "output/ts/output_t" + str(t * 100) + ".nc"
    with nc.Dataset(filename, "r") as ds:
        P = fullDistribution(ds, grid)
    err[ts] = np.linalg.norm(P - P_full_ode[ts], ord=2)
    errba[ts] = np.linalg.norm(P_best_approximation_ode[ts] - P_full_ode[ts], ord=2)

plt.figure(4, figsize=(8, 5))
plt.plot(times, err)
plt.plot(times, errba)
plt.legend(["DLR approximation", "best-approximation"])
plt.xlabel("$t$")
plt.ylabel("2-norm error")
# plt.ylim([0, 0.005])

plt.show()