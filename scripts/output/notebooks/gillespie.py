from numba import njit
import numpy as np
import random

def performTimestep(state, reaction_system):
    r1 = random.random()
    r2 = random.random()
    alpha = np.ones(reaction_system.size())

    for i, reaction in enumerate(reaction_system.reactions):
        reactants = reaction.propensity.keys()
        for reactant in reactants:
            alpha[i] *= reaction.propensity[reactant](state[reactant])

    alpha_0 = np.sum(alpha)
    tau = -np.log(r1) / alpha_0

    alpha_i = 0
    for i in range(reaction_system.size()):
        alpha_i_next = alpha_i + alpha[i]
        if r2 >= alpha_i / alpha_0 and r2 < alpha_i_next / alpha_0:
            break
        alpha_i = alpha_i_next

    dn = np.array(reaction_system.reactions[i].nu, dtype="int64")

    return tau, dn


def calculateTrajectory(sample_times, t_stop, initial_state, reaction_system):
    t = 0
    t_idx = 0
    state = np.copy(initial_state)
    trajectory = np.zeros((sample_times.size, initial_state.size), dtype="int")

    while t < t_stop:
        for st in sample_times[t_idx:]:
            if t >= st:
                trajectory[t_idx, :] = state
                t_idx += 1
            else:
                break

        tau, dn = performTimestep(state, reaction_system)
        t += tau
        state += dn

    return trajectory

def calculateFullTrajectory(t_stop, initial_state, reaction_system):
    t0 = 0
    t = t0
    time = [t0]
    state = np.copy(initial_state)
    trajectory = np.copy(initial_state)

    while t < t_stop:
        tau, dn = performTimestep(state, reaction_system)
        t += tau
        state += dn
        trajectory = np.vstack((trajectory, state))
        time.append(t)

    return time, trajectory


def runSimulation(sample_times, t_stop, initial_state, n_runs, reaction_system):
    result = np.zeros((n_runs, sample_times.size, initial_state.size), dtype="int64")

    for i in range(n_runs):
        trajectory = calculateTrajectory(sample_times, t_stop, initial_state, reaction_system)
        result[i, :, :] = trajectory
    return result


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from lambda_phage import reaction_system

    # random.seed(42)
    sample_times = np.linspace(0.0, 10.0, 10+1)
    initial_state = np.zeros(reaction_system.d(), dtype="int64")
    trajectory = runSimulation(sample_times, 10.0, initial_state, 50, reaction_system)

    plt.plot(sample_times, trajectory[:, :, 1].T)
    plt.show()
    # fig, axs = plt.subplots(2, 1)
    # for i in range(5):
    #     time, trajectory = calculateFullTrajectory(100.0, initial_state, reaction_system)
    #     axs[0].plot(time, trajectory[:, 0])
    #     axs[1].plot(time, trajectory[:, 1])
    # plt.show()

