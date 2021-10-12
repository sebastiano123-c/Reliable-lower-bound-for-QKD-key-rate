
"""
EB simple BB84

@author: Sebastiano Cocchi
"""
from src import qkd
import numpy as np
import matplotlib.pyplot as plt
import time

start_time = time.time()

# new simulation
sim = qkd.QKD(dim_a=2, dim_b=2, n_of_signal_states=4,
    list_states_a=[qkd.zero, qkd.one, qkd.plus, qkd.minus], list_of_prob_a=[0.25, 0.25, 0.25, 0.25],
    list_states_b=[qkd.zero, qkd.one, qkd.plus, qkd.minus], list_of_prob_b=[0.25, 0.25, 0.25, 0.25])

# define qber interval
qber, key_th, key_primal, key_dual = np.linspace(0., 0.12, 15), [], [], []

# iteration
for ii in qber:
    print("\n QBER =", ii)

    # theorical
    hp = qkd.binary_entropy(ii)

    # apply quantum channel
    sim.apply_quantum_channel(qkd.depolarizing_channel(2*ii))

    # set constraints
    gamma = []
    for jj in sim.povm:
        gamma.append(np.trace( jj @ sim.rho_ab))
    sim.set_constraints(gamma, sim.povm)

    # compute primal and dual problem
    sim.compute_primal(epsilon=1e-11, maxit=1000, finesse=100, solver_name="MOSEK")
    sim.compute_dual(solver_name="MOSEK")

    # append
    key_th.append(1 - 2*hp)
    key_primal.append( sim.primal_sol - hp)
    key_dual.append( sim.dual_sol - hp)

    # print result
    print("--- --- --- --- --- --- --- --- ---")
    print(" step 1 =", sim.primal_sol)
    print(" step 2 =", sim.dual_sol)

# CPU time
print( "\n CPU time: ", time.time() - start_time, "s")

# plot
fig, ax = plt.subplots(figsize=(20, 11))
ax.plot(qber, key_th, "--", linewidth=1.2, alpha=0.5, label="theorical")
ax.plot(qber, key_primal, "o", alpha=0.5, label="step 1")
ax.plot(qber, key_dual, ".", alpha=0.5, label="step 2")
plt.xlabel("QBER")
plt.ylabel("Secret key rate")
plt.title("Reliable lower bound EB BB84 with public announcement and sifting")
plt.ylim([0., None])
plt.xlim([0., None])
plt.legend(loc='best')
plt.grid()
# plt.savefig("analysis/bb84"+str(100*Pz)+".png")
plt.show()