from src import qkd
import numpy as np
import matplotlib.pyplot as plt
import time

start_time = time.time()

# dimensions
da = 2
db = 2
nst = 4
dtot = da*db

# parameters
epsilon = 1e-10
Pz = 0.5
start, stop, step = 0., 0.12, 15
maxit = 1000
finesse = 20
solver_name = "MOSEK"

# define states
states = [qkd.zero, qkd.one, qkd.plus, qkd.minus]

# ALICE probabilities
Px = (1. - Pz)
ProbAlice = [Pz/2., Pz/2., Px/2, Px/2]
if (np.sum(ProbAlice) != 1): print("ProbAlice != 1")

# BOB porbabilities
BS = [0.7, 0.3] # beamsplitter
ProbBob = [BS[0]/2., BS[0]/2., BS[1]/2., BS[1]/2.]
if (np.sum(ProbBob) != 1): print("ProbBob != 1")

# new simulation
sim = qkd.QKD(da, db, nst, states, ProbAlice, states, ProbBob)

# define qber interval
qber = np.linspace(start, stop, step)

# result arrays
key_th      = []
key_primal  = []
key_dual    = []

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
    sim.compute_primal(epsilon, maxit, finesse, "MOSEK")
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
plt.title("Reliable lower bound P&M BB84 with public announcement and sifting")
plt.ylim([0., None])
plt.xlim([0., None])
plt.legend(loc='best')
plt.grid()
# plt.savefig("analysis/bb84"+str(100*Pz)+".png")
plt.show()