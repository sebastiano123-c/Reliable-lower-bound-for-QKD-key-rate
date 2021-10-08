from src import qkd
import numpy as np
import matplotlib.pyplot as plt
import time

start_time = time.time()

# parameters
epsilon = 1e-11
Pz = 0.1
start, stop, step = 0., 0.12, 15
maxit = 1000
finesse = 100
solver_name = "MOSEK"

# define qubits
zero = qkd.zero
one  = qkd.one
states = [zero, one, (zero+one)/np.sqrt(2.), (zero-one)/np.sqrt(2.)]

# ALICE probabilities
Px = (1. - Pz)
ProbAlice = [Pz/2., Pz/2., Px/2., 0.]
if (np.sum(ProbAlice) != 1): print("ProbAlice != 1")

# BOB porbabilities
BS = [0.5, 0.5] # beamsplitter
ProbBob = [BS[0]/2., BS[0]/2., BS[1]/2., BS[1]/2.]
if (np.sum(ProbBob) != 1): print("ProbBob != 1")

#  post selection and probability of passing the post selection process
postselect_prob = [Pz*BS[0], Px*BS[1]]
ppass = sum(postselect_prob)

# local measurments (|0><0|, |1><1|, |+><+| and |-><-|)
sigma_00 = np.outer( states[0], np.conj(states[0]) )
sigma_11 = np.outer( states[1], np.conj(states[1]) )
sigma_pp = np.outer( states[2], np.conj(states[2]) )
sigma_mm = np.outer( states[3], np.conj(states[3]) )

# new simulation
sim = qkd.QKD(2,2,4,states, ProbAlice, states, ProbBob)

# possible outcome measurments
ZA = np.zeros((2, 4, 4))*1j
ZA[0] = np.kron(sigma_00, np.eye(2)) #( |0><0| + |+><+| ) .o. id_2 --> bit 0
ZA[1] = np.kron(sigma_11, np.eye(2)) #( |1><1| + |-><-| ) .o. id_2 --> bit 1
sim.set_pinching_channel(ZA)

# Gamma for constraints
Ez = np.kron(sigma_00, sigma_11) + np.kron(sigma_11, sigma_00)
Ex = np.kron(sigma_pp, sigma_mm) + np.kron(sigma_mm, sigma_pp)
Gamma = [Ez, Ex]
Gamma = [G*p for G, p in zip(Gamma, postselect_prob)]

qber = np.linspace(start, stop, step)

key_th     = []
key_primal = []
key_dual   = []

for ii in qber:
    print("\n QBER =", ii)
    # theorical
    hp = qkd.binary_entropy(ii)

    # apply quantum channel (in this case depolarization)
    sim.apply_quantum_channel(qkd.depolarizing_channel(2*ii))

    # set contraints
    gamma = np.array([ii, ii])*postselect_prob
    sim.set_constraints(gamma, Gamma)

    # compute primal and dual problem
    sim.compute_primal(epsilon, maxit, finesse)
    sim.compute_dual()

    key_th.append( 1 - 2*hp )
    key_primal.append( sim.primal_sol - hp )
    key_dual.append( sim.dual_sol - hp )
    
    print("--- --- --- --- --- --- --- --- ---")
    print(" step 1 =", sim.primal_sol)
    print(" step 2 =", sim.dual_sol)

print( "\n CPU time: ", time.time() - start_time , "s")

fig, ax = plt.subplots(figsize=(20, 11))
ax.plot(qber, key_th, "--", linewidth=1.2, alpha=0.5, label="theorical")
ax.plot(qber, key_primal, "o", alpha=0.5, label="step 1")
ax.plot(qber, key_dual, ".", alpha=0.5, label="step 2")
plt.xlabel("QBER")
plt.ylabel("Secret key rate")
plt.title("Reliable lower bound P&M BB84")
plt.ylim([0., None])
plt.xlim([0., None])
plt.legend(loc='best')
plt.grid()
plt.savefig("analysis/simple_bb84_"+str(100*Pz)+".png")
plt.show()