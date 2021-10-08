import numpy as np
from scipy.linalg import sqrtm
from src import qkd
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
BS = [0.5, 0.5] # beamsplitter
ProbBob = [BS[0]/2., BS[0]/2., BS[1]/2., BS[1]/2.]
if (np.sum(ProbBob) != 1): print("ProbBob != 1")

# local measurments
k0b0 = np.outer( states[0], np.conj(states[0]) ) # |0><0|
k1b1 = np.outer( states[1], np.conj(states[1]) ) # |1><1|
kpbp = np.outer( states[2], np.conj(states[2]) ) # |+><+|
kmbm = np.outer( states[3], np.conj(states[3]) ) # |-><-|

# define identities for convinience
id_2 = np.eye(2)
id_4 = np.eye(4)
id_a = np.eye(da)
id_b = np.eye(db)
id_tot = np.eye(dtot)

# After the qubit sending, Alice can measure A using the POVM
POVMA = [
    0.5*k0b0,
    0.5*k1b1,
    0.5*kpbp,
    0.5*kmbm
] # which have dimension da-by-da and satisfy POVM properties
qkd.is_povm(POVMA, "POVMA")

# On the other hand, Bob can measure using the POVM
POVMB = [
    0.5*k0b0,
    0.5*k1b1,
    0.5*kpbp,
    0.5*kmbm
] # which have dimension db-by-db and satisfy POVM properties
qkd.is_povm(POVMB, "POVMB")

# The POVM of the entire system is given by
POVM = []
for ii in POVMA:
        for jj in POVMB:
            POVM.append( np.kron( ii, jj ) )
qkd.is_povm(POVM, "POVM")

# 1) PUBLIC ANNOUNCEMENT: kraus operators of A dim = da*4-by-4
KA = [ 
    np.kron( sqrtm(POVMA[0]), np.kron(qkd.zero[:, np.newaxis], qkd.zero[:, np.newaxis])) +
    np.kron( sqrtm(POVMA[1]), np.kron(qkd.zero[:, np.newaxis], qkd.one[:, np.newaxis])),
    np.kron( sqrtm(POVMA[2]), np.kron(qkd.one[:, np.newaxis] , qkd.zero[:, np.newaxis])) +
    np.kron( sqrtm(POVMA[3]), np.kron(qkd.one[:, np.newaxis] , qkd.one[:, np.newaxis]))
] #   which satisfy Kraus property
qkd.is_kraus(KA, "KA")

#   kraus operators of B dim = db*4-by-4
KB = [ 
    np.kron(sqrtm(POVMB[0]), np.kron(qkd.zero[:, np.newaxis], qkd.zero[:, np.newaxis])) +
    np.kron(sqrtm(POVMB[1]), np.kron(qkd.zero[:, np.newaxis], qkd.one[:, np.newaxis])),
    np.kron(sqrtm(POVMB[2]), np.kron(qkd.one[:, np.newaxis] , qkd.zero[:, np.newaxis])) +
    np.kron(sqrtm(POVMB[3]), np.kron(qkd.one[:, np.newaxis] , qkd.one[:, np.newaxis]))
] #   which satisfy Kraus property
qkd.is_kraus(KB, "KB")

#   The total Kraus representation of the Public Announcement is
K = []
for ii in KA:
        for jj in KB:
            K.append( np.kron(ii, jj))
qkd.is_kraus(K, "K")

# 2) SIFTING PHASE: acts like a projector with dimension [da*4*db*4]-by-[da*4*db*4]
proj = np.kron( id_a, np.kron( k0b0, np.kron( np.kron(id_2, id_b), np.kron( k0b0, id_2 ))) ) +\
       np.kron( id_a, np.kron( k1b1, np.kron( np.kron(id_2, id_b), np.kron( k1b1, id_2 ))) )

# 3) KEY MAP: is a isometry which creates a new register R which stores the information on the bits and it is a [2*da*4*db*4]-by-[2*da*4*db*4] matrix
V = np.kron( qkd.zero[:, np.newaxis], np.kron( np.kron(id_a, id_2), np.kron( k0b0, np.kron(id_b, id_4)) )) +\
    np.kron( qkd.one[:, np.newaxis] , np.kron( np.kron(id_a, id_2), np.kron( k1b1, np.kron(id_b, id_4)) ))

# 4) PINCHING CHANNEL: decohere the register R. It has the effect of making R a classical register
pinching = [ np.kron( k0b0 , np.kron( np.kron(id_a, np.kron( id_4, id_b )), id_4 )),
             np.kron( k1b1 , np.kron( np.kron(id_a, np.kron( id_4, id_b )), id_4 )) ]

# new simulation
sim = qkd.QKD(da, db, nst, states, ProbAlice, states, ProbBob)

# set protocol parts
sim.set_povm(POVM)
sim.set_public_string_announcement(K)
sim.set_sifting_phase_postselection(proj)
sim.set_key_map(V)
sim.set_pinching_channel(pinching)

# construct the operative basis with povms
sim.get_full_hermitian_operator_basis()

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
    for jj in sim.orth_set_a:
        gamma.append(np.trace( jj @ sim.rho_ab))
    for jj in sim.povm:
        gamma.append(np.trace( jj @ sim.rho_ab))
    sim.set_constraints(gamma, np.concatenate([sim.orth_set_a, sim.povm]))

    # compute primal and dual problem
    sim.compute_primal(epsilon, maxit, finesse, "MOSEK")
    sim.compute_dual("MOSEK")

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
plt.savefig("analysis/bb84"+str(100*Pz)+".png")
plt.show()
