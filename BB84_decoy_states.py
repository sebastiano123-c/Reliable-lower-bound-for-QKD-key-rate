
"""
EB BB84 decoy states protocol

@author: Sebastiano Cocchi
"""
from src import qkd
import numpy as np
import matplotlib.pyplot as plt

# basis = np.eye(4)

# # povm a
# povma = []
# for ii in basis:
#     povma.append( np.outer(ii, np.conj(ii)) )

# decoy states
decoy_basis  = np.eye(6)

# states
#  Z basis
vac        = decoy_basis[0] + np.zeros_like(decoy_basis[0])*1j  # |vac>_Z
one_zero_z = decoy_basis[1] + np.zeros_like(decoy_basis[0])*1j  # |1,0>_Z
two_zero_z = decoy_basis[2] + np.zeros_like(decoy_basis[0])*1j  # |2,0>_Z
zero_one_z = decoy_basis[3] + np.zeros_like(decoy_basis[0])*1j  # |0,1>_Z
zero_two_z = decoy_basis[4] + np.zeros_like(decoy_basis[0])*1j  # |0,2>_Z
one_one_z  = decoy_basis[5] + np.zeros_like(decoy_basis[0])*1j  # |1,1>_Z

#  X basis
one_zero_x = (one_zero_z + zero_one_z)/np.sqrt(2)               # |1,0>_X
two_zero_x = (two_zero_z + 2 * one_one_z + zero_two_z)/2        # |2,0>_X
zero_one_x = (one_zero_z - zero_one_z)/np.sqrt(2)               # |0,1>_X
zero_two_x = (two_zero_z - 2 * one_one_z + zero_two_z)/2        # |0,2>_X
one_one_x  = (two_zero_z - zero_two_z)/2                        # |1,1>_X

# PNS attack
PNS_attack = [ 
    (np.outer(vac, np.conj(vac))),                                                           # Pnot
    (np.outer(one_zero_z, np.conj(one_zero_z)) + np.outer(zero_one_z, np.conj(zero_one_z))), # P0
    (np.outer(two_zero_z, np.conj(two_zero_z)) + np.outer(zero_two_z, np.conj(zero_two_z))), # P1
    (np.outer(one_one_z, np.conj(one_one_z)))                                                # P01
]

# povmb
povmb = [ 
    # np.outer(vac, np.conj(vac)),                                         # |vac><vac|
    np.outer(one_zero_z + two_zero_z, np.conj(one_zero_z + two_zero_z)), # |j,0><j,0|_z
    np.outer(zero_one_z + zero_two_z, np.conj(zero_one_z + zero_two_z)), # |0,j><0,j|_z
    # np.outer(one_one_z, np.conj(one_one_z)),                             # |1,1><1,1|_z
    np.outer(one_zero_x + two_zero_x, np.conj(one_zero_x + two_zero_x)), # |j,0><j,0|_x
    np.outer(zero_one_x + zero_two_x, np.conj(zero_one_x + zero_two_x)), # |0,j><0,j|_x
    # np.outer(one_one_x, np.conj(one_one_x)),                             # |1,1><1,1|_x
]

# qber
qber = np.linspace(0., 0.12, 6)

for mu in [1e-2, 1e-1, 5e-1, 1, 5]:
    alpha = np.sqrt(mu)

    # weak coherent states
    decoy_states = [
        np.exp(-mu/2)*( vac + alpha * one_zero_z + alpha**2/np.sqrt(2) * two_zero_z ), # |alpha; 0>_Z
        np.exp(-mu/2)*( vac + alpha * zero_one_z + alpha**2/np.sqrt(2) * zero_two_z ), # |alpha; 1>_Z
        np.exp(-mu/2)*( vac + alpha * one_zero_x + alpha**2/np.sqrt(2) * two_zero_x ), # |alpha; 0>_X
        np.exp(-mu/2)*( vac + alpha * zero_one_x + alpha**2/np.sqrt(2) * zero_two_x )  # |alpha; 1>_X
    ] 

    # new simulation
    sim = qkd.QKD(
        dim_a=2,
        dim_b=6,
        n_of_signal_states=4, 
        list_states_a=[qkd.zero, qkd.one, qkd.plus, qkd.minus], list_of_prob_a=[0.25, 0.25, 0.25, 0.25],
        list_states_b=decoy_states, list_of_prob_b=[0.25, 0.25, 0.25, 0.25],
        #povm_a=povma,
        povm_b=povmb,
        simulation_name="BB84 decoy"
    )

    th, pr, dl = [], [], []

    for ii in qber:
        print("\nQBER =", ii)
        hp = qkd.binary_entropy(ii)

        sim.apply_eve_strategy(PNS_attack)
        sim.apply_quantum_channel(qkd.depolarizing_channel(2*ii, 6))

        gamma = []
        # simulate measurments
        for jj in sim.povm:
            gamma.append(np.trace(jj @ sim.rho_ab))
        # for jj in sim.orth_set_a:
        #     gamma.append(np.trace(jj @ sim.rho_ab))
        
        # sim.set_constraints(gamma, np.concatenate([sim.povm, sim.orth_set_a]))
        sim.set_constraints(gamma, sim.povm)

        sim.compute_primal()
        sim.compute_dual()

        pr.append((sim.primal_sol) - hp)
        dl.append((sim.dual_sol) - hp)
        print(pr[-1])
        print(dl[-1])

    # plt.plot(qber, pr, "o", alpha=0.5, label="step 1")
    plt.plot(qber, dl, ".-", alpha=0.5, label="mu="+str(mu))
plt.plot(qber, [1-2*qkd.binary_entropy(i) for i in qber], "-", alpha=0.5, linewidth=1.2, label="theoretical")
plt.title("Reliable lower bound for EB BB84 PNS attacks")
plt.xlabel("QBER")
plt.ylabel("Secret key rate")
plt.ylim([0.,None])
plt.legend(loc="best")
plt.grid()
# plt.savefig("analysis/three_states_bb84.png")
plt.show()