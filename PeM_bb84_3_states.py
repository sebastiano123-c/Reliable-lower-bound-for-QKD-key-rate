from src import qkd
import numpy as np
import matplotlib.pyplot as plt

basis = np.eye(4)

# povm a
povma = []
for ii in basis:
    povma.append( np.outer(ii, np.conj(ii)) )

# new simulation
sim = qkd.QKD(dim_a=4, dim_b=2, n_of_signal_states=4,
list_states_a=np.eye(4), list_of_prob_a=[0.25,0.25,0.5,0],
list_states_b=[qkd.zero, qkd.one, qkd.plus, qkd.minus], list_of_prob_b=[0.7/2,0.7/2,0.3/2,0.3/2],
povm_a=povma)

qber = np.linspace(0., 0.12, 15)
th, pr, dl = [], [], []

for ii in qber:
    print("\QBER =", ii)
    hp = qkd.binary_entropy(ii)

    sim.apply_quantum_channel(qkd.depolarizing_channel(2*ii))

    gamma = []
    for jj in sim.povm:
        gamma.append(np.trace(jj @ sim.rho_ab))
    for jj in sim.orth_set_a:
        gamma.append(np.trace(jj @ sim.rho_ab))
    
    sim.set_constraints(gamma, np.concatenate([sim.povm, sim.orth_set_a]))

    sim.compute_primal()
    sim.compute_dual()

    th.append(1-2*hp)
    pr.append(sim.primal_sol - hp)
    dl.append(sim.dual_sol - hp)

plt.plot(qber, th, ".-", alpha=0.5, label="theorical")
plt.plot(qber, th, ".-", alpha=0.5, label="step 1")
plt.plot(qber, th, ".-", alpha=0.5, label="step 2")
plt.title("reliable lower bound for P&M 3 states BB84")
plt.grid()
plt.legend(loc="best")
plt.show()
