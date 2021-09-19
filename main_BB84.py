#---------------------------------------------------------------------
# QKD: RELIABLE LOWER BOUND FOR P&M - BB84 PROTOCOL SIMULATION
#---------------------------------------------------------------------
    #
    #   AUTHOR:
    #       Sebastiano Cocchi
    #
    #   LENGUAGE:
    #       Python 3
    #
    #   DESCRIPTION:
    #      "THIS SCRIPT FINDS A RELIABLE LOWER BOUND FOR THE SECRET 
    #       KEY RATE OF THE BB84 PROTOCOL WITH TWO MUBs IN THE PREPARE 
    #       AND MEASURE (P&M) SCHEME.
    #       TO RECOVER THE ENTANGLED-BASED (EB) SCHEMES SECURITY
    #       PROOF, SOURCE-REPLACEMENT SCHEME IS PERFORMED (A 
    #       BRIEF DESCRIPTION IS PRESENTED HERE BELOW).
    #       LOWER BOUND IS FOUND FIRSTLY CALCULATING THE MINIMUM OF THE 
    #       RELATIVE ENTROPY BETWEEN THE STATE (SHARED BY ALICE AND
    #       BOB) AND THE STATE AFTER THE ALICE HOLDS HIS RAW KEY.
    #       THIS WILL FIND THE WORST CASE SCENARIO IN WHICH EVE HOLD A 
    #       PURIFICATION OF THE QUBIT.
    #       THUS, THE PROBLEM IS CONVERTED INTO A MAXIMIZATION PROBLEM OF 
    #       FORMER THE FUNCTION.
    #       THE CLASS OF THESE MINIMIZATIONS/MAXIMIZATIONS IS SEMIDEFINITE
    #       PROBLEMS (SDP).
    #       THE TWO RESULTS MAY BE VERY CLOSE TO EACH OTHER AND PROVIDE A
    #       RELIABLE LOWER BOUND ON THE KEY RATE."
    #
    #   FURTHER READINGS CAN BE FOUND AT:
    #       https://doi.org/10.22331/q-2018-07-26-77
    #       https://doi.org/10.1103/PhysRevResearch.3.013274
    #
    #   PACKAGES REQUIRED:
    #    * numpy
    #    * scipy
    #    * cvxpy (with solvers CVXOPT, it can be obtained by 'pip install cvxopt')
    #    * matplotlib
    #
    #   SYMBOLS INDEX:
    #    1)  .x. = tensor product
    #    2)  .+. = direct sum
    #    3)  **+ = hermitian conjugate
    #    4)  .T  = matrix transpose
    #    5)  s.t.= such that
    #
#---------------------------------------------------------------------
# EB AND P&M SCHEMES:
#---------------------------------------------------------------------
    # " An entangled state composed of two photons is created.
    #   One particle is given to Alice, one to Bob.
    #   Both perform a measurment in one of the two bases."
    #
    #              |Psi>          
    #   Alice <______|______>   Bob    
    #   
    # Prepare and Measure (P&M) scheme:
    # " Alice prepares a qubit and sends it to Bob."
    #
#---------------------------------------------------------------------
# SOURCE-REPLACEMENT SCHEME for P&M BB84:
#---------------------------------------------------------------------
    # P&M schemes can be seen as EB schemes using the
    # so-called Source-Replacement scheme in which Alice
    # can choose between {|0>,|1>} and {|+>,|->}.
    # She create an entangled state psi_AA', where A is a
    # register in H^4 and A' is the register in which is encoded the bit,
    # so it is a H^2.
    #
    #
    #   SCHEME:
    #
    #                                            ........Eve............
    #   ..............Alice............          : U_E=I_A.0.E_A'(rho) :         .......Bob.......
    #   :               |PHI>_AA'     :          :........../\.........:         :    A' --> B   :
    #   :   ________________/\________:____________________/  \__________________:_______________
    #   :   A:{0,1,2,3}H^2      A':{0,1,+,-}H^2                                  :B:{0,1,+,-}H^2 :
    #   :       dim=4           dim=2 :                                          :     dim=2     :
    #   :        |                    :                                          :       |       :
    #   :        V                    :                                          :       V       :
    #   :      POVM_A                 :                                          :    POVM_B     :
    #   :.............................:                                          :...............:
    #
    #   STATE A    : STATE A'(B') : BASIS CHOICE : BIT VALUE 
    #       |0>    :     |0>      :      Z       :     0
    #       |1>    :     |1>      :      Z       :     1
    #       |2>    :     |+>      :      X       :     0
    #       |3>    :     |->      :      X       :     1
    #
#---------------------------------------------------------------------
# PROCEDURE, CALCULATION AND ALGORITHM
#---------------------------------------------------------------------
    #   PROCEDURE:
    #    1- Alice creates the states:
    #           |psi>_AA' = sum_i p_i |i>_A .o. |phi_i>_A'
    #    2- Alice sends the qubit A' via quantum channel;
    #    3- ANNOUNCEMENT PHASE:
    #       POVMA and POVMB are the POVM of Alice and Bob, by which are
    #       defined Bob announcement
    #            KB_b = \sum_{\beta_b} \sqrt{POVMB^{b,\beta_b}} .o. |b>_{\tilde{B}} .o. |\beta_b>_{\bar{B}}
    #       where \tilde{B} and \bar{B} are the announcement, and the second index denotes a particular
    #       element associated with that announcement.
    #       Alice announcement is
    #            KA_a = \sum_{\alpha_a} \sqrt{POVMA^{a,\alpha_a}} .o. |a>_{\tilde{A}} .o. |\alpha_a>_{\bar{A}}
    #       The action of KA and Kb on the state \rho is
    #           rho_2 = \sum_{a,b} ( KA_a .o. KB_b ).\rho.( KA_a .o. KB_b )**+
    #                 = A(\rho)
    #   4- SIFTING PHASE:
    #       viewed as a projector of the form
    #           Proj = \sum_{a,b} |a><a|_{\tilde{A}} .o. |b><b|_{\tilde{B}} 
    #       and identities elsewhere.
    #   5- KEY MAP:
    #       write the key map as the function g(a, \alpha_a, b). We deﬁne an isometry V that stores the key
    #       information in a register system R, as follows
    #           V = \sum_{a, \alpha_a, b} |g(a, \alpha_a, b)>_R .0. |a><a|_{\tilde{A}} .o. |\alpha_a><\alpha_a|_{\bar{A}} .o. |b><b|_{\tilde{B}}
    #   6- PINCHING CHANNEL:
    #       decohere the register in order to obtain a classical register R
    #           Z(\sigma) = \sum_j (|j><j|_R .o. I ).\sigma.(|j><j|_R .o. I )
    #       where I denotes the identity acting on the other subsystems.
    #   
    #   CALCULATION:
    #    1- construct the state rho_AB;
    #    2- quantum channel acting on rho_AB (like depolarization, etc.);
    #    3- define the CP map as:
    #       G(\rho) = V . proj . A(\rho) . proj . V**+
    #    4- define the pinching channel as
    #       Z(\rho) = \sum_j (|j><j|_R .o. I ).\sigma.(|j><j|_R .o. I )
    #    5- construct the basis of operators
    #           Gamma_i = POVM_i
    #       and find the Gram-Schmidt process for them.
    #    6- enlarge the basis using THeta_j complete set for A s.t.
    #           Tr[ \Theta_j \rho_AB] = \theta_j      
    #    7- calculation of the constraints:
    #           p_i = Tr[POVM_i . \rho_AB]
    #
    #   ALGORITHM:
    #    1- set counter=0, epsilon, maxit and rho_0 as starting point using the Gamma_tilde
    #           rho_0 = \sum_i gamma_tilde_i Gamma_tilde_i  +  \sum_j  omega_j  Omega_j
    #   STEP 1
    #    2- calculate
    #           f(\rho) = D( G(\rho_0) || Z(G(\rho_0)) )
    #       where 
    #           D(\rho||\sigma) = Tr[ \rho. log(\rho) ] - Tr[ \rho. log(\sigma) ]
    #       is the Relative Entropy.
    #    3- calculate the gradient of this function
    #           grad_f(\rho).T = G**+(log(G(\rho))) + G**+(log(Z(G(\rho))))
    #    4- minimize  : \Delta rho = arg min _{\delta rho} [ (\delta rho).T grad_f(\rho_0) ]
    #       subject to: constraints p_i.
    #    IF: Tr[ (\Delta rho).T . grad_f(\rho_0) ] < \epsilon and GOTO STEP 2
    #    ELSE: find tt \in [0, 1] s.t. minimize f( \rho_0 + tt*\Delta rho )
    #    5- encrease the counter by 1 and set \rho_0 = \rho_0 + tt*\Delta rho
    #       for the next iteration and repeat from 2-.
    #    IF: counter == maxit GOTO STEP 2.
    #           
    #   STEP 2
    #    0- from STEP 1 we know \rho and its gradient
    #    1- maximize : gamma_i.y 
    #      subject to: sum_j y_j Gamma_j <= grad_f(\rho)
    #   
    #   LOWER BOUND:
    #    the result is:
    #      f(rho) - Tr( rho . grad_f(rho) ) + max{gamma_i.y }
    #
#---------------------------------------------------------------------
# The program is diveded in two parts:
#   1) explanation of the conceptual steps of the procedure and
#      the declaration of the usefule operators;
#   2) the algorithm procedure implementing the operators defined in 
#      the previous point and SDP minimization procedure.
# N.B.: the algorithm is a recursive iteration incrementing the depolarization
#       probability.
#---------------------------------------------------------------------
# PART 1): definition of the operators
#---------------------------------------------------------------------

import numpy as np
from scipy.linalg import logm#, sqrtm
import cvxpy as cp
from cvxopt import matrix, solvers
from mosek.fusion import *
import matplotlib.pyplot as plt
import time

# CONSTANTS
# parameters
d_a = 4
d_b = 2
epsilon = 1e-10
start, stop, step = 0., .3, 10
maxit = 20
finesse = 10
Pz = 0.5
Prob = [Pz/2, Pz/2, (1. - Pz)/2, (1. - Pz)/2]
solver_name = "MOSEK"
solver_verbosity = False

# define ids
id_2 = np.eye(2)
id_4 = np.eye(4)
id_8 = np.eye(8)
id_64 = np.eye(64)
id_128 = np.eye(128)

# pauli matrices
pauli = [[[1,0],[0,1]], [[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]]

# define qubits 
zero = np.array([1,0])
one = np.array([0,1])
states = [zero, one, (zero+one)/np.sqrt(2.), (zero-one)/np.sqrt(2.)]
basis = [[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]

# alice measurments
sigma_00 = np.outer(basis[0], np.conj(basis[0]))
sigma_11 = np.outer(basis[1], np.conj(basis[1]))
sigma_pp = np.outer(basis[2], np.conj(basis[2]))
sigma_mm = np.outer(basis[3], np.conj(basis[3]))

# bob measurments
tau_00 = np.outer(states[0], np.conj(states[0]))
tau_11 = np.outer(states[1], np.conj(states[1]))
tau_pp = np.outer(states[2], np.conj(states[2]))
tau_mm = np.outer(states[3], np.conj(states[3]))

# FUNCTIONS
def complex_to_real_isometrym(mat):
    re = np.real(mat)
    im = np.imag(mat)
    return 1/2 * np.block( [[re, -im], [im, re]] )

def inv_complex_to_real_isometrym(mat):
    dims = [len(ii) for ii in mat]
    dim = dims[0]
    sep = int(dim/2)
    re = mat[:sep, :sep]
    im = mat[sep:, :sep]
    return 2 * (re + 1j* im)

def gram_schmidtm(V):
    """Orthogonalize a set of matrices stored as elements of V(:,:,:)."""
    n, k, l = np.shape(V)
    U = np.zeros((n, k, l))*1j
    U[0] = V[0] / np.trace(np.conj(V[0]).T @ V[0])
    for ii in range(1, n): # remember the Hilbert-schmidt norm of two operators A, B is NORM = Tr[ A @ B ]
        U[ii] = V[ii]
        for jj in range(0, ii-1):
            U[ii] = U[ii] - (np.trace( np.conj(U[jj]).T @ U[ii]) / np.trace( np.conj(U[jj]).T @ U[jj] ) ) * U[jj]
        if( abs(np.trace( np.conj(U[ii]).T @ U[ii] )) <= 1e-8 ):
            print( "gram_schmidtm:: stopped at ", ii, ". 'norm == inf'" )
            U = np.delete(U, [kk for kk in range(ii,n)], axis=0)
    return U

def extend_basism(V, B):
    """Extend a set of j.leq.k orthogonal matrices and k.leq.m orthonormal matrices to basis of the space H^m.
        ---------------------------------------------
        Keyword arguments:
        V -- the new set of matrices (size (j, n, n) )
        B -- the orthonormal set to be extended (size (k, n, n))
        return C = (B,V^{orth})
    """
    C = B
    for ii in range(np.shape(V)[0]): # remember the Hilbert-schmidt norm of two operators A, B is NORM = Tr[ A @ B ]
        U = V[ii]
        for jj in C:
            U = U - (np.trace( np.conj(jj).T @ U) / np.trace( np.conj(jj).T @ jj ) ) * jj
        if( abs(np.trace( np.conj(U).T @ U )) >= 1e-8 ):
            C = np.append(C, [U], axis=0)
    return C

def relative_entropy(rho, sigma):
    # avlr = np.real(np.linalg.eigvals(rho))
    # avls = np.real(np.linalg.eigvals(sigma))
    # res = 0.
    # for ii in range(len(avlr)):
    #     if(abs(avlr[ii]) >= 1e-10):
    #         if(abs(avls[ii]) <= 1e-10):#if q_i is zero --> -p_i log(0) = +inf
    #             res = + np.inf
    #         else:                   # sum_i p_i log(p_i/q_i)
    #             res = res + avlr[ii] * np.log(avlr[ii]/avls[ii])
    #     #if(res < 0): print("relative_entropy:: is negative.")
    res = np.trace(rho @ (logm(rho) - logm(sigma)))
    return res/np.log(2)

def vonneumann_entropy(rho):
    return -np.trace( rho @ logm(rho) )

def partial_trace(rho, qubit_2_keep):
    """ Calculate the partial trace for qubit system
    Parameters
    ----------
    rho: np.ndarray
        Density matrix
    qubit_2_keep: list
        Index of qubit to be kept after taking the trace (NB: start from 0!)
    Returns
    -------
    rho_res: np.ndarray
        Density matrix after taking partial trace
    """
    num_qubit = int(np.log2(rho.shape[0]))
    qubit_axis = [(i, num_qubit + i) for i in range(num_qubit)
                  if i not in qubit_2_keep]
    minus_factor = [(i, 2 * i) for i in range(len(qubit_axis))]
    minus_qubit_axis = [(q[0] - m[0], q[1] - m[1])
                        for q, m in zip(qubit_axis, minus_factor)]
    rho_res = np.reshape(rho, [2, 2] * num_qubit)
    qubit_left = num_qubit - len(qubit_axis)
    for i, j in minus_qubit_axis:
        rho_res = np.trace(rho_res, axis1=i, axis2=j)
    if qubit_left > 1:
        rho_res = np.reshape(rho_res, [2 ** qubit_left] * 2)
    return rho_res

def binary_entropy(p):
    if p==0: return 0
    elif p==1: return 0
    else: return - p*np.log(p)/np.log(2) - (1-p)*np.log(1-p)/np.log(2)

def sdp_solver(d, rho, grad_f, Gam, gam, solver_name='MOSEK', solver_verbosity=False):

    # minimize: X.T @ Grad_f
    X = cp.Variable((d, d), complex=True)

    # subject to:
    constraints = [] # The operator >> denotes matrix inequality.
    constraints = constraints + [ X + rho >> 0 ] # Pos(H_AB)
    constraints = constraints + [ cp.real(cp.trace( X + rho )) == 1. ] # trace sum up to 1
    constraints = constraints + [ X + rho == cp.conj(X + rho).T ] # is a Hermitian
    for ii, elm in enumerate(Gam):
        constraints = constraints + [ cp.trace( (X + rho) @ elm)  == gam[ii]]

    # solve
    obj = cp.Minimize( cp.real(cp.trace( X.T @  grad_f )) )      
    prob = cp.Problem( obj, constraints )
    try:
        prob.solve(solver=solver_name, verbose=solver_verbosity)
    except:
        print("\n",solver_name + " failed.")
        return np.eye(d)

    # check if there is any problem in the minimization procedure
    if prob.status in ["infeasible", "unbounded"]:# Otherwise, prob.value is inf or -inf, respectively.
        print("Status problem: %s" % prob.status)
        exit()
    
    # solution
    sol = X.value

    # clear variables
    del X, prob, obj, constraints

    return sol

def compute_primal(rho_0, ZA, Gamma, gamma, epsilon = 1e-10, maxit = 20, finesse = 10, solver_name = 'MOSEK', solver_verbosity = False):

    # set 
    counter = 1

    # start algorithm
    while(counter <= maxit):
        d = np.shape(rho_0)[0]
        ZrhoZ = sum( [ ii @ rho_0 @ ii for ii in ZA ] )
        ZrhoZ = ZrhoZ / np.trace(ZrhoZ)
        if( np.abs( np.trace(ZrhoZ) - 1. ) >= 1e-8): print("Tr[ZrhoZ] != 1", np.trace(ZrhoZ))
        if( np.allclose( np.conj(ZrhoZ).T, ZrhoZ) == False ): print("ZrhoZ NOT hermitian", ZrhoZ)
        if( np.all( np.linalg.eigvals(ZrhoZ) < - 1e-8)): print("ZrhoZ is NEGATIVE")

        # gradient
        bb84_frho = np.real(relative_entropy( rho_0, ZrhoZ ))
        bb84_grad = (logm(rho_0) + logm(ZrhoZ)).T/np.sqrt(2)
        print("iteration", counter, " f(rho) =", bb84_frho)

        Delta_rho_0 = sdp_solver(d, rho_0, bb84_grad, Gamma, gamma, solver_name, solver_verbosity)

        if(abs( np.trace( Delta_rho_0.T @ bb84_grad)) <= epsilon):
            print("Algorithm exited at:", counter, "step")
            break
        elif (counter == maxit):
            print("algorithm reached maxit =", maxit)
            break

        # find l\in (0,1) s.t. min f(rho_i + l* delta_rho)
        tt = 0.
        f_1 = bb84_frho
        for ii in np.linspace(0., 1., finesse):
            rho_temp = rho_0 + Delta_rho_0 * ii
            ZrhoZ = sum( [ ii @ rho_temp @ ii for ii in ZA] )
            ZrhoZ = ZrhoZ / np.trace(ZrhoZ)
            f_2 = np.real(relative_entropy(rho_temp, ZrhoZ))
            if (f_2 < f_1):
                tt = ii
                f_1 = f_2

        # if f_1 == f_rho the next step will be the same
        if(abs(f_1-bb84_frho) <= 1e-8): break

        # assign rho_{i+1} for the next iteration
        rho_0 = rho_0 + Delta_rho_0 * tt
        counter = counter + 1

    return bb84_frho, bb84_grad

def compute_dual(f_rho, grad_f, Gamma, gamma, solver_name = 'MOSEK'):
    # maximize: y.gamma
    n = len(gamma)
    Y = cp.Variable(n)
    # subject to: positivity and tr == 1
    dual_constraints = [ cp.real( sum( [Y[ii]*Gamma[ii].T for ii in range(n)] ) - grad_f ) << 0. ] # The operator >> denotes matrix inequality.
    dual_obj = cp.Maximize( cp.real( gamma @ Y ) )
    dual_prob = cp.Problem( dual_obj, dual_constraints)
    dual_prob.solve(solver=solver_name)
    if dual_prob.status in ["infeasible", "unbounded"]:
        # Otherwise, prob.value is inf or -inf, respectively.
        print("Status problem: %s" % dual_prob.status)
    step2_bound = np.real( f_rho - np.trace( rho_0 @ grad_f ) + dual_prob.value)

    return step2_bound

# MAIN
start_time = time.time()

# total dimension
d_tot = d_a * d_b

# probabilities
if (np.sum(Prob) != 1):
    print("Prob != 1")

# Alice POVM
POVMA = [ np.outer(basis[0],np.conj(basis[0])),
          np.outer(basis[1],np.conj(basis[1])),
          np.outer(basis[2],np.conj(basis[2])),
          np.outer(basis[3],np.conj(basis[3]))
]
# which have dimension 4-by-4 and satisfy POVM properties
if( np.allclose( sum([ np.conj(ii).T @ ii for ii in POVMA]), id_4 ) == False): print("sum POVMA**+ POVMA != 1", sum([ np.real(np.conj(ii).T @ ii) for ii in POVMA]) )
for ii in POVMA:
    if(np.allclose( np.conj(ii).T, ii) == False ): print("POVMA NOT hermitian")
    if(np.all( np.linalg.eigvals(ii) < -1e-8)): print("POVMA is NEGATIVE")

# Bob POVM
POVMB = [ 1/np.sqrt(2)*np.outer( states[0], np.conj(states[0]) ),
          1/np.sqrt(2)*np.outer( states[1], np.conj(states[1]) ),
          1/np.sqrt(2)*np.outer( states[2], np.conj(states[2]) ),
          1/np.sqrt(2)*np.outer( states[3], np.conj(states[3]) )
]
# which have dimension 2-by-2 and satisfy POVM properties
if ( np.allclose(np.real(sum([ np.conj(ii).T @ ii for ii in POVMB])), id_2 ) == False): print("sum POVMB**+ POVMB != 1", sum([ np.conj(ii).T @ ii for ii in POVMB]) )
for ii in POVMB:
    if(np.allclose( np.conj(ii).T, ii) == False ): print("POVMB NOT hermitian")
    if(np.all( np.linalg.eigvals(ii) < - 1e-8)): print("POVMB is NEGATIVE")

# The POVM of the entire system
POVM = []
for ii in POVMA:
        for jj in POVMB:
            temp = np.kron( ii, jj )
            POVM.append( temp )
# which have dimension 8-by-8 and satisfy POVM properties
if ( np.allclose(sum([ np.conj(ii).T @ ii for ii in POVM]), id_8 ) == False): print("sum POVM**+ POVM != 1", sum([np.real(np.conj(ii).T @ ii) for ii in POVM]))
for ii in POVM:
    if(np.allclose( np.conj(ii).T, ii) == False ): print("POVM NOT hermitian")
    if(np.all( np.linalg.eigvals(ii) < - 1e-8)): print("POVM is NEGATIVE")

# Theta.o.I_B 
Theta = []
for ii in pauli:
    for jj in pauli:
        Theta.append( np.kron( ii, np.kron( jj, id_2 ) ) )

# Omega_j complete set of variables
Omega = []
for ii in pauli:
    for jj in pauli:
        for kk in pauli:
            Omega.append( np.kron( ii, np.kron( jj, kk ) ) )

# define a complete orthonormal set of operators
Gammat = gram_schmidtm(POVM)
Gammat = extend_basism(Gammat, Theta)
kk = np.shape(Gammat)[0]
Gammat = gram_schmidtm(Gammat)
orthbasis = extend_basism(Gammat, Omega)
jj = np.shape(orthbasis)[0]
Omega = orthbasis[kk:]

# possible outcome measurments
ZA = np.zeros((2, 8, 8))*1j
ZA[0] = np.kron(np.outer(basis[0], np.conj(basis[0])) + np.outer(basis[2], np.conj(basis[2])), id_2) #( |0><0| + |+><+| ) .o. id_2 --> bit 0
ZA[1] = np.kron(np.outer(basis[1], np.conj(basis[1])) + np.outer(basis[3], np.conj(basis[3])), id_2) #( |1><1| + |-><-| ) .o. id_2 --> bit 1

# which have dimension 4-by-4 and satisfy POVM properties
if(np.allclose(sum([ np.conj(ii).T @ ii for ii in ZA]), id_8) == False): print("sum ZA**+ ZA != 1", sum([ np.real(np.conj(ii).T @ ii) for ii in ZA]) )
for ii in ZA:
    if(np.allclose( np.conj(ii).T, ii) == False ): print("ZA NOT hermitian")
    if(np.all( np.linalg.eigvals(ii) < -1e-8)): print("ZA is NEGATIVE")

postselect_prob = np.array([Pz/2, (1-Pz)/2])  # assume they choose pz=px=1/2
Ez = np.kron(sigma_00, tau_11) + np.kron(sigma_11, tau_00)
Ex = np.kron(sigma_pp, tau_mm) + np.kron(sigma_mm, tau_pp)
Gamma = [Ez, Ex]
Gamma = [G*p for G, p in zip(Gamma, postselect_prob)]

# orthogonal basis
Omega = []
for ii in pauli:
    for jj in pauli:
        for kk in pauli:
            Omega.append(np.kron(ii, np.kron(jj, kk)))
Gammat = gram_schmidtm(Gamma)
kk = np.shape(Gammat)[0]

print(np.shape(Gamma))
print(np.shape(Omega))
orth_normal_basis = extend_basism(Gammat, Omega)
Omega = orth_normal_basis[kk:]

# define the state psi_AA'
psi_aa = 0.
for ii in range(d_a):
    psi_aa = psi_aa + np.sqrt(Prob[ii])*np.kron(basis[ii], states[ii] )
if(np.sqrt(psi_aa @ np.conj(psi_aa))-1. >= 1e-8): print("||psi_aa|| ==", np.sqrt(psi_aa @ np.conj(psi_aa)))

# define rho_aa and check if it is physical
rho_aa = np.outer( psi_aa, np.conj(psi_aa) )
if( np.abs( np.trace(rho_aa) - 1.) >= 1e-8): print("Tr[rho_aa] != 1 (", np.trace(rho_aa),")")
if( np.allclose( np.conj(rho_aa).T, rho_aa) == False ): print("rho_aa NOT hermitian")
if( np.all( np.linalg.eigvals(rho_aa) < - 1e-8)): print("rho_aa is NEGATIVE")

# Alice keeps A and send A' which is subjected to the quantum channel
# action, such as the depolarization.
# Here is the check for polarization probabilities going from 0 to 1.
depolarization_probability, simple_step1_bounds, simple_step2_bounds, theoric_bounds = [], [], [], []
for uu in np.linspace(start, stop, step):
    # constructing the 2-by-2 polarization operator
    depo = [ np.sqrt(1-3./4.*uu)*np.array( pauli[0] ),
                 np.sqrt(uu/4.)*np.array( pauli[1] ),
                 np.sqrt(uu/4.)*np.array( pauli[2] ),
                 np.sqrt(uu/4.)*np.array( pauli[3] )
    ]
    if(np.allclose( sum([np.conj(ii).T @ ii for ii in depo]), id_2 ) == False): print("depo is not Kraus")
    
    # depolarization = id_4.o.depo
    rho_ab, sumkt = 0., 0.
    for jj in depo:
        kraus_temp = np.kron( id_4, jj)
        rho_ab = rho_ab + kraus_temp @ rho_aa @ np.conj( kraus_temp ).T
        sumkt = sumkt + np.conj( kraus_temp ).T @ kraus_temp
    # check the depolarization satisfy Kraus representation property
    if( np.allclose( sumkt, id_8) == False): print("Depo4 is not kraus.")
    # normalize
    rho_ab = rho_ab / np.trace(rho_ab)
    # check
    if( np.abs( np.trace(rho_ab) - 1.) >= 1e-8): print("Tr[rho_ab] != 1", np.trace(rho_ab))
    if( np.allclose( np.conj(rho_ab).T, rho_ab) == False): print("rho_ab NOT hermitian")
    if( np.all( np.linalg.eigvals(rho_ab) < - 1e-8) ): print("rho_ab is NEGATIVE")

    # purity
    purity = np.real(np.trace( rho_ab @ rho_ab ))
    print("- - - - - - - - - - - - -")
    print("(depolarization prob.=", uu,". purity of rho_0     =", purity, ")")

# SDP PROBLEM
    # definitions
    counter = 1   # 1) set counter to 0
    rho_0 = 0.    # 2) find rho_0
    for ii in orth_normal_basis:
        rho_0 = rho_0 + np.trace(rho_aa @ ii) * ii
    rho_0 = rho_0 / np.real(np.trace(rho_0))

    # for ii in Gammat:
    #     rho_0 = rho_0 + np.trace(rho_ab @ ii) * ii
    # for ii in Omega:
    #     rho_0 = rho_0 + np.trace(rho_ab @ ii) * ii
    # # normalize
    # rho_0 = rho_0 / np.trace( rho_0 )
    
    # rho_0 = rho_ab

    # #constraints
    # # povm constraints
    # p_j, p_tilde, theta_j, omega_j, gamma_tilde = [], [], [], [], []
    # for ii in POVM:
    #     tmp = np.trace( ii @ rho_ab )
    #     if( abs(np.imag(tmp)) >= 1e-10 ): print("ERROR: complex mean value")
    #     p_j.append( tmp )
    # for ii in Theta:
    #     tmp = np.trace( ii @ rho_ab )
    #     if( abs(np.imag(tmp)) >= 1e-10 ): print("ERROR: complex mean value")
    #     theta_j.append( tmp )

    #gamma , Gamma
    # gam, Gam = np.concatenate((p_j, theta_j)),  np.concatenate((POVM, Theta))
    gamma = np.array([uu, uu]) * postselect_prob

    # check rho_0 is physical
    if( np.abs( np.trace(rho_0) - 1. ) >= 1e-8): print("Tr[rho_0] != 1", np.trace(rho_0))
    if( np.allclose( np.conj(rho_0).T, rho_0) == False ): print("rho_0 NOT hermitian", rho_0)
    if( np.all( np.linalg.eigvals(rho_0) < - 1e-8)): print("rho_0 is NEGATIVE")

#  STEP 1
    # step 1 result
    bb84_frho, bb84_grad = compute_primal(rho_0, ZA, Gamma, gamma, epsilon, maxit, finesse, 'CVXOPT', False)

    # store result
    step1_bound = np.real( bb84_frho ) 
    print("Step 1 result is    = ", step1_bound)
    simple_step1_bounds.append(step1_bound)

#  STEP 2
    # step 2 result
    step2_bound = compute_dual(bb84_frho, bb84_grad, Gamma, gamma, 'MOSEK')

    #store result
    print("Step 1 result is    = ", step2_bound)
    simple_step2_bounds.append(step2_bound)

    # and theoretic values
    depolarization_probability.append(uu)
    thrc = 1 - 2*binary_entropy(uu/2)
    if thrc <= 0: thrc = 0.
    theoric_bounds.append(thrc)

print(" --- ", time.time() - start_time, "s --- ")

# plot
plt.plot(depolarization_probability, simple_step1_bounds, "-.", label="simple bb84 step 1")
# plt.plot(depolarization_probability, simple_step2_bounds, "--", label="simple bb84 step 2")
plt.plot(depolarization_probability, theoric_bounds, "--", label="BB84 th.")
plt.title(" BB84 simulation ")
plt.xlabel("depolarization probability")
plt.ylabel(r'$f(\rho_{AB})$')
plt.legend(loc='best')
plt.grid(True)
plt.show()



# minimize: X.T @ Grad_f
        # re_r = np.real(rho_0)
        # im_r = np.imag(rho_0)
        # rho_cmplx = complex_to_real_isometrym(rho_0)
        # m = d_tot
        # X = cp.Variable((2*m, 2*m))
        # # subject to:
        # constraints = [] # The operator >> denotes matrix inequality.
        # constraints = constraints + [ X + rho_cmplx >> 0 ] # Pos(H_AB)
        # constraints = constraints + [ cp.real(cp.trace( X + rho_cmplx )) == 2.] # trace sum up to 1
        # # constraints = constraints + [ X + rho_cmplx == cp.conj(X + rho_cmplx).T] # is a Hermitian
        # for ii, elm in enumerate(Gam):
        #     constraints = constraints + [ cp.trace( (X + rho_cmplx) @ complex_to_real_isometrym(elm) ) == gam[ii]]

        # # solve
        # obj = cp.Minimize( cp.trace( X.T @  complex_to_real_isometrym(bb84_grad) ) )      
        # prob = cp.Problem( obj, constraints )
        # prob.solve(solver=solver_name, verbose=solver_verbosity)
        # # check if there is any problem in the minimization procedure
        # if prob.status in ["infeasible", "unbounded"]:# Otherwise, prob.value is inf or -inf, respectively.
        #     print("Status problem: %s" % prob.status)
        #     break
        # # the solution X is X.value
        # Delta_rho_0 = inv_complex_to_real_isometrym(X.value)
        

        
# def sdp_solver1(d, rho, grad_f, Gamma, gamma, model_name="CMPLX_SDP"):#using mosek library
#     with Model(model_name) as M:
        
#         # Setting up the variables
#         X = M.variable("X", Domain.inPSDCone(2*d))

#         # Setting up constant coefficient matrices
#         C  = Matrix.dense ( complex_to_real_isometrym( grad_f ) )

#         # Objective
#         M.objective(ObjectiveSense.Minimize, Expr.dot(C, X))

#         # Constraints
#         rho = complex_to_real_isometrym(rho)
#         for ii, val in enumerate(Gamma):
#             M.constraint("c"+str(ii+1), Expr.sum( Expr.mulElm( Expr.add(X, rho), complex_to_real_isometrym(val) ) ), Domain.equalsTo(np.real(gamma[ii])))
#         M.constraint("c0",  Expr.sum( Expr.mulElm( Expr.add(X, rho), np.eye(2*d) ) ), Domain.equalsTo(1))

#         # solve
#         M.solve()

#         print(X.level())


# def sdp_solver(d, rho, grad_f, Gam, gam, solver_name='MOSEK', solver_verbosity=False):

#     # minimize: X.T @ Grad_f
#     rho_dbl = complex_to_real_isometrym(rho)
#     X = cp.Variable((2*d, 2*d))

#     # subject to:
#     constraints = [] # The operator >> denotes matrix inequality.
#     constraints = constraints + [ X + rho_dbl >> 0 ] # Pos(H_AB)
#     constraints = constraints + [ cp.trace( X + rho_dbl ) == 1.] # trace sum up to 1
#     # constraints = constraints + [ X + rho_dbl == cp.conj(X + rho_dbl).T] # is a Hermitian
#     for ii, elm in enumerate(Gam):
#         constraints = constraints + [ cp.trace( (X + rho_dbl) @ complex_to_real_isometrym(elm) ) == gam[ii]]

#     # solve
#     obj = cp.Minimize( cp.trace( X.T @  complex_to_real_isometrym(grad_f) ) )      
#     prob = cp.Problem( obj, constraints )
#     prob.solve(solver=solver_name, verbose=solver_verbosity)

#     # check if there is any problem in the minimization procedure
#     if prob.status in ["infeasible", "unbounded"]:# Otherwise, prob.value is inf or -inf, respectively.
#         print("Status problem: %s" % prob.status)
#         exit()
#     sol = X.value

#     del X, prob, obj, constraints

#     return inv_complex_to_real_isometrym(sol)


# def sdp_solver3(d_tot, rho, grad_f, Omega, Gam, gam):

#     # minimize: C^T @ X = C_1 X_1 + C_2 X_2 + ...
#     # define C as the sequence made of C = { Tr[Omega_1^T @ grad_f_j], Tr[Omega_2 @ grad_f_j], ... }
#     m = np.shape(Omega)[0]
#     c = [np.trace( ii.T @ grad_f ) for ii in Omega]
#     re_c = np.real(c)
#     im_c = np.imag(c)
#     c = matrix( re_c ) # to convert ndarray to array use .tolist()
#     # c = matrix( c, tc='z' )

#     # subject to: sum_i X_i Omega_i >= - rho_i
#     # G stores row-wise the matrices as rows
#     G = []
#     for mat in Omega:
#         tmp = complex_to_real_isometrym(mat)
#         G = np.append(G, tmp.reshape(-1))
#     G = [matrix( G.tolist(), (int(len(G)/m), m), tc='d' )]
#     # G = [matrix( G, tc= 'z' )]
#     # define -rho_i
#     h = complex_to_real_isometrym(rho)
#     h = [ matrix( h.tolist(), (d_tot*2, d_tot*2), tc='d' ) ]
#     # h = [matrix( -rho, tc='z' )]

#     # solve
#     sol = solvers.sdp(c, Gs=G, hs=h)
    
#     # construct Delta_rho
#     Delta_rho = 0.
#     for ii, val in enumerate(sol['x']):
#         Delta_rho = Delta_rho + val * Omega[ii]
    
#     return Delta_rho