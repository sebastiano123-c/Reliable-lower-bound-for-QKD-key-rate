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
    #       LOWER BOUND IS FOUND FIRSTLY CALCULATING THE MINIMUM OF A
    #       SEMIDEFINITE POGRAM (SDP); THUS, THE PROBLEM IS INVERTED,
    #       AND NOW IT A MAXIMIZATION OF A FUNCTION.
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
# Entangled Based (EB) scheme:
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
    # P&M schemes can be seen as EB schemes using the
    # so-called Source-Replacement scheme in which Alice
    # can choose between {|0>,|1>} and {|+>,|->}.
    # She create an entangled state psi_AA', where A is a
    # register in H^4 and A' is the register in which is encoded the bit,
    # so it is a H^2.
    #
    #   QUANTUM CHANNEL: 
    #    * DEPOLARIZING CHANNEL;
    #
#---------------------------------------------------------------------
# SOURCE-REPLACEMENT SCHEME for P&M BB84:
#---------------------------------------------------------------------
    #
    #   SCHEME:
    #
    #                                            ........Eve............
    #   ..............Alice............          : U_E=I_A.0.E_A'(rho) :         .......Bob.......
    #   :               |PHI>_AA'     :          :........../\.........:         :    A' --> B   :
    #   :   ________________/\________:____________________/  \__________________:_______________
    #   :   A:{0,1,+,-}H^2      A':{0,1,+,-}H^2                                  :B:{0,1,+,-}H^2 :
    #   :       dim=2           dim=2 :                                          :     dim=2     :
    #   :        |                    :                                          :       |       :
    #   :        V                    :                                          :       V       :
    #   :      POVM_A                 :                                          :    POVM_B     :
    #   :.............................:                                          :...............:
    #
    #   STATE A    : STATE A'(B') : BASIS CHOICE : BIT VALUE 
    #       |0>    :     |0>      :      Z       :     0
    #       |1>    :     |1>      :      Z       :     1
    #       |+>    :     |+>      :      X       :     0
    #       |->    :     |->      :      X       :     1
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
from numpy.lib.type_check import real
from scipy.linalg import logm
import cvxpy as cp
import matplotlib.pyplot as plt
import time

# CONSTANTS
# parameters
epsilon = 1e-10
Pz = 0.99
start, stop, step = 0, 0.14, 10
maxit = 20
finesse = 2000
solver_name = 'MOSEK'
solver_verbosity = False

# probabilities
Px = (1. - Pz)
postselect_prob = np.array([Pz*Pz, Px*Px]) 
Prob = np.array([Pz, Pz, Px, Px])
Prob = Prob / np.linalg.norm(Prob)

# pauli matrices
pauli = [[[1,0],[0,1]], [[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]]

# define qubits 
zero = np.array([1,0])
one = np.array([0,1])
states = [zero, one, (zero+one)/np.sqrt(2.), (zero-one)/np.sqrt(2.)]

# define ids
id_2 = np.eye(2)
id_4 = np.eye(4)
id_8 = np.eye(8)
id_64 = np.eye(64)
id_128 = np.eye(128)

# local proj
sigma_00 = np.outer(states[0], np.conj(states[0]))
sigma_11 = np.outer(states[1], np.conj(states[1]))
sigma_pp = np.outer(states[2], np.conj(states[2]))
sigma_mm = np.outer(states[3], np.conj(states[3]))

# FUNCTIONS
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
    avlr = np.real(np.linalg.eigvals(rho))
    avls = np.real(np.linalg.eigvals(sigma))
    res = 0.
    for ii in range(len(avlr)):
        if(abs(avlr[ii])>=1e-8):
            if(abs(avls[ii])<=1e-8):#if q_i is zero --> -p_i log(0) = +inf 
                res = + np.inf
            else:                   # sum_i p_i log(p_i/q_i)
                res = res + avlr[ii] * np.log(avlr[ii]/avls[ii])/np.log(2)
        #if(res < 0): print("relative_entropy:: is negative.")
    # res = np.trace(rho @ (logm(rho)/np.log(2) - logm(sigma)/np.log(2)))
    return res

def binary_entropy(p):
    if p==0: return 0
    elif p==1: return 0
    else: return - p*np.log(p)/np.log(2) - (1-p)*np.log(1-p)/np.log(2)

def CP_map(rho, Kraus, sifting, isometry):
    '''G map'''
    rho_temp = 0.
    for ii in Kraus:
        rho_temp = rho_temp + ii @ rho @ np.conj(ii).T
    rho_temp = sifting @ rho_temp @ sifting
    Ppass = np.real(np.trace( rho_temp ))
    rho_temp = rho_temp / Ppass
    rho_temp = isometry @ rho_temp @ np.conj(isometry).T
    #rho_temp = rho_temp / np.real(np.trace(rho_temp))
    return rho_temp

def CP_inverse_map(rho, Kraus, sifting, isometry):
    '''inveser CP map'''
    rho_temp = 0.
    rho_temp = np.conj(isometry).T @ rho @ isometry
    rho_temp = sifting @ rho_temp @ sifting
    rho_fin = 0.
    for ii in Kraus:
        rho_fin = rho_fin + np.conj( ii ).T @ rho_temp @ ii
    #rho_fin = rho_fin / np.real(np.trace(rho_fin))
    return rho_fin

def sdp_solver(d, rho, grad_f, Gamma, gamma, solver_name = 'MOSEK', solver_verbosity = False):

    # minimize: X.T @ Grad_f
    X = cp.Variable((d, d), complex=True)

    # subject to:
    constraints = [] # The operator >> denotes matrix inequality.
    constraints = constraints + [ X + rho >> 0 ] # Pos(H_AB)
    constraints = constraints + [ cp.real(cp.trace( X + rho )) == 1 ] # trace sum up to 1
    constraints = constraints + [ X + rho == cp.conj(X + rho).T ] # is a Hermitian
    for ii, elm in enumerate(Gamma):
        constraints = constraints + [ cp.trace( (X + rho) @ elm)  == gamma[ii]]

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

def compute_primal(rho_0, Kraus, sifting, isometry, ZA, Gamma, gamma, epsilon = 1e-10, maxit = 20, finesse = 10, solver_name = 'MOSEK', solver_verbosity = False):

    # set 
    counter = 1

    # start algorithm
    while(counter <= maxit):
        
        # dimension
        d = np.shape(rho_0)[0]

        # define the two states
        rho_4 = CP_map(rho_0, Kraus, sifting, isometry)
        rho_5 = sum( [ ii @ rho_4 @ ii for ii in ZA ] )
        rho_5 = rho_5 / np.real(np.trace(rho_5))
        if( np.abs( np.trace(rho_5) - 1. ) >= 1e-8): print("Tr[rho_5] != 1", np.trace(rho_5))
        if( np.allclose( np.conj(rho_5).T, rho_5) == False ): print("rho_5 NOT hermitian", rho_5)
        if( np.all( np.linalg.eigvals(rho_5) < - 1e-8)): print("rho_5 is NEGATIVE")

        # f_rho and gradient
        bb84_frho = np.real(relative_entropy( rho_4, rho_5 ))
        bb84_grad = (CP_inverse_map( logm( rho_4 )/np.sqrt(2), Kraus, sifting, isometry) - CP_inverse_map( logm( rho_5 )/np.sqrt(2), Kraus, sifting, isometry ) ).T

        # print iteration result
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

def compute_dual(grad_f, Gamma, gamma, solver_name = 'MOSEK', solver_verbosity = False):

    # maximize: y.gamma
    n = len(gamma)
    Y = cp.Variable(n)

    # subject to: positivity and tr == 1
    dual_constraints = [ sum( [Y[ii]*Gamma[ii] for ii in range(n)] ) <= cp.real(grad_f) ] # The operator >> denotes matrix inequality.
    dual_constraints += [ sum([Y[ii]*Gamma[ii] for ii in range(n)]) >> 0 ]

    # write
    dual_obj = cp.Maximize( cp.real(gamma) @ Y )
    dual_prob = cp.Problem( dual_obj, dual_constraints )

    # solve
    dual_prob.solve(solver=solver_name, verbose=solver_verbosity)

    # check
    if dual_prob.status in ["infeasible", "unbounded"]:
        # Otherwise, prob.value is inf or -inf, respectively.
        print("Status problem: %s" % dual_prob.status)
    
    return dual_prob.value

# MAIN
start_time = time.time()

# possible outcome measurments
ZA = np.zeros((2, 4, 4))*1j
ZA[0] = np.kron( sigma_00, id_2) #( |0><0| + |+><+| ) .o. id_2 --> bit 0
ZA[1] = np.kron( sigma_11, id_2) #( |1><1| + |-><-| ) .o. id_2 --> bit 1

# normalization
ZA = ZA / sum([np.conj(ii).T @ ii for ii in ZA])[0][0]

# which have dimension 4-by-4 and satisfy POVM properties
if ( np.allclose( sum([ np.conj(ii).T @ ii for ii in ZA] ), id_4 ) == False): print("sum ZA**+ ZA != 1", sum([ np.conj(ii).T @ ii for ii in ZA]) )
for ii in ZA:
    if(np.allclose( np.conj(ii).T, ii) == False ): print("ZA NOT hermitian")
    if(np.all( np.linalg.eigvals(ii) < -1e-8)): print("ZA is NEGATIVE")

# Gamma for constraints
Ez = np.kron(sigma_00, sigma_11) + np.kron(sigma_11, sigma_00)
Ex = np.kron(sigma_pp, sigma_mm) + np.kron(sigma_mm, sigma_pp)
Gamma = [Ez, Ex]
Gamma = [G*p for G, p in zip(Gamma, postselect_prob)]

# orthogonal set of observables
Omega = []
for ii in pauli:
    for jj in pauli:
        Omega.append(np.kron(ii, jj))
Gammat = gram_schmidtm(Gamma)
orth_normal_basis = extend_basism(Gammat, Omega)

# Alice prepare psi_AA' \in H^8
psi_aa = 0.
for ii, elm in enumerate(states):
    psi_aa = psi_aa + Prob[ii]*np.kron( elm, elm )
psi_aa = psi_aa / np.sqrt(np.dot(np.conj(psi_aa), psi_aa))
# define rho_aa and check if it is physical
rho_aa = np.outer( psi_aa, np.conj(psi_aa) )
if( np.abs( np.trace(rho_aa) - 1.) >= 1e-8): print("Tr[rho_aa] != 1 (", np.trace(rho_aa),")")
if( np.allclose( np.conj(rho_aa).T, rho_aa)  == False ): print("rho_aa NOT hermitian")
if( np.all( np.linalg.eigvals(rho_aa) < - 1e-8)): print("rho_aa is NEGATIVE")

# Alice keeps A and send A' which is subjected to the quantum channel
# action, such as the depolarization.
# Here is the check for polarization probabilities going from 0 to 1.
th_key, step1_key, step2_key = [], [], []

p_pass = np.sum(postselect_prob)
QBER = np.linspace(start, stop, step)

for qber in QBER:
    print("\n [ QBER =", qber, "]")
# SDP PROBLEM

    # set
    counter = 1   # set counter to 0
    rho_0 = 0.*1j # set rho_0
    for ii in orth_normal_basis:
        rho_0 = rho_0 + np.trace(ii@rho_aa)*ii
    rho_0 = rho_0 / np.real(np.trace(rho_0))

    # constraints
    gamma = np.array([qber, qber]) * postselect_prob

    # binary entropy
    hp = binary_entropy(qber)

#  STEP 1: primal problem

    # step 1 result
    bb84_frho, bb84_grad = compute_primal(rho_0, [id_4], id_4, id_4, ZA, Gamma, gamma, epsilon, maxit, finesse, 'CVXOPT', solver_verbosity)
    step1_bound = p_pass*(np.real( bb84_frho ) - hp) #leak_EC

    # store result
    print("Step 1 result is    = ", step1_bound)
    if step1_bound < 0: step1_bound = 0.
    step1_key.append(step1_bound)

#  STEP 2: dual problem

    # step 2 result
    step2_bound = compute_dual(bb84_grad, Gamma, gamma, 'MOSEK')
    step2_bound = np.real( bb84_frho - np.trace( rho_0 @ bb84_grad ) + step2_bound)

    #store result
    print("Step 1 result is    = ", step2_bound)
    step2_key.append(step2_bound)

#  theoretic values
    thrc = p_pass*(1 - 2*hp)
    if thrc <= 0: thrc = 0.
    th_key.append(thrc)


print(" --- ", time.time() - start_time, "s --- ")

# plot
plt.plot(QBER, step1_key, "o", label="step 1")
# plt.plot(QBER, step2_key, "o", label="simple bb84 step 2")
plt.plot(QBER, th_key, "--", label="theory")
plt.title(" BB84 simulation ")
plt.xlabel("QBER")
plt.ylabel('secret key rate')#(r'$f(\rho_{AB})$')
plt.legend(loc='best')
plt.grid(True)
plt.show()
