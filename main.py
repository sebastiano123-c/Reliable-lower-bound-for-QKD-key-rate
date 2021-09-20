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
from scipy.linalg import logm, sqrtm
import cvxpy as cp
import matplotlib.pyplot as plt
import time

np.seterr('raise')
start_time = time.time()

# parameters
epsilon = 1e-10
Pz = 0.5
start, stop, step = 0, 0.15, 20
maxit = 20
finesse = 10

# pauli matrices
pauli = [[[1,0],[0,1]], [[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]]

# define qubits 
zero = np.array([1,0])
one = np.array([0,1])
states = [zero, one, (zero+one)/np.sqrt(2.), (zero-one)/np.sqrt(2.)]

# probabilities
Px = (1. - Pz)
postselect_prob = np.array([Pz/2, Px/2]) 
Prob = [Pz/2., Pz/2., Px/2., Px/2]
if (np.sum(Prob) != 1): print("Prob != 1")

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

# functions
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
    #     if(abs(avlr[ii])>=1e-10):
    #         if(abs(avls[ii])>=1e-10):#if q_i is zero --> -p_i log(0) = +inf 
    #             # sum_i p_i log(p_i/q_i)
    #             try:
    #                 res = res + avlr[ii] * np.log(avlr[ii]/avls[ii])/np.log(2)
    #             except:
    #                 print(avlr[ii], avls[ii])
    #         else:
    #             print("Q=0 notimply P=0")
    res = np.trace(rho @ (logm(rho)/np.log(2) - logm(sigma)/np.log(2)))
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

def sdp_solver(d, rho, grad_f, Gamma, gamma, solver_name='MOSEK', solver_verbosity=False):

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
            rho_4 = CP_map(rho_0 + Delta_rho_0 * ii, Kraus, sifting, isometry)
            rho_5 = sum( [ ii @ rho_4 @ ii for ii in ZA] )
            rho_5 = rho_5 / np.trace(rho_5)
            f_2 = np.real(relative_entropy(rho_4, rho_5))
            if (f_2 <= f_1):
                tt = ii
                f_1 = f_2

        # if f_1 == f_rho the next step will be the same
        if(abs(f_1-bb84_frho) <= 1e-10): break

        # assign rho_{i+1} for the next iteration
        rho_0 = rho_0 + Delta_rho_0 * tt
        counter = counter + 1

    return np.real(bb84_frho), bb84_grad

def compute_dual(grad_f, Gamma, gamma, solver_name = 'MOSEK', solver_verbosity = False):

    # maximize: y.gamma
    n = len(gamma)
    Y = cp.Variable(n)

    # subject to: positivity and tr == 1
    dual_constraints = [ sum( [Y[ii]*Gamma[ii] for ii in range(n)] ) << cp.real(grad_f) ] # The operator >> denotes matrix inequality.

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

#---------------------------------------------------------------------
#   for SIMPLE BB84
#---------------------------------------------------------------------
# possible outcome measurments
ZA = np.zeros((2, 4, 4))*1j
ZA[0] = np.kron( sigma_00, id_2) #( |0><0| + |+><+| ) .o. id_2 --> bit 0
ZA[1] = np.kron( sigma_11, id_2) #( |1><1| + |-><-| ) .o. id_2 --> bit 1

# which have dimension 4-by-4 and satisfy POVM properties
if ( np.allclose( sum([ np.conj(ii).T @ ii for ii in ZA] ), id_4 ) == False): print("sum POVMA**+ POVMA != 1", sum([ np.conj(ii).T @ ii for ii in ZA]) )
for ii in ZA:
    if(np.allclose( np.conj(ii).T, ii) == False ): print("POVMA NOT hermitian")
    if(np.all( np.linalg.eigvals(ii) < -1e-8)): print("POVMA is NEGATIVE")

# Gamma for constraints
Ez = np.kron(sigma_00, sigma_11) + np.kron(sigma_11, sigma_00)
Ex = np.kron(sigma_pp, sigma_mm) + np.kron(sigma_mm, sigma_pp)
Gamma = [Ez, Ex]
Gamma = [G*p for G, p in zip(Gamma, postselect_prob)]

#---------------------------------------------------------------------
#   for RECONCILIATION SCHEME
#---------------------------------------------------------------------
# After the qubit sending, Alice can measure A using the POVM 
POVMA = [
    1/np.sqrt(2)*np.outer( states[0], np.conj(states[0]) ),
    1/np.sqrt(2)*np.outer( states[1], np.conj(states[1]) ),
    1/np.sqrt(2)*np.outer( states[2], np.conj(states[2]) ),
    1/np.sqrt(2)*np.outer( states[3], np.conj(states[3]) )
]
# which have dimension 4-by-4 and satisfy POVM properties
if ( np.allclose(sum([ np.conj(ii).T @ ii for ii in POVMA]), id_2 ) == False ): print("sum POVMA**+ POVMA != 1", sum([ np.conj(ii).T @ ii for ii in POVMA]) )
for ii in POVMA:
    if(np.allclose( np.conj(ii).T, ii) == False ): print("POVMA NOT hermitian")
    if(np.all( np.linalg.eigvals(ii) < -1e-8)): print("POVMA is NEGATIVE")
# On the other hand, Bob can measure using the POVM
POVMB = [
    1/np.sqrt(2)*np.outer( states[0], np.conj( states[0] ) ),
    1/np.sqrt(2)*np.outer( states[1], np.conj( states[1] ) ),
    1/np.sqrt(2)*np.outer( states[2], np.conj( states[2] ) ),
    1/np.sqrt(2)*np.outer( states[3], np.conj( states[3] ) )
]
# which have dimension 2-by-2 and satisfy POVM properties
if ( np.allclose(sum([ np.conj(ii).T @ ii for ii in POVMB]), id_2 ) == False ): print("sum POVMB**+ POVMB != 1", sum([ np.conj(np.transpose(ii)) @ ii for ii in POVMB]) )
for ii in POVMB:
    if(np.allclose( np.conj(ii).T, ii) == False ): print("POVMB NOT hermitian")
    if(np.all( np.linalg.eigvals(ii) < - 1e-8)): print("POVMB is NEGATIVE")

# PUBLIC ANNOUNCEMENT:
#   kraus operators of A dim = 16-by-4
KA = [ np.kron( sqrtm(POVMA[0]) , np.kron( zero[:, np.newaxis], zero[:, np.newaxis])) + 
       np.kron( sqrtm(POVMA[1]) , np.kron( zero[:, np.newaxis], one[:, np.newaxis])),
       np.kron( sqrtm(POVMA[2]) , np.kron( one[:, np.newaxis], zero[:, np.newaxis])) + 
       np.kron( sqrtm(POVMA[3]) , np.kron( one[:, np.newaxis], one[:, np.newaxis]))
]
#   which satisfy Kraus property
if ( np.allclose( np.sum([ np.conj(ii).T @ ii for ii in KA]), id_2) ): print("sum KA**+ KA != 1", sum([ np.conj(ii).T @ ii for ii in KA]) )
#   kraus operators of B dim = 8-by-4
KB = [ np.kron( sqrtm(POVMB[0]) , np.kron( zero[:, np.newaxis], zero[:, np.newaxis])) + 
       np.kron( sqrtm(POVMB[1]) , np.kron( zero[:, np.newaxis], one[:, np.newaxis])),
       np.kron( sqrtm(POVMB[2]) , np.kron( one[:, np.newaxis], zero[:, np.newaxis])) + 
       np.kron( sqrtm(POVMB[3]) , np.kron( one[:, np.newaxis], one[:, np.newaxis]))
]
#   which satisfy Kraus property
if ( np.allclose(np.sum([ np.conj(ii).T @ ii for ii in KB]), id_2 ) ): print("sum KB**+ KB != 1", sum([ np.conj(ii).T @ ii for ii in KB]) )
#   The total Kraus representation of the Public Announcement is
K = []
for ii in KA:
        for jj in KB:
            K.append( np.kron(ii, jj))
#   which satisfy Kraus property
if ( np.allclose(np.sum([ np.conj(ii).T @ ii for ii in K]), id_4 ) ): print("sum K**+ K != 1", sum([ np.conj(ii).T @ ii for ii in K]) )

# SIFTING PHASE:
k0b0 = np.outer([1,0],[1,0]) # |0><0|
k1b1 = np.outer([0,1],[0,1]) # |1><1|
#   acts like a projector with dimension 128-by-128
proj = np.kron( id_2, np.kron( k0b0, np.kron( id_4, np.kron( k0b0, id_2  ))) ) +\
       np.kron( id_2, np.kron( k1b1, np.kron( id_4, np.kron( k1b1, id_2  ))) )
       
# KEY MAP: 
#   is a isometry which creates a new register R which stores the information on the bits
#   and it is a 258-by-258 matrix
# V = np.kron( zero[:, np.newaxis], np.kron( id_4, np.kron( k0b0, np.kron( k0b0, np.kron( id_2, np.kron( k0b0, id_2 ) ) ) ) ) ) +\
#     np.kron(  one[:, np.newaxis], np.kron( id_4, np.kron( k0b0, np.kron( k1b1, np.kron( id_2, np.kron( k0b0, id_2 ) ) ) ) ) )
V = np.kron( zero[:, np.newaxis], np.kron( id_4, np.kron( k0b0, id_8))) +\
    np.kron(  one[:, np.newaxis], np.kron( id_4, np.kron( k1b1, id_8)))

# PINCHING CHANNEL:
#   decohere the register R. It has the effect of making R a classical register
pinching = [ np.kron( k0b0 , id_64),
             np.kron( k1b1 , id_64)]

# The POVM of the entire system is given by
POVM = []
for ii in POVMA:
        for jj in POVMB:
            temp = np.kron( ii, jj )
            POVM.append( temp )
# which have dimension 8-by-8 and satisfy POVM properties
if ( np.all(sum([ np.conj(ii).T @ ii for ii in POVM]) != id_4 ) ): print("sum POVM**+ POVM != 1",[np.conj(ii).T @ ii for ii in POVM])
for ii in POVM:
    if(np.allclose( np.conj(ii).T, ii) == False ): print("POVM NOT hermitian")
    if(np.all( np.linalg.eigvals(ii) < - 1e-8)): print("POVM is NEGATIVE")
# to fully define the contraints we need also
Theta = []
# pauli are a tomographically complete set of A
for ii in pauli:
    Theta.append( np.kron( ii, id_2) )

# FIND HERM(H_AB) basis
# contruct the set of orthogonal elements in respectot the Hilbert-Schmidt norm
orth_set = []
orth_set.append(POVM[0])
for ii, elm in enumerate(POVM):
    flg = 0#set flag
    for jj, val in enumerate(orth_set):
        trn = np.trace(np.conj(elm).T @ val)
        if(abs(trn)>=1e-8):#if it is NOT null elm and val are NOT orthogonal
            flg = 1
    if(flg == 0):#if elm is orthogonal to all
        orth_set=np.append(orth_set, [elm], axis=0)
# Gamma = (POVM, Theta)
orth_set = extend_basism(Theta, orth_set)

# Gram-Schmidt decomoposition
POVM_tilde = gram_schmidtm(orth_set)
# check if satisfy POVM_tilde properties
if ( np.allclose(np.sum([ np.conj(ii).T @ ii for ii in POVM_tilde]), id_4 ) ): print("sum POVM_tilde**+ POVM_tilde != 1",[np.conj(ii).T @ ii for ii in POVM_tilde])
for ii in POVM_tilde:
    if(np.allclose( np.conj(ii).T, ii) == False ): print("POVM_tilde NOT hermitian")
    if(np.all( np.linalg.eigvals(ii) < -1e-8)): print("POVM_tilde is NEGATIVE")

# Set used to extend POVM_tilde to a basis
Omega = []
for ii in range(4):
    for jj in range(4):
        M = 0.5 * np.kron(pauli[ii], pauli[jj])
        #check for hermiticity
        if(np.allclose(M, np.conj(M).T) == False): print("Gamma", ii, jj, "not hermitian.")
        #check if is a Tr[G_mu G_mu]==1
        if( abs(np.trace(M @ M) - 1.) >= 1e-8): print("Tr[G_mu G_mu] != 1", np.trace(M @ M))
        Omega.append( M )

#Gamma_tilde = {Gamma_tilde_k} U {Omega_j}
kk = np.shape(POVM_tilde)[0]
Gamma_tilde_k = POVM_tilde
Gamma_tilde = extend_basism(Omega, POVM_tilde)
jj = np.shape(Gamma_tilde)[0] - kk
Omega = Gamma_tilde[kk:]
#check that the elements of the basis are orthogonal
for ii, elm in enumerate(Gamma_tilde):
    for jj, val in enumerate(Gamma_tilde):
        if (ii!=jj):
            if(np.trace(np.conj(elm).T @ val )):
                print("MAIN:: element", ii, "and", jj, "are NOT orthogonal.")
#check if the basis contains 16 elements
if(np.shape(Gamma_tilde)[0] != np.shape(Gamma_tilde)[1]*np.shape(Gamma_tilde)[2]):
    print("MAIN:: the basis contains", np.shape(Gamma_tilde)[0], "!=", np.shape(Gamma_tilde)[1]*np.shape(Gamma_tilde)[2])

#---------------------------------------------------------------------
# PART 2): algorithm
#---------------------------------------------------------------------

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
qber = np.linspace(start, stop, step)
step1_bounds, step2_bounds, theoric_bounds = [], [], []
simple_step1_bounds, simple_step2_bounds = [], []
for uu in qber:
    print("\n[ QBER =", uu)
    # constructing the 2-by-2 polarization operator
    depo_prob = 2*uu
    depo_sq = [ np.sqrt(1-3/4*depo_prob)*np.array( pauli[0] ),
                    np.sqrt(depo_prob/4)*np.array( pauli[1] ),
                    np.sqrt(depo_prob/4)*np.array( pauli[2] ),
                    np.sqrt(depo_prob/4)*np.array( pauli[3] )
    ]
    rho_ab, sumkt = 0., 0.
    for jj in depo_sq:
        kraus_temp = np.kron( id_2, jj)
        rho_ab = rho_ab + kraus_temp @ rho_aa @ np.conj( kraus_temp ).T
        sumkt = sumkt + np.conj( kraus_temp ).T @ kraus_temp
    # check the depolarization satisfy Kraus representation property
    if( np.allclose( sumkt, id_4) == False): print("Depo4 is not kraus.")
    rho_ab = rho_ab / np.trace(rho_ab)
    if( np.abs( np.trace(rho_ab) - 1.) >= 1e-8): print("Tr[rho_ab] != 1", np.trace(rho_ab))
    if( np.allclose( np.conj(rho_ab).T, rho_ab) == False): print("rho_ab NOT hermitian")
    if( np.all( np.linalg.eigvals(rho_ab) < - 1e-8) ): print("rho_ab is NEGATIVE")
    print("purity of rho_0     =", np.real(np.trace( rho_ab @ rho_ab )))

    #constraints
    # povm constraints
    p_j, p_tilde, theta_j, omega_j, gamma_tilde = [], [], [], [], []
    for ii in POVM :
        p_j.append( np.trace( ii @ rho_ab ) )
    for ii in range(np.shape( Theta )[0]):
        theta_j.append( np.trace( Theta[ii] @ rho_ab ) )
    # mean values of operators forming rho_ab
    for ii in POVM_tilde:
        p_tilde.append(np.real(np.trace( ii @ rho_ab)))
    for ii in Omega:
        omega_j.append(np.real(np.trace( ii @ rho_ab)))
    for ii in Gamma_tilde:
        gamma_tilde.append(np.trace( ii @ rho_ab))     

    #gamma , Gamma
    gam, Gam = np.concatenate((p_j, theta_j)),  np.concatenate((POVM, Theta))

    #### SDP PROBLEM ###
    counter = 1   # set counter to 0
    #find rho_0
    # rho_0 = 0.    
    # #gamma tilde
    # for ii, val in enumerate( POVM_tilde ):
    #     rho_0 = rho_0 + p_tilde[ii]*val
    # #omega
    
    # rho_0 = rho_0 / np.trace( rho_0 )
    # rho_ab = rho_0
    rho_0 = rho_ab
    if( np.abs( np.trace(rho_0) - 1. ) >= 1e-8): print("Tr[rho_0] != 1", np.trace(rho_0))
    if( np.allclose( np.conj(rho_0).T, rho_0) == False ): print("rho_0 NOT hermitian", rho_0)
    if( np.all( np.linalg.eigvals(rho_0) < - 1e-8)): print("rho_0 is NEGATIVE")

#---------------------------------------------------------------------
#   for SIMPLE BB84
#---------------------------------------------------------------------
    print("*** simple BB84")
    rho_0 = rho_ab
    gamma = np.array([uu, uu]) * postselect_prob
    hp = binary_entropy(uu)
#  STEP 1
    f_rho, grad_f = compute_primal(rho_0, [id_4], id_4, id_4, ZA, Gamma, gamma, epsilon, maxit, finesse, 'MOSEK')
    step1_bound = np.real( f_rho ) - hp
    if step1_bound<0: step1_bound = 0.
    simple_step1_bounds.append(step1_bound)    
    print("Step 1 result is    = ", step1_bound)
#  STEP 2: lower bound
    step2_bound = compute_dual(grad_f, Gam, gam)
    step2_bound = np.real( f_rho - np.trace( rho_0 @ grad_f ) + step2_bound) - hp
    if step2_bound<0: step2_bound = 0.
    print("Step 2 result is    = ", step2_bound)
    simple_step2_bounds.append(step2_bound)
    th = 1-2*hp
    if th<0: th = 0.
    theoric_bounds.append(th)
    
#---------------------------------------------------------------------
#   for RECONCILIATION SCHEME
#---------------------------------------------------------------------
    print("*** BB84 with reconc.")
    counter = 1
    rho_0 = rho_ab
    hp = binary_entropy(uu)
#  STEP 1
    f_rho, grad_f = compute_primal(rho_0, K, proj, V, pinching, Gam, gam, epsilon, maxit, finesse, 'MOSEK')
    step1_bound = np.real( f_rho - hp)
    if step1_bound<0: step1_bound = 0.
    step1_bounds.append(step1_bound)
    print("Step 1 result is    = ", step1_bound)
#  STEP 2: lower bound
    step2_bound = compute_dual(grad_f, Gam, gam)
    step2_bound = np.real( f_rho - np.trace( rho_0 @ grad_f ) + step2_bound) - hp
    if step2_bound<0: step2_bound = 0.
    print("Step 2 result is    = ", step2_bound)
    step2_bounds.append(step2_bound)

print(" --- ", time.time() - start_time, "s --- ")

# plot
plt.plot(qber, simple_step1_bounds, "*", label="simple bb84 step 1")
plt.plot(qber, simple_step2_bounds, "+", label="simple bb84 step 2")
plt.plot(qber, step1_bounds, "*", label="bb84 step 1")
plt.plot(qber, step2_bounds, "+", label="bb84 step 2")
plt.plot(qber, theoric_bounds, "--", label="theoric")
plt.title(" Reliable key rate lower bound for BB84 (Pz="+str(round(Pz,3))+")")
plt.xlabel("QBER")
plt.ylabel("asymptotic key rate")
plt.legend(loc='best')
plt.grid(True)
plt.show()