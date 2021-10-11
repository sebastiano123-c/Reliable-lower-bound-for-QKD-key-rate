"""
@author: Sebastiano Cocchi
"""
import numpy as np
from scipy.linalg import logm, block_diag, sqrtm
import cvxpy as cp

# pauli matrices
pauli    = [
    [[1,0],[0,1]],
    [[0,1],[1,0]], 
    [[0,-1j],[1j,0]], 
    [[1,0],[0,-1]]
]
pauli_4D = [
    [[1,0,0], [0,1,0], [0,0,1]],
    [[0,1,0], [1,0,0], [0,0,0]],
    [[0,-1j,0], [1j,0,0], [0,0,0]],
    [[1,0,0], [0,-1,0], [0,0,0]],
    [[0,0,1], [0,0,0], [1,0,0]],
    [[0,0,-1j], [0,0,0], [1j,0,0]],
    [[0,0,0], [0,0,1], [0,1,0]],
    [[0,0,0], [0,0,-1j], [0,1j,0]],
    [[1,0,0], [0,1,0], [0,0,-2]]/np.sqrt(3.)
]

# qubits
zero   = np.array([1,0])
one    = np.array([0,1])
plus   = (zero+one)/np.sqrt(2.)
minus  = (zero-one)/np.sqrt(2.)
plusi  = (zero+1j*one)/np.sqrt(2.)
minusi = (zero-1j*one)/np.sqrt(2.)

# local measurements
k0b0   = np.outer( zero  , np.conj(zero)   ) # |0><0|
k1b1   = np.outer( one   , np.conj(one)    ) # |1><1|
kpbp   = np.outer( plus  , np.conj(plus)   ) # |+><+|
kmbm   = np.outer( minus , np.conj(minus)  ) # |-><-|
kpibpi = np.outer( plusi , np.conj(plusi)  ) # |+i><+i|
kmibmi = np.outer( minusi, np.conj(minusi) ) # |-i><-i|

def is_physical(rho, tol=1e-10, name="rho"):
    """
    check if a density matrix is physical
    """
    if( np.abs( np.trace(rho) - 1.) >= tol): print("is_physical:: Tr["+name+"] != 1", np.trace(rho))
    if( np.allclose( np.conj(rho).T, rho) == False): print("is_physical:: "+name+" NOT hermitian")
    if( np.all( np.linalg.eigvals(rho) < - tol) ): print("is_physical:: "+name+" is NEGATIVE")

def is_povm(povm, name="povm"):
    if ( np.allclose(sum([ii for ii in povm]), np.eye(len(povm[1])) ) == False ): print("is_povm::sum "+name+"**+ "+name+" != 1", sum([ ii for ii in povm]) )
    for ii in povm:
        if(np.allclose( np.conj(ii).T, ii) == False ): print("is_povm:: "+name+" NOT hermitian")
        if(np.all( np.linalg.eigvals(ii) < - 1e-8)): print("is_povm:: "+name+" is NEGATIVE")

def is_kraus(kraus, name="kraus"):
    if ( np.allclose( np.sum([ np.conj(ii).T @ ii for ii in kraus]), np.eye(len(kraus[1]))) ): print("sum "+name+"**+ "+name+" != 1", sum([ np.conj(ii).T @ ii for ii in kraus]) )

def get_complete_matrix_basis(space_dimension):

    # initialize with the lower dimensio
    if space_dimension in [2,4,8,16]:
        new_dim = 2
        basis = pauli
    else:
        new_dim = 3
        basis = pauli_4D
   
    while(new_dim < space_dimension):

        new_dim   = new_dim*2  # space dimension increase is doubled
        n_elemnt  = new_dim**2 # # of elemnts of a matrix basis is the square of the space dimension

        # istantiate a new basis
        temp = 1j*np.zeros([n_elemnt, new_dim, new_dim])

        # create new basis with dimension doubled
        for ii in range(len(basis)):
            for jj in range(len(pauli)): 

                # kron product
                tmp = np.kron(basis[ii], pauli[jj])
                if( abs(np.trace(tmp @ tmp)) !=0 ): # normalization)
                    tmp = tmp /  np.sqrt(np.trace(tmp @ tmp)) # normalization

                #check for hermiticity
                if(np.allclose(tmp, np.conj(tmp).T) == False): print("tmp ", ii, jj, "not hermitian.")

                #check if is a Tr[G_mu G_mu]==1
                if( abs(np.trace(tmp @ tmp) - 1.) >= 1e-8): print("Tr[tmp tmp] != 1", np.trace(tmp @ tmp))

                temp[jj+(ii)*len(pauli)] = tmp

        basis = temp

    # check orth
    for ii, elm in enumerate(basis):
        for jj, val in enumerate(basis):
            if (ii != jj):
                t = np.trace(np.conj(elm).T @ val)
                if(t != 0):
                    print("Not orthogonal ", ii, jj, t)

    return basis

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
    if(np.shape(V)[1]!=np.shape(B)[1]):
        print("extend_basism:: size V=[", np.shape(V)[1], ",", np.shape(V)[2], "] is different from size B=[", np.shape(B)[1], ",", np.shape(B)[2], "].")
        exit()

    C = B
    for ii in range(np.shape(V)[0]): # remember the Hilbert-schmidt norm of two operators A, B is NORM = Tr[ A @ B ]
        U = V[ii]
        for jj in C:
            U = U - (np.trace( np.conj(jj).T @ U) / np.trace( np.conj(jj).T @ jj ) ) * jj
        if( abs(np.trace( np.conj(U).T @ U )) >= 1e-8 ):
            C = np.append(C, [U], axis=0)
    return C

def binary_entropy(p):
    '''
    computes the binary entropy of a given probabilty p
    ---------------------------------------------------
    args:
    p = probability
    return -p*log(p)/log(2)-(1-p)*log(1-p)/log(2)
    '''
    if p==0: return 0
    elif p==1: return 0
    else: return - p*np.log(p)/np.log(2) - (1-p)*np.log(1-p)/np.log(2)

def von_neumann_entropy(rho):
    """
    Computes the von Neumann entropy of a given quantum state.
    --------------------------------------
    args:
    rho: density matrix
    return tr[rho@log[rho]]
    """
    fudge = 1e-16
    sh = rho.shape
    new_rho = (1 - fudge) * rho + fudge * np.eye(sh[0])
    return -1 * np.trace(np.matmul(rho, logm(new_rho)))

def relative_entropy(rho, sigma):
    '''
    relative_entropy: computes D(a||b) = a@log(a)-a@log(b)
    ---------------------------------------------
    Keyword arguments:
    rho: first density matrix
    sigma: second density matrix
    return D(rho||sigma)
    '''

    # tr[ rho @ logm(rho)]
    res = - von_neumann_entropy(rho)

    # tr[ rho @ logm(sigma)]
    fudge = 1e-16
    sh = sigma.shape

    # to avoid the calculation of of the log a singular matrix
    new_rho = (1 - fudge) * sigma + fudge * np.eye(sh[0])
    res = res - np.real(np.trace(rho @ logm(new_rho)))
    return res/np.log(2)

def depolarizing_channel(p, dim=2):
    """
    returns the Kraus representation of the depolaring channel with probability p
    """
    # initialize with the lower dimensio
    if dim == 2:
        new_dim = 2
        depo = [
            np.sqrt(1-3/4*p)*np.array( pauli[0] ),
            np.sqrt(p/4)*np.array( pauli[1] ),
            np.sqrt(p/4)*np.array( pauli[2] ),
            np.sqrt(p/4)*np.array( pauli[3] )
        ]
    elif dim == 3:
        new_dim = 3
        depo = [ 
            np.sqrt(1-8/9*p)*np.array( pauli_4D[0] ),
            np.sqrt(p/9)*np.array( pauli_4D[1] ),
            np.sqrt(p/9)*np.array( pauli_4D[2] ),
            np.sqrt(p/9)*np.array( pauli_4D[3] ),
            np.sqrt(p/9)*np.array( pauli_4D[4] ),
            np.sqrt(p/9)*np.array( pauli_4D[5] ),
            np.sqrt(p/9)*np.array( pauli_4D[6] ),
            np.sqrt(p/9)*np.array( pauli_4D[7] ),
            np.sqrt(p/9)*np.array( pauli_4D[8] )
        ]
    else:
        print("depolarizing_channel:: dimension", dim," not yet support please try with dim = 2 or dim = 3")
        exit()

    return depo

def coherent_state(alpha):
    """
    coherent state if alpha is small and 2 order and above can be neglected
    """
    return np.exp(-np.abs(alpha)**2/2)*( zero + alpha * one )

def direct_sum(A, B):
    """
    computes the direct sum of two matrices A and B
    """
    return block_diag(A, B)

class QKD:
    def __init__(self, dim_a, dim_b, n_of_signal_states, list_states_a, list_of_prob_a, list_states_b, list_of_prob_b,
    povm_a=None,
    povm_b=None,
    kraus_a=None,
    kraus_b=None,
    sifting_phase=None,
    key_map=None,
    pinching_channel=None,
    key_map_povm=None,
    announcement_register_a=[zero, one], announcement_register_b=[zero, one]):

        # dimensions
        self.da     = dim_a
        self.db     = dim_b
        self.dtot   = dim_a*dim_b
        self.nst    = n_of_signal_states

        # states
        self.list_states_a  = list_states_a
        self.list_of_prob_a = list_of_prob_a
        self.list_states_b  = list_states_b
        self.list_of_prob_b = list_of_prob_b
        self.rho_aa         = []
        self.__get_state_aa()

        # rho_AB
        self.rho_ab = []

        # -- communication --
        # povm
        self.povm_a = povm_a
        self.povm_b = povm_b
        self.povm   = []
        if povm_a is None:
            self.__get_povm_a()
        if povm_b is None:
            self.__get_povm_b()
        self.__get_povm()

        # public string announcement
        self.kraus_a        = kraus_a
        self.register_a     = announcement_register_a # list of possbile bits in a_tilde and a_bar
        self.dim_register_a = len(announcement_register_a)

        self.kraus_b        = kraus_b
        self.register_b     = announcement_register_b # list of possbile bits in a_tilde and a_bar
        self.dim_register_b = len(announcement_register_b)
        
        self.public_string_announcement = [] # kraus operators
        if kraus_a is None:
            self.__get_kraus_a()
        if kraus_b is None:
            self.__get_kraus_b()
        self.__get_kraus()

        # sifting phase postselection
        self.sifting_phase_postselection = sifting_phase
        if sifting_phase is None:
            self.__get_sifting_phase_projector()
        self.ppass = 0

        # key map (isometry by adding a new register R)
        self.key_map = key_map
        if key_map is None:
            self.__get_key_map()

        # pinching channel to dechoere the register into a classical register
        self.pinching_channel = pinching_channel
        if pinching_channel is None:
            self.__get_pinching_channel()

        # key map povm
        if key_map_povm is not None:
            self.set_key_map_povm(key_map_povm)

        # -- constraints estimation --
        self.rho_0 = []
        # complete set of orthonormal observables in A
        self.orth_set_a = []
        self.set_operator_basis_a ( get_complete_matrix_basis(self.da) )

        # complete set of orthonormal observables in AB
        self.orth_set_ab = []
        self.set_operator_basis_ab( get_complete_matrix_basis(self.dtot) )

        # tilde{Gamma} and Omega
        self.Gamma_t = []
        self.Omega   = []
        self.hermitian_operator_basis = []
        self.__get_full_hermitian_operator_basis()

        # constraints
        self.gamma = []
        self.Gamma = []

        # -- algorithm parameters --
        self.epsilon     = 0
        self.maxit       = 0
        self.finesse     = 0
        self.solver_name = 0

        # algorithm estimations
        self.primal_sol   = 0
        self.dual_sol     = 0
        self.grad_f       = []
        self.delta_rho_i  = []

# private

    def __get_state_aa(self):
        """
        private: compute rho_aa'
        """
        # compute psi_aa'
        psi_aa = 0
        for ii in range(self.nst):
            psi_aa = psi_aa + np.sqrt(self.list_of_prob_a[ii])*np.kron( self.list_states_a[ii], self.list_states_b[ii] )
        psi_aa = psi_aa / np.sqrt(np.dot(np.conj(psi_aa), psi_aa)) # normalize

        # define rho_aa and check if it is physical
        self.rho_aa = np.outer( psi_aa, np.conj(psi_aa) )

        # check if it is physical
        if( np.abs( np.trace(self.rho_aa) - 1.) >= 1e-8): print("Tr[rho_aa] != 1 (", np.trace(self.rho_aa), ")")
        if( np.allclose( np.conj(self.rho_aa).T, self.rho_aa)  == False ): print("rho_aa NOT hermitian")
        if( np.all( np.linalg.eigvals(self.rho_aa) < - 1e-8)): print("rho_aa is NEGATIVE")

    def __get_povm_a(self):
        """
        private: automatically generates Alice POVM
        """
        # Alice
        self.povm_a = [ 2*j*np.outer( i, np.conj(i) ) for i, j in zip(self.list_states_a, self.list_of_prob_a)]
        is_povm(self.povm_a, "POVMA private::")

    def __get_povm_b(self):
        """
        private: automatically generates Bob POVM
        """
        # Bob
        self.povm_b = [ 2*j*np.outer( i, np.conj(i) ) for i, j in zip(self.list_states_b, self.list_of_prob_b)]
        is_povm(self.povm_b, "POVMB private::")

    def __get_povm(self):
        """
        private: automatically generates AB POVM
        """
        # AB: the POVM of the entire system is given by
        POVM = []
        for ii in self.povm_a:
                for jj in self.povm_b:
                    POVM.append( np.kron( ii, jj ) )
        is_povm(POVM, "POVM private::")
        self.povm = POVM

    def __get_kraus_a(self):
        """
         private: PUBLIC STRING ANNOUNCEMENT: kraus operators A
        """
        #   kraus operators of A dim = da*4-by-4
        self.kraus_a = []
        for ii, a_tilde in enumerate(self.register_a):
            temp = 0
            for jj, a_bar in enumerate(self.register_a):
                temp += np.kron( sqrtm(self.povm_a[jj+ii*len(self.register_a)]), np.kron(a_tilde[:, np.newaxis], a_bar[:, np.newaxis]))
            self.kraus_a.append(temp)
        is_kraus(self.kraus_a, "self.kraus_a private::")

    def __get_kraus_b(self):
        """
         private: PUBLIC STRING ANNOUNCEMENT: kraus operators B
        """
        #   kraus operators of B dim = db*4-by-4
        self.kraus_b = []
        for ii, b_tilde in enumerate(self.register_b):
            temp = 0
            for jj, b_bar in enumerate(self.register_b):
                temp += np.kron( sqrtm(self.povm_b[jj+ii*len(self.register_b)]), np.kron(b_tilde[:, np.newaxis], b_bar[:, np.newaxis]))
            self.kraus_b.append(temp)
        is_kraus(self.kraus_b, "self.kraus_b private::")

    def __get_kraus(self):
        """
         private: PUBLIC STRING ANNOUNCEMENT: kraus operators AB
        """ 
        #   The total Kraus representation of the Public Announcement is
        K = []
        for ii in self.kraus_a:
                for jj in self.kraus_b:
                    K.append( np.kron(ii, jj))
        is_kraus(K, "K private::")
        self.public_string_announcement = K

    def __get_sifting_phase_projector(self):
        """
         private: SIFTING PHASE: acts like a projector with dimension [da*4*db*4]-by-[da*4*db*4]
        """ 
        self.sifting_phase_postselection = np.kron( np.eye(self.da), np.kron( k0b0, np.kron( np.kron(np.eye(2), np.eye(self.db)), np.kron( k0b0, np.eye(2) ))) ) +\
            np.kron( np.eye(self.da), np.kron( k1b1, np.kron( np.kron(np.eye(2), np.eye(self.db)), np.kron( k1b1, np.eye(2) ))) )

    def __get_key_map(self):
        """
         private: KEY MAP: is a isometry which creates a new register R which stores the information on the bits and it is a [2*da*4*db*4]-by-[2*da*4*db*4] matrix
        """ 
        self.key_map = np.kron( zero[:, np.newaxis], np.kron( np.kron(np.eye(self.da), np.eye(2)), np.kron( k0b0, np.kron(np.eye(self.db), np.eye(4))) )) +\
                       np.kron( one[:, np.newaxis] , np.kron( np.kron(np.eye(self.da), np.eye(2)), np.kron( k1b1, np.kron(np.eye(self.db), np.eye(4))) ))

    def __get_pinching_channel(self):
        """
         private: PINCHING CHANNEL: decohere the register R. It has the effect of making R a classical register
        """ 
        self.pinching_channel = [
            np.kron( k0b0 , np.kron( np.kron(np.eye(self.da), np.kron( np.eye(4), np.eye(self.db) )), np.eye(4) )),
            np.kron( k1b1 , np.kron( np.kron(np.eye(self.da), np.kron( np.eye(4), np.eye(self.db) )), np.eye(4) ))
        ]

    def __get_full_hermitian_operator_basis(self):
        """
        computes the complete set of hermitian operators for AB
        -----------------------------------------------------------------
        fills:
            self.Gamma_t
            self.Omega
            self.hermitian_operator_basis
        """
        # construct the basis with POVM, Theta and Omega
        orth_set = []
        orth_set.append(self.povm[0])
        for ii, elm in enumerate(self.povm):
            flg = 0#set flag
            for jj, val in enumerate(orth_set):
                trn = np.trace(np.conj(elm).T @ val)
                if(abs(trn)>=1e-8):#if it is NOT null elm and val are NOT orthogonal
                    flg = 1
            if(flg == 0):#if elm is orthogonal to all
                orth_set=np.append(orth_set, [elm], axis=0)
        # Gamma = (POVM, Theta)
        orth_set = extend_basism(self.orth_set_a, orth_set)
    
        # Gram-Schmidt decomoposition
        POVM_tilde = gram_schmidtm(orth_set)
        # check if satisfy POVM_tilde properties
        # is_povm(POVM_tilde)

        # Gamma_tilde = {Gamma_tilde_k} U {Omega_j}
        kk = np.shape(POVM_tilde)[0]
        self.Gamma_t = POVM_tilde
        Gamma_tilde = extend_basism(self.orth_set_ab, POVM_tilde)
        jj = np.shape(Gamma_tilde)[0] - kk
        self.Omega = Gamma_tilde[kk:]
        # check that the elements of the basis are orthogonal
        for ii, elm in enumerate(Gamma_tilde):
            for jj, val in enumerate(Gamma_tilde):
                if (ii!=jj):
                    if(abs(np.trace(np.conj(elm).T @ val ))>=1e-11):
                        print("__get_full_hermitian_operator_basis:: element", ii, "and", jj, "are NOT orthogonal. ",abs(np.trace(np.conj(elm).T @ val )))
        #check if the basis contains 16 elements
        if(np.shape(Gamma_tilde)[0] != np.shape(Gamma_tilde)[1]*np.shape(Gamma_tilde)[2]):
            print("__get_full_hermitian_operator_basis:: the basis contains", np.shape(Gamma_tilde)[0], "!=", np.shape(Gamma_tilde)[1]*np.shape(Gamma_tilde)[2])
        self.hermitian_operator_basis = Gamma_tilde

    def __get_rho_0(self):
        """
        compute rho_0 as first step of step 1
        """
        # needs operators to find the orthonormal basis
        self.__get_full_hermitian_operator_basis()
        
        # basis = {Gamma_tilde} U {Omega}
        self.rho_0 = 0
        for ii in self.hermitian_operator_basis: 
            self.rho_0 = self.rho_0 + np.trace( ii @ self.rho_ab ) * ii
        self.rho_0 =  self.rho_0 /np.trace( self.rho_0)
        is_physical(self.rho_0, name="self.rho_0")

# public

    def set_key_map_povm(self, key_map_povm):
        """
        useful to set Z  s.t.  D(rho_AB|| Z @ rho_AB @ Z)  
        """
        self.public_string_announcement = [np.eye(self.dtot)]
        self.sifting_phase_postselection = np.eye(self.dtot)
        self.key_map = np.eye(self.dtot)
        self.pinching_channel = key_map_povm

    def set_povm_a(self, povma):
        """
        sets Alice POVM
        """
        # Alice
        self.povm_a = povma
        is_povm(self.povm_a, "set_povm_a::")

        # re-calculate povm
        self.__get_povm()

    def set_povm_b(self, povmb):
        """
        sets Bob POVM
        """
        # Bob
        self.povm_b =povmb
        is_povm(self.povm_b, "set_povm_b::")

        # re-calculate povm
        self.__get_povm()

    def set_povm(self, povm_ab):
        """
        set POVM AB
        -----------------------------------------------------------------
        fills:
            self.povm
        """
        self.povm = povm_ab
        is_povm(povm_ab, "set_povm::")

    def set_kraus_a(self, kraus_a):
        """
         sets kraus operators A
        """
        self.kraus_a = kraus_a
        is_kraus(self.kraus_a, "self.kraus_a::")

        # re-calculate kraus
        self.__get_kraus()

    def set_kraus_b(self, kraus_b):
        """
         sets kraus operators B
        """
        #   kraus operators of B dim = db*4-by-4
        self.kraus_b =  kraus_b
        is_kraus(self.kraus_b, "self.kraus_b::")

        # re-calculate kraus
        self.__get_kraus()

    def set_public_string_announcement(self, public_string_announcement):
        """
        set Kraus operators for the public string announcement
        -----------------------------------------------------------------
        fills:
            self.public_string_announcement
        """
        is_kraus(public_string_announcement, "public_string_announcement::")
        self.public_string_announcement = public_string_announcement # kraus operators

    def set_sifting_phase_postselection(self, sifting):
        """
        set the projector of the sifting phase
        -----------------------------------------------------------------
        fills:
            self.sifting_phase_postselection
        """
        self.sifting_phase_postselection = sifting

    def set_key_map(self, key_map):
        """
        set the list of operators that will try out the key
        -----------------------------------------------------------------
        fills:
            self.key_map
        """
        self.key_map = key_map

    def set_pinching_channel(self, pinching_channel):
        """
        set the list of operators that decohere tha register R
        -----------------------------------------------------------------
        fills:
            self.pinching_channel
        """
        self.pinching_channel = pinching_channel

    def set_operator_basis_a(self, orth_set_of_observables_a):
        """
        (use only for P&M schemes) operator basis in A
        input must have dimension dim A
        -----------------------------------------------------------------
        fills:
            self.orth_set_a
        """
        orth_set_of_observables_a_in_ab = []

        # enlarge to the entire AB system
        for ii in orth_set_of_observables_a:
            orth_set_of_observables_a_in_ab.append( np.kron(ii, np.eye(self.db)) )

        self.orth_set_a = orth_set_of_observables_a_in_ab

    def set_operator_basis_ab(self, orth_set_of_observables_ab):
        """
        set operator basis in AB
        -----------------------------------------------------------------
        fills:
            self.orth_set_ab
        """
        self.orth_set_ab = orth_set_of_observables_ab

    def apply_quantum_channel(self, channel_kraus_operator):
        """
        public: transforms system A' into B providing the Kraus representation (channel_kraus_operator) of the quantum channel
        """
        rho_ab, sumkt = 0., 0.

        # apply operators
        for jj in channel_kraus_operator:
            kraus_temp = np.kron( np.eye(self.da), jj)
            rho_ab = rho_ab + kraus_temp @ self.rho_aa @ np.conj( kraus_temp ).T
            sumkt = sumkt + np.conj( kraus_temp ).T @ kraus_temp
        
        # normalize
        # rho_ab = rho_ab / np.trace(rho_ab)

        # check if it is physical
        is_physical(rho_ab, 1e-8, "rho_ab")
    
        # check the depolarization satisfy Kraus representation property
        if( np.allclose( sumkt, np.eye(self.dtot)) == False): print("Depo4 is not kraus.")
        
        # add variable
        self.rho_ab = rho_ab

        # print purity
        print("purity of rho_i     =", np.real(np.trace( rho_ab @ rho_ab )))

    def set_constraints(self, gamma, Gamma):
        """
        set contraints where Tr[rho_ab @ Gamma_j] = gamma_j
        -----------------------------------------------------------------
        fills:
            self.gamma
            self.Gamma
        """
        self.gamma = gamma
        self.Gamma = Gamma
        self.__get_rho_0()

    def G_map(self, rho):
        '''
        G map defined as (isometry @ sifting @ Kraus) @rho@ (Kraus @ sifting @ isometry)**+
        '''
        rho_temp = 0.
        for ii in self.public_string_announcement:
            rho_temp = rho_temp + ii @ rho @ np.conj(ii).T
        rho_temp = self.sifting_phase_postselection @ rho_temp @ self.sifting_phase_postselection
        Prob_pass = np.real(np.trace( rho_temp ))
        rho_temp = rho_temp / Prob_pass
        rho_temp = self.key_map @ rho_temp @ np.conj(self.key_map).T
        # rho_temp = rho_temp / np.real(np.trace(rho_temp))
        return rho_temp, Prob_pass

    def Z_map(self, rho):
        """
        apply pinching channel
        """
        return sum( [ ii @ rho @ ii for ii in self.pinching_channel ] )

    def G_inverse_map(self, rho):
        '''the inverse G map'''
        rho_temp = 0.
        rho_temp = np.conj(self.key_map).T @ rho @ self.key_map
        rho_temp = self.sifting_phase_postselection @ rho_temp @ self.sifting_phase_postselection
        rho_fin = 0.
        for ii in self.public_string_announcement:
            rho_fin = rho_fin + np.conj( ii ).T @ rho_temp @ ii
        #rho_fin = rho_fin / np.real(np.trace(rho_fin))
        return rho_fin

    def set_primal_problem(self):
        '''
        sets the semidefinite program (SDP)
        -----------------------------------------------------------------
        minimize   : tr[ X**T @ grad_f ]    
        subject to :  X + rho >> 0                            (positivity)
                        Tr[X + rho] == 1                        (unitary trace)
                        X + rho == (X + rho)**+                 (hermiticity)
                        sum_j Tr[(X + rho) Gamma_j] = gamma_j   (measurments)
        where X is the solution of the SDP.

        return X
        '''
        # dimension of the matrix
        d = self.dtot

        # minimize: X.T @ Grad_f
        X = cp.Variable((d, d), complex=True)

        # subject to:
        constraints = [] # The operator >> denotes matrix inequality.
        constraints = constraints + [ X + self.rho_i >> 0 ] # Pos(H_AB)
        constraints = constraints + [ cp.real(cp.trace( X )) == 0 ] # trace sum up to 1
        constraints = constraints + [ X == cp.conj(X).T ] # is a Hermitian
        for ii, elm in enumerate(self.Gamma):
            constraints = constraints + [ cp.trace( (X + self.rho_i) @ elm)  == self.gamma[ii]]

        # construct problem
        obj = cp.Minimize( cp.real(cp.trace( X.T @  self.grad_f )) )
        primal_problem = cp.Problem( obj, constraints )
        return primal_problem, X

    def compute_primal(self, epsilon = 1e-10, maxit = 20, finesse = 10, solver_name = 'MOSEK', solver_verbosity = False):
        # epsilon=1e-10, maxit=20, finesse=10, solver_name='CVXOPT'
        ''''
        computes the primal problem with compute_primal_problem
        -----------------------------------------------------------------
        '''
        # set
        counter = 1
        self.rho_i = self.rho_ab

        # start algorithm
        while(counter <= maxit):

            # define the two states
            rho_4, self.ppass = self.G_map(self.rho_i)
            rho_5 = self.Z_map(rho_4)

            # primal_sol and gradient
            # compute f(rho) = D(rho_4||rho_5)
            self.primal_sol = np.real(relative_entropy( rho_4, rho_5 )) 

            # gradf(rho)**T = G**+ ( log[G(rho_i)] ) - G**+ ( log[Z(G(rho_i))] ) 
            self.grad_f = (self.G_inverse_map( logm( rho_4 )/np.sqrt(2) ) - self.G_inverse_map( logm( rho_5 )/np.sqrt(2) ))

            # print iteration result
            print("iteration", counter, " f(rho) =", self.primal_sol)

            # solve problem
            prob, X = self.set_primal_problem()
            
            try:
                prob.solve(solver=solver_name, verbose=solver_verbosity)
            except:
                print(solver_name," failed")
                if solver_name == 'MOSEK':
                    print("trying with CVXOPT")
                    prob.solve(solver='CVXOPT', verbose=solver_verbosity)
                elif solver_name == "CVXOPT":
                    print("trying with MOSEK")
                    prob.solve(solver='MOSEK', verbose=solver_verbosity)
                
            # check if there is any problem in the minimization procedure
            if prob.status in ["infeasible", "unbounded"]:# Otherwise, prob.value is inf or -inf, respectively.
                print("Status problem: %s" % prob.status)
                return 0
            else:
                # solution
                self.delta_rho_i = X.value

            # check if the problem is converged
            if(abs( np.trace( self.delta_rho_i.T @ self.grad_f)) <= epsilon):
                print("Algorithm exited at:", counter, "step")
                break

            # find l\in (0,1) s.t. min f(rho_i + l* delta_rho)
            tt = 0.
            f_1 = self.primal_sol
            for ii in np.linspace(0., 1., finesse):
                rho_4, _ = self.G_map(self.rho_i + self.delta_rho_i * ii)
                rho_5 = self.Z_map(rho_4)
                f_2 = np.real(relative_entropy(rho_4, rho_5))
                if (f_2 <= f_1):
                    tt = ii
                    f_1 = f_2

            # if f_1 == primal_sol the next step will be the same
            if(abs(f_1 - self.primal_sol) <= 1e-10): break

            # assign rho_{i+1} for the next iteration
            self.rho_i = self.rho_i + self.delta_rho_i * tt

            # update counter
            counter = counter + 1

            if (counter == maxit):
                print("algorithm reached maxit =", maxit)
                break

        #return
        return self.primal_sol

    def set_dual_problem(self):
        Z = cp.Variable(len(self.gamma))
        obj_exact = cp.sum(Z @ self.gamma)
        sum_Z = np.sum(ii * G for ii, G in zip(Z, self.Gamma))

        Y = cp.Variable()
        dual_obj = cp.Maximize(cp.real(Y + obj_exact))
        dual_constr = [( sum_Z + Y * np.eye(self.dtot) << self.grad_f)]
        dual_problem = cp.Problem(dual_obj, dual_constr)
        return dual_problem

    def compute_dual(self, solver_name = 'MOSEK', solver_verbosity = False):
        '''
        computes the dual problem of the SDP
        -----------------------------------------------------------------
        maximize   : Y @ gamma = (Y_1, ..., Y_n) @ (gamma_1, ..., gamma_n)
        subject to : sum_i Y_i*Gamma_i << grad_f
        where Y is the solution of the SDP.
        
        return Y = (Y_1, ..., Y_n)
        '''
        # set problem
        dual_problem = self.set_dual_problem()

        # solve
        try:
            dual_problem.solve(solver=solver_name, verbose=solver_verbosity)
        except:
            print(solver_name," failed, try CVXOPT")
            dual_problem.solve(solver='CVXOPT', verbose=solver_verbosity)


        # check
        if dual_problem.status in ["infeasible", "unbounded"]:
            print("Status problem: %s" % dual_problem.status)

        self.dual_sol =  np.real(self.primal_sol - np.trace( self.rho_i @ self.grad_f ) + dual_problem.value)

        # print difference between setp 1 and step 2
        print("precision (step1-step2) =", abs(self.dual_sol - self.primal_sol))

        # return
        return self.dual_sol