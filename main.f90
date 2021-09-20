!---------------------------------------------------------------------
! QKD: RELIABLE LOWER BOUND FOR P&M - BB84 PROTOCOL SIMULATION
!---------------------------------------------------------------------
    !
    !   AUTHOR:
    !    Sebastiano Cocchi
    !
    !   LENGUAGE:
    !    FORTRAN 90
    !
    !   DESCRIPTION:
    !   "FINDS A RELIABLE LOWER BOUND ESTIMATION OF THE SECRET
    !    KEY RATE FOR P&M BB84 PROTOCOL WITH TWO MUBs.
    !    TO RECOVER THE ENTANGLED-BASED (EB) SCHEMES SECURITY
    !    PROOF, SOURCE-REPLACEMENT SCHEME IS PERFORMED (A 
    !    BRIEF DESCRIPTION IS PRESENTED HERE)."
    !
    !   FURTHER READINGS CAN BE FOUND AT:
    !       https://doi.org/10.22331/q-2018-07-26-77
    !       https://doi.org/10.1103/PhysRevResearch.3.013274
    !
    !   SYMBOLS INDEX:
    !    1)  .x. = tensor product
    !    2)  .+. = direct sum
    !    3)  **+ = hermitian conjugate
    !
    !   COMPILE SEQUENCE:
    !    gfortran -c QKDsimulation.f90 -Ofast -L'C:/lapack-3.9.0/build/bin' -llapack -lblas -unroll=n -g
    !    gcc -c sdplrlib.c
    !    gfortran -o QKDsimulation QKDsimulation.o sdplrlib.o
!---------------------------------------------------------------------
! Entangled Based (EB) scheme:
!---------------------------------------------------------------------
    ! " An entangled state composed of two photons is created.
    !   One particle is given to Alice, one to Bob.
    !   Both performs a measurment. "
    !
    !              |Psi>          
    !   Alice <______|______>   Bob    
    !   
    ! Prepare and Measure (P&M) scheme:
    ! " Alice send to Bob a qubit."
    !
    ! P&M schemes can be seen as EB schemes using the
    ! so-called Source-Replacement scheme in which Alice
    ! can choose between {|0>,|1>} and {|+>,|->}.
    ! She create an entangled state psi_AA', where A is a
    ! register in H^4 and A' is the register in which is encoded the bit,
    ! so it is a H^2.
    !
    !   QUANTUM CHANNEL: 
    !    * DEPOLARIZING CHANNEL;
    !
!---------------------------------------------------------------------
! SOURCE-REPLACEMENT SCHEME for P&M BB84:
!---------------------------------------------------------------------
    !
    !   SCHEME:
    !
    !                                            ........Eve............
    !   ..............Alice............          : U_E=I_A.0.E_A'(rho) :         .......Bob.......
    !   :               |PHI>_AA'     :          :........../\.........:         :    A' --> B   :
    !   :   ________________/\________:____________________/  \__________________:_______________
    !   :   A:{0,1,+,-}H^2      A':{0,1,+,-}H^2                                  :B:{0,1,+,-}H^2 :
    !   :       dim=2           dim=2 :                                          :     dim=2     :
    !   :        |                    :                                          :       |       :
    !   :        V                    :                                          :       V       :
    !   :      POVM_A                 :                                          :    POVM_B     :
    !   :.............................:                                          :...............:
    !
    !   STATE A    : STATE A'(B') : BASIS CHOICE : BIT VALUE 
    !       |0>    :     |0>      :      Z       :     0
    !       |1>    :     |1>      :      Z       :     1
    !       |+>    :     |+>      :      X       :     0
    !       |->    :     |->      :      X       :     1
    !   
!---------------------------------------------------------------------
! The program is diveded in two parts:
!   1) explanation of the conceptual steps of the procedure and
!      the declaration of the usefule operators;
!   2) the algorithm procedure implementing the operators defined in 
!      the previous point and SDP minimization procedure.
! N.B.: the algorithm is a recursive iteration incrementing the depolarization
!       probability.

include "src/debugging.f90"
include "src/matrices.f90"
include "src/QKD.f90"

program main
    use debugging
    use matrices
    use QKD
    implicit none

! Declarations
    complex(8), dimension(4)                    :: psi_aa
    complex(8), dimension(2,2)                  :: Z, X, id, null_matrix, sigma_00,sigma_11,&
                                                &sigma_pp, sigma_mm
    complex(8), dimension(4,4)                  :: rho_ab, id_4
    complex(8), dimension(:,:,:), ALLOCATABLE   :: POVM_A, POVM_B, POVM, K_A, K_B, Kraus, &
                                                & Gamma_i, Gamma_tilde, key_map, Theta, &
                                                & tttemp, Omega_j, Gamma_tilde_k
    complex(8), dimension(:,:), ALLOCATABLE     :: rho_temp, sifting, &
                                                & isometry, rho_0, grad_f
    real(8), dimension(:), ALLOCATABLE          :: p_i, p_tilde, theta_j
    complex(8)                                  :: signal(4,2), Omega(16,4,4), pauli(4,2,2), depo(4,4,4)
    real(8)                                     :: f_rho, epsilon, uu, sp
    real                                        :: Pz, Pz_A(4), Pz_B(4), trn
    real                                        :: start, finish, N_start, N_stop, hp, th
    real                                        :: PI
    INTEGER                                     :: ios, ii, jj, kk, ll, mm, sz, Maxit,&
                                                & finesse, N_steps
    INTEGER(8)                                  :: KeyLength
    CHARACTER(200)                              :: ps_file

    external trace_fn

    ! start computing
    call cpu_time(start)

!---------------------------------------------------------------------
!   PART 1): definition of the operators
!---------------------------------------------------------------------
! variables of the protocol
    Pz = 0.5        ! Z-basis probability
    N_start = 0.; N_stop = 0.12; N_steps = 5    ! # of steps 
    Maxit = 10      ! maximum # of iterations
    epsilon = 1E-10 ! convergence tightness
    finesse = 10    ! finesse need for searching tt

! constants
    PI = 4.0*ATAN(1.0)
    KeyLength = 10
    null_matrix(1,:) = [cmplx(0.,0.), cmplx(0.,0.)]
    null_matrix(2,:) = [cmplx(0.,0.), cmplx(0.,0.)]
    id = identity(2_8)+cmplx(0.,0.); id_4 = identity(4_8)+cmplx(0.,0.)!; id_64 = identity(64_8)+cmplx(0.,0.)

! Pauli matrices
    pauli(1,:,:)=id;pauli(2,:,:)=sigma("x");pauli(3,:,:)=sigma("y");pauli(4,:,:)=sigma("z")

! define probabilities P_i
    Pz_A=[Pz/2,Pz/2,(1-Pz)/2,(1-Pz)/2]
    Pz_B=Pz_A

! encoding qubits
    Z(1,:)=[1,0];Z(2,:)=[0,1];X(1,:)=[1,1]/sqrt(2.);X(2,:)=[1,-1]/sqrt(2.)

! signal states
    signal(1,:)=Z(1,:);signal(2,:)=Z(2,:);signal(3,:)=(Z(1,:)+Z(2,:))/sqrt(2.);signal(4,:)=(Z(1,:)-Z(2,:))/sqrt(2.)

! states siftingectors
    sigma_00 = Outer_product(signal(1,:), conjg(signal(1,:)))
    sigma_11 = Outer_product(signal(2,:), conjg(signal(2,:)))
    sigma_pp = Outer_product(signal(3,:), conjg(signal(3,:)))
    sigma_mm = Outer_product(signal(4,:), conjg(signal(4,:)))

! measurments outcomes
!  1) Alice POVM 
    allocate(POVM_A(4,2,2))
    POVM_A(1,:,:) = 2*Pz_A(1)*Outer_product(signal(1,:), conjg(signal(1,:)))
    POVM_A(2,:,:) = 2*Pz_A(2)*Outer_product(signal(2,:), conjg(signal(2,:)))
    POVM_A(3,:,:) = 2*Pz_A(3)*Outer_product(signal(3,:), conjg(signal(3,:)))
    POVM_A(4,:,:) = 2*Pz_A(4)*Outer_product(signal(4,:), conjg(signal(4,:)))
    ! check if it is a POVM
    do kk=1,size(POVM_A,2)
        do jj=1,size(POVM_A,3)
            call checkpoint( abs(POVM_A(1,kk,jj)+POVM_A(2,kk,jj)+POVM_A(3,kk,jj)+POVM_A(4,kk,jj)&
                                        & - id(kk,jj))<1e-3,text="Alice POVM is not a POVM",&
                                        &var=POVM_A(1,kk,jj)+POVM_A(2,kk,jj)+POVM_A(3,kk,jj)+POVM_A(4,kk,jj))
        enddo
    enddo
    ! check hermiticity and positivity
    sz = size(POVM_A, 2)
    do jj = 1, size(POVM_A, 1)
        call checkpoint(is_hermitian(POVM_A(jj,:,:),sz), text="POVM_A not HERMITIAN", var=jj) !hermiticity
        call checkpoint(is_positive(POVM_A(jj,:,:),sz), text="POVM_A not POSITIVE", var=jj) !positivity
    enddo

!    POVM Bob 
    allocate(POVM_B(4, 2, 2))
    POVM_B = POVM_A
    ! check if they sum to identity
    do kk=1,size(POVM_B,2)
        do jj=1,size(POVM_B,3)
            call checkpoint( abs(POVM_B(1,kk,jj)+POVM_B(2,kk,jj)+POVM_B(3,kk,jj)+POVM_B(4,kk,jj)&
                                        & - id(kk,jj))<1e-3,text="Bob POVM is not a POVM",&
                                        &var=POVM_B(1,kk,jj)+POVM_B(2,kk,jj)+POVM_B(3,kk,jj)+POVM_B(4,kk,jj))
        enddo
    enddo
    ! check hermiticity and positivity
    sz = size(POVM_B,2)
    do jj = 1, size(POVM_B,1)
        call checkpoint(is_positive(POVM_B(jj,:,:),sz),text= "POVM_B not POSITIVE" ,var=jj) !positivity
        call checkpoint(is_hermitian(POVM_B(jj,:,:),sz),text="POVM_B not HERMITIAN",var=jj) !hermiticity
    enddo

!    POVM AB 
    allocate(POVM(size(POVM_A, 1)*size(POVM_B, 1), size(POVM_A, 2)*size(POVM_B, 2), size(POVM_A, 3)*size(POVM_B, 3)))
    do kk = 1, size(POVM_B, 1)
        do jj = 1, size(POVM_A, 1)
            POVM( jj + (kk -1)*size(POVM_A, 1), :, :) = Kronecker_product( POVM_A(jj,:,:), POVM_B(kk,:,:) )
        enddo
    enddo
    ! check if they sum to identity
    if(all( (real(POVM(1,:,:)+POVM(2,:,:)+POVM(3,:,:)+POVM(4,:,:) - id_4)) >= 1e-4)) stop 'AB POVM is not a POVM'
    ! check hermiticity and positivity
    sz = size(POVM,2)
    do jj = 1, size(POVM,1)
        call checkpoint(is_positive(POVM(jj,:,:),sz),text= "POVM not POSITIVE" ,var=jj) !positivity
        call checkpoint(is_hermitian(POVM(jj,:,:),sz),text="POVM not HERMITIAN",var=jj) !hermiticity
    enddo

! 4* Theta (i.e. Eve cannot interact with A) 
    !this defines the fine-grained constraints
    !     Tr[Theta_j .o. id_2] = theta_j 
    allocate(Theta(4, 4, 4)); allocate(theta_j(size(Theta, 1)))
    do jj = 1, 4
        Theta(jj, :, :) = Kronecker_product( pauli(jj, :, :), id )
    enddo
    
! 5* Kraus representation of announcements
!    Alice KRAUS:
    allocate(K_A(2,8,2))
    sz = size(POVM_A(1,:,:),2)  
    K_A(1,:,:) = Kronecker_product(sqrtM(POVM_A(1,:,:), sz) ,&  ! A
                    & Kronecker_product( column(Z(1,:)),&       ! tilde{A}
                    & column(Z(1,:)))) &                        ! bar{A}
                &+&
                & Kronecker_product(sqrtM(POVM_A(2,:,:), sz) ,& ! A
                    & Kronecker_product( column(Z(1,:)), &      ! tilde{A}
                    & column(Z(2,:))))                          ! bar{A}
    K_A(2,:,:) = Kronecker_product(sqrtM(POVM_A(3,:,:), sz) ,&  ! A
                    & Kronecker_product( column(Z(1,:)),&       ! tilde{A}
                    & column(Z(1,:)))) &                        ! bar{A}
                &+&
                & Kronecker_product(sqrtM(POVM_A(4,:,:), sz) ,& ! A
                    & Kronecker_product( column(Z(1,:)), &      ! tilde{A}
                    & column(Z(2,:))))                          ! bar{A}

!    Bob KRAUS:    
    allocate(K_B(2,8,2))
    sz = size(POVM_B(1,:,:),2)
    K_B(1,:,:) = Kronecker_product(Kronecker_product( sqrtM(POVM_B(1,:,:),sz) , column(Z(1,:))), column(Z(1,:)))+&
                &Kronecker_product(Kronecker_product( sqrtM(POVM_B(2,:,:),sz) , column(Z(1,:))), column(Z(2,:)))
    K_B(2,:,:) = Kronecker_product(Kronecker_product( sqrtM(POVM_B(3,:,:),sz) , column(Z(2,:))), column(Z(1,:)))+&
                &Kronecker_product(Kronecker_product( sqrtM(POVM_B(4,:,:),sz) , column(Z(2,:))), column(Z(2,:)))
!    TOTAL KRAUS Alice tensor product Bob gives the sifting operator in Kraus representation
    allocate(Kraus( size(K_A,1)*size(K_B,1), size(K_A,2)*size(K_B,2), size(K_A,3)*size(K_B,3) ))
    do kk = 1, size(K_B,1)
        do jj = 1, size(K_A,1)
            Kraus( jj + (kk -1)*size(K_A, 1), :, :) = Kronecker_product(K_A(jj,:,:), K_B(kk,:,:))
        enddo
    enddo

! 6* calculate the projector of the sifting-phase
    allocate(sifting(size(K_A,2)*size(K_B,2), size(K_A,2)*size(K_B,2)))
    sifting = Kronecker_product( id , &                                        ! A
            & Kronecker_product( sigma_00 , &   ! tilde{A}
            & Kronecker_product( id, &                                      ! bar{A}
            & Kronecker_product( id, &                                      ! B 
            & Kronecker_product( sigma_00, &      ! tilde{B}
            & id))))) &                                                     ! bar{B}
        & + &                                                               !+                        
            & Kronecker_product( id , &                                     ! A
            & Kronecker_product( sigma_11 , &     ! tilde{A}
            & Kronecker_product( id, &                                      ! bar{A}
            & Kronecker_product( id, &                                      ! B
            & Kronecker_product( sigma_11, &      ! tilde{B}
            & id)))))                                                       ! bar{B}

! 7* consider the isometry V such that Alice stores 0 (1) in her key when she obtains outcome P1 or P3 (P2 or P4)
    !  V = |0>_R .o. |0><0|_t{A} .o. |0><0|_b{A} .o. |0><0|_t{B} +
    !    + |1>_R .o. |0><0|_t{A} .o. |1><1|_b{A} .o. |0><0|_t{B}
    allocate(isometry(size(z,1)*size(K_A,2)*size(K_B,2), size(K_A,2)*size(K_B,2) )); isometry = cmplx(0.,0.)
    isometry = Kronecker_product( column(Z(1,:)) , &                            ! R
                & Kronecker_product( id, &                                      ! A
                & Kronecker_product( sigma_00 , &   ! tilde{A}
                & Kronecker_product( sigma_00 , &   ! bar{A}
                & Kronecker_product( id, &                                      ! B
                & Kronecker_product( sigma_00, &    ! tilde{B}
                & id &                                                          ! bar{B}
    &))))))&
            & + &
                & Kronecker_product( column(Z(2,:)), &                          ! R
                & Kronecker_product( id, &                                      ! A
                & Kronecker_product( sigma_00 , &   ! tilde{A}
                & Kronecker_product( sigma_11 , &   ! bar{A}
                & Kronecker_product( id, &                                      ! B
                & Kronecker_product( sigma_00, &    ! tilde{B}
                & id &                                                          ! bar{B}
    &))))))

! 8* key_map channel Z: decohere R in his basis, which turns R into a classical register denoted Z^R
    allocate(key_map(2, size(isometry,1), size(isometry,1))); key_map = cmplx(0.,0.)
    key_map(1,:,:) =  Kronecker_product( Outer_product(Z(1,:),Z(1,:)), identity(64_8))
    key_map(2,:,:) =  Kronecker_product( Outer_product(Z(2,:),Z(2,:)), identity(64_8))

! 9* define Gamma for constraints
    allocate(Gamma_i(size(POVM, 1)*size(Theta, 1), size(POVM, 2), size(POVM, 3) ))
    Gamma_i(:size(POVM, 1), :, :) = POVM
    Gamma_i(size(POVM, 1) + 1:, :, :) = Theta

! 10* define the set of Hermitian operators for tomography
    sz = 1
    allocate(Gamma_tilde(sz, size(POVM,2), size(POVM,3)))
    Gamma_tilde(1,:,:) = POVM(1,:,:)
    do ii = 1, size(POVM, 1)
        ios = 0!set flag
        do jj = 1, sz
            trn = real(mat_trace(matmul(transpose(conjg(POVM(ii,:,:))), Gamma_tilde(jj,:,:))),4)
            if(abs(trn)>=1e-8) then!#if it is NOT null elm and val are NOT orthogonal
                ios = 1
            endif
        enddo
        if(ios == 0)then!#if elm is orthogonal to all
            allocate(tttemp(sz, size(Gamma_tilde,2), size(Gamma_tilde,3)))
            tttemp = Gamma_tilde
            DEALLOCATE(Gamma_tilde)
            ALLOCATE(Gamma_tilde( sz+1, size(tttemp,2), size(tttemp,3) ))
            do kk = 1, sz          
                Gamma_tilde(kk,:,:) = tttemp(kk,:,:)
            enddo
            Gamma_tilde(sz+1,:,:) = POVM(ii,:,:)
            DEALLOCATE(tttemp)
            sz = sz + 1
        endif
    enddo
    ! Gamma = (POVM, Theta)
    call extend_basism(Theta, Gamma_tilde, size(Theta, 1), size(Theta, 2))
    ! Gram-Schmidt process =>obtain k orthonormal Gamma_tildes
    kk = size(Gamma_tilde, 1)
    call gram_schmidtm(Gamma_tilde, kk,  size(Gamma_tilde,2),  size(Gamma_tilde,3) )
    ! orthogonality check
    do ll = 1, kk
        do mm = 1, kk
            if(ll/=mm) then
                if(abs(mat_trace(matmul(Gamma_tilde(ll,:,:),Gamma_tilde(mm,:,:))))>=1e-8) stop "Gamma_tilde is NOT orthogonal."
            endif
        enddo
    enddo
    Gamma_tilde_k = Gamma_tilde

! 11* Omega, set of tomographycally observables to complete the basis
    do ll = 1, size(pauli, 1)
        do mm = 1, size(pauli, 1)
            allocate(rho_temp(size(pauli,2)*size(pauli,2), size(pauli,3)*size(pauli,3)))
            rho_temp = 0.5 * Kronecker_product( pauli(ll,:,:), pauli(mm,:,:) )
            !check for hermiticity
            if(all( abs(rho_temp - transpose(conjg(rho_temp))) <= 1e-8 ) .eqv. .false.) print*, "Gamma", ii, ll, "not hermitian."
            !check if is a Tr[G_mu G_mu]==1
            if( abs(mat_trace(matmul(rho_temp, rho_temp)) - 1.) >= 1e-8) print*, "Tr[G_mu G_mu] != 1", &
                                                                & real(mat_trace(matmul(rho_temp, rho_temp)))
            Omega(mm+(ll-1)*(size(pauli, 1)),:,:) = rho_temp
            DEALLOCATE(rho_temp)
        enddo
    enddo
    ! extend basis using j Omegas
    call extend_basism(Omega, Gamma_tilde, size(Omega, 1), size(Omega, 2))
    jj = size(Gamma_tilde, 1) - kk
    ! which are
    allocate(Omega_j(jj, size(Omega, 2), size(Omega, 3)))
    Omega_j = Gamma_tilde( kk+1 :, :, : )
    ! orthogonality check
    do ll = 1, size(Gamma_tilde,1)
        do mm = 1, size(Gamma_tilde,1)
            if(ll/=mm) then
                if(abs(mat_trace(matmul(Gamma_tilde(ll,:,:),Gamma_tilde(mm,:,:))))>=1e-8) stop "Gamma_tildes are NOT orthogonal."
            endif
        enddo
    enddo
    ! Gamma_tilde = Gamma_tilde_k + Omega_j

!---------------------------------------------------------------------
!   PART 2): algorithm
!---------------------------------------------------------------------
! gnuplot output file for graphics
    ps_file = "QKD.ps"
    open(unit=10, file=ps_file)
    write(10, '("$data<< EOD")')

! 1* construct the state between A and B
!    create Psi_AA'
    psi_aa = cmplx(0., 0.)
    do jj = 1, 4
        psi_aa = psi_aa + Pz_A(jj)*row(kronecker_product( column(signal(jj,:)), column(signal(jj,:)) ))
    enddo
!    creation of the density matrix rho_{AB} given by (eqn.1)
    rho_ab = Outer_product( psi_aa, conjg(psi_aa) )
    !!!! PAY ATTENTION THIS PASS NORMALIZE IN ANY CASE !!!!
    rho_AB = rho_AB / real(mat_trace(rho_AB))
    ! check if the density operator is physical
    sz = size(rho_AB,1)
    call checkpoint(real(mat_trace(rho_AB))-1 <= 1e-8, text="Tr(rho_AB)/=1",var=mat_trace(rho_AB))
    call checkpoint(is_hermitian(rho_AB,sz), text="rho_AB is not hermitian")
    call checkpoint(is_positive(rho_AB,sz), text="rho_AB is not positive")
    print*, "rho_ab purity", real(mat_trace(matmul( rho_ab, rho_ab )))

! 2* starting point of the algorithm
    do ii = 0, N_steps
        ! qber
        uu = dble(N_stop - N_start)*dble(ii)/dble(N_steps)
        write(*,*) ""
        write(*,'(A, f6.3, A)') " [ QBER = ", uu, " ]"
        ! binary entropy
        hp = real(binary_entropy(uu), 4)
        ! theorical value
        th = 1. - 2.* hp
        if (th <= 1e-10) th = 0.

! 3* Depolarization channel:
        depo(1,:,:) = sqrt(1-3./4.*uu)*kronecker_product( id, pauli(1,:,:) )
        depo(2,:,:) =      sqrt(uu/4.)*kronecker_product( id, pauli(2,:,:) )
        depo(3,:,:) =      sqrt(uu/4.)*kronecker_product( id, pauli(3,:,:) )
        depo(4,:,:) =      sqrt(uu/4.)*kronecker_product( id, pauli(4,:,:) )

! 4* action on the state rho
        sz = size(rho_AB, 1)
        allocate(rho_temp(sz, sz)); rho_temp = cmplx(0.,0.)

        ! apply the depolarization to the state rho_AB
        do jj = 1, size(depo, 1)
            rho_temp = rho_temp + matmul( depo(jj,:,:), matmul( rho_AB, conjg(transpose( depo(jj,:,:))) ) )
        enddo
        ! check if the state after the depolarization is still physical
        rho_AB = rho_temp / mat_trace(rho_temp)
        call checkpoint(real(mat_trace(rho_AB))-1 <= 1e-8, text="Tr(rho_AB)/=1", var=mat_trace(rho_AB))
        call checkpoint(is_hermitian(rho_AB,sz), text="rho_AB is not hermitian")
        call checkpoint(is_positive(rho_AB,sz), text="rho_AB is not positive")
        deallocate(rho_temp)
        print*, "purity after quantum channel", real(mat_trace(matmul(rho_ab,rho_ab)))
        
! 6* Calculate the constraints
        allocate(p_i(size(Gamma_i, 1)))
        allocate(p_tilde(size(Gamma_tilde_k, 1)))
        do jj = 1, size(p_i)
            p_i(jj) = real(mat_trace( matmul( Gamma_i(jj,:,:), rho_ab) ))
        enddo
        do jj = 1, size(p_tilde)
            p_tilde(jj) = real(mat_trace( matmul( Gamma_tilde_k(jj,:,:), rho_ab) ))
        enddo
        ! create rho_0
        allocate(rho_0(size(rho_ab, 1), size(rho_ab, 2)))
        rho_0 = cmplx(0., 0.)
        do jj = 1, size(Gamma_tilde_k, 1)
            rho_0 = rho_0 + p_tilde(jj) * Gamma_tilde_k(jj,:,:)
        enddo
        rho_0 = rho_0 / mat_trace(rho_0) ! renormalization
        !rho_0 = rho_ab
        ! check if rho_0 is physical
        sz = size(rho_0, 1)
        call checkpoint(real(mat_trace(rho_0))-1 <= 1e-8, text="Tr(rho_0)/=1",var=mat_trace(rho_0))
        call checkpoint(is_hermitian(rho_0,sz), text="rho_0 is not hermitian", mat=rho_0)
        call checkpoint(is_positive(rho_0,sz), text="rho_0 is not positive")

! ALGORITHM 
! 7* STEP 1
        ALLOCATE(grad_f(size(rho_0, 1), size(rho_0,2)))
        call compute_primal(rho_0, f_rho, grad_f, Kraus, sifting, isometry, key_map, Omega_j, Maxit, finesse, epsilon)

! 8* STEP 2
        sp = f_rho - hp
        if(trn<1E-10)trn=0.
        write(10,'(es20.10," ",es20.10," ",es20.10," ",es20.10)')uu, sp, th, f_rho
        deallocate(rho_0, grad_f, p_i, p_tilde, stat=ios)
        call checkpoint(ios==0, text="rho_0 deallocation failed")
    enddo
    write(10,'("EOD")')
    write(10,'(a)')"set terminal wxt 0 size 3000,1600"        
    ! write(10,'(a)')"set size 1,1"        
    write(10,'("set xlabel ''depolarization probability''")')
    write(10,'("set ylabel ''secret key rate''")')
    write(10,'("set title ''reliable lower bound for BB84 protocol''")')
    write(10,'("set grid")')
    write(10,'("set key outside")')    
    write(10,'("#set xr [:]")')    
    write(10,'(A)')"p $data u 1:2 w lp t 'step 1', $data u 1:3 w lp t 'theoric', $data u 1:4 w lp t 'f_{\rho}'"
    close(10)!close file
!   END PROCEDURE
!---------------------------------------------------------------------------------------------------------------------

! deallocation
    DEALLOCATE(POVM_A, POVM_B, POVM, Gamma_i, Gamma_tilde, K_A, K_B, Kraus, isometry,&
                & sifting, key_map, Omega_j, stat=ios)
    call checkpoint(ios==0, text="final deallocation failed")
    
! end computation
    call cpu_time(finish)
    WRITE(*, '(f10.3, A)') finish-start, " s"
    call execute_command_line(("gnuplot -p "//trim(ps_file)))

endprogram main