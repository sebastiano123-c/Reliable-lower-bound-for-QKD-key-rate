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

! DECLARATIONS
! dimensions
    INTEGER                                     :: da, db, dtot, nst                        
    parameter                                   (da=4, db=2, dtot=da*db, nst=4)             ! nst = # signal states
    complex(8)                                  :: basis(nst,da)
    complex(8)                                  :: signal(nst,db)
! probabilities
    real                                        :: pz, px, P_A(nst), P_B(nst)
    parameter                                   (pz=0.5, px=1-pz, p_a=[pz,pz,px,px]/2., p_b=[pz,pz,px,px]/2.)
! identities
    complex(8)                                  :: id_a(da,da), id_b(db,db), id_t(dtot,dtot)! id_A, id_B
! single qubits Z basis and X basis
    complex(8), DIMENSION(2,2)                  :: Z, X
    parameter                                   (Z=reshape([1,0,0,1],[2,2]))
    parameter                                   (X=reshape([1,1,1,-1]/sqrt(2.),[2,2]))
! classical register projectors
    complex(8), DIMENSION(2,2)                  :: sigma_00, sigma_11
    parameter                                   (sigma_00=reshape([1,0,0,0],[2,2]))
    parameter                                   (sigma_11=reshape([0,0,0,1],[2,2]))
! pauli matrices basis
    complex(8)                                  :: pauli(4,2,2)
    parameter                                   (pauli=reshape([ cmplx(1.,0.),cmplx(0.,0.),cmplx(0.,0.),cmplx(1.,0.),&
                                                                &cmplx(0.,0.),cmplx(1.,0.),cmplx(0.,1.),cmplx(0.,0.),&
                                                                &cmplx(0.,0.),cmplx(1.,0.),cmplx(0.,-1.),cmplx(0.,0.),&
                                                                &cmplx(1.,0.),cmplx(0.,0.),cmplx(0.,0.),-cmplx(1.,0.)],[4,2,2]))
! states
    complex(8)                                  :: psi_aa(dtot), rho_ab(dtot,dtot)          ! states
! dimension of povms
    INTEGER                                     :: npovma, npovmb, npovm                    ! number of povms
    parameter                                   (npovma=nst,npovmb=nst,npovm=npovma*npovmb)
    complex(8)                                  :: axa(npovma,da,da),bxb(npovmb,db,db)      ! projectors
    complex(8)                                  :: POVM_A(npovma,da,da),POVM_B(npovmb,db,db),POVM(npovm,dtot,dtot)!povms
! kraus operators
    integer, PARAMETER                          :: nka = 2, nkb = 2, nk = nka*nkb
    complex(8)                                  :: K_A(nka,4*da,da), K_B(nkb,4*db,db), Kraus(nk,16*dtot,dtot)
! sifting phase
    integer, PARAMETER                          :: nsi=16*dtot                              ! sifting dimension
    complex(8)                                  :: sifting(nsi,nsi)
! isometry
    complex(8)                                  :: isometry(2*nsi,nsi)
! key map
    INTEGER, PARAMETER                          :: nkm = 2, dkm = 2*nsi
    complex(8)                                  :: key_map(nkm, dkm, dkm)
! complete set of operators in A
    integer, PARAMETER                          :: ntheta = da**2
    complex(8)                                  :: Theta(ntheta, dtot, dtot)
! complet set of observables in AB
    complex(8)                                  :: Omega(da**2*db**2,dtot,dtot)
! hermitian operators for constraints in AB
    integer, PARAMETER                          :: ngamma = npovm+ntheta
    complex(8)                                  :: Gamma_i(ngamma,dtot,dtot)                ! measurments
    complex(8), allocatable                     :: Gamma_tilde(:,:,:)                       ! gram-schimdt decomposition
    complex(8), allocatable                     :: Omega_j(:,:,:), Gamma_tilde_k(:,:,:)
! quantum channel
    complex(8)                                  :: depo(db**2,dtot,dtot) ! depolarization
! variables of the protocol
    real, PARAMETER                             :: N_start = 0.1, N_stop = 0.13
    integer, PARAMETER                          :: N_steps = 10
    integer, PARAMETER                          :: Maxit = 10
    integer, PARAMETER                          :: finesse = 10
    real(8), PARAMETER                          :: epsilon = 1E-10
! constants
    real, parameter                             :: PI = 4.0*ATAN(1.0)
! others
    complex(8)                                  :: id(2,2), id_4(4,4), id_1_4(1,4,4)
    parameter                                   (id = reshape([1,0,0,1],[2,2]))
    parameter                                   (id_4 = reshape([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],[4,4]))
    complex(8), ALLOCATABLE                     :: tttemp(:,:,:)
    complex(8), ALLOCATABLE                     :: rho_temp(:,:), rho_0(:,:), grad_f(:,:)
    real(8)                                     :: f_rho, uu, sp
    real                                        :: start, finish, hp, th, trn
    INTEGER                                     :: ios, ii, jj, kk, ll, mm, sz
    CHARACTER(200)                              :: ps_file

    ! start computing
    call cpu_time(start)

!---------------------------------------------------------------------
!   PART 1): definition of the operators
!---------------------------------------------------------------------
! constants
    id_a = identity(int8(da))
    id_b = identity(int8(db))
    id_b = identity(int8(dtot))
    id_1_4(1,:,:)=id_4

! encoding qubit basis
    basis=identity(int8(nst))

! signal states
    signal(1,:)=Z(1,:);signal(2,:)=Z(2,:);signal(3,:)=X(1,:);signal(4,:)=X(2,:)

! define Theta (i.e. Eve cannot interact with A) 
!this defines the fine-grained constraints
!     Tr[Theta_j .o. id_b] = theta_j 
    do jj = 1, 4
        do kk = 1, 4
            Theta(kk+(jj-1)*4, :, :) = Kronecker_product(pauli(jj, :, :),&
            & Kronecker_product( pauli(kk, :, :), id_b ))
        enddo
    enddo

! measurments outcomes
!  1) Alice POVM 
    axa = cmplx(0.,0.)
    do ii = 1, npovma
        axa(ii,:,:) = Outer_product(basis(ii,:), conjg(basis(ii,:)))
    enddo
    POVM_A = axa
    ! check if it is a POVM
    do kk=1,da
        do jj=1,da
            call checkpoint( abs(POVM_A(1,kk,jj)+POVM_A(2,kk,jj)+POVM_A(3,kk,jj)+POVM_A(4,kk,jj)&
            &- id_a(kk,jj))<1e-4,text="Alice POVM is not a POVM",&
            &var=POVM_A(1,kk,jj)+POVM_A(2,kk,jj)+POVM_A(3,kk,jj)+POVM_A(4,kk,jj))
        enddo
    enddo
    ! check hermiticity and positivity
    do jj = 1, npovma
        call checkpoint(is_hermitian(POVM_A(jj,:,:),da), text="POVM_A not HERMITIAN", var=jj) !hermiticity
        call checkpoint(is_positive(POVM_A(jj,:,:),da), text="POVM_A not POSITIVE", var=jj) !positivity
    enddo

!    POVM Bob 
    bxb = cmplx(0.,0.)
    do ii = 1, npovma
        bxb(ii,:,:) = Outer_product(signal(ii,:), conjg(signal(ii,:)))
    enddo
    POVM_B = bxb/2
    ! check if they sum to identity
    do kk=1,db
        do jj=1,db
            call checkpoint( abs(POVM_B(1,kk,jj)+POVM_B(2,kk,jj)+POVM_B(3,kk,jj)+POVM_B(4,kk,jj)&
                                        & - id_b(kk,jj))<1e-3,text="Bob POVM is not a POVM",&
                                        &var=POVM_B(1,kk,jj)+POVM_B(2,kk,jj)+POVM_B(3,kk,jj)+POVM_B(4,kk,jj))
        enddo
    enddo
    ! check hermiticity and positivity
    do jj = 1, npovmb
        call checkpoint(is_positive(POVM_B(jj,:,:),db),text= "POVM_B not POSITIVE" ,var=jj) !positivity
        call checkpoint(is_hermitian(POVM_B(jj,:,:),db),text="POVM_B not HERMITIAN",var=jj) !hermiticity
    enddo

!    POVM AB 
    do kk = 1, npovmb
        do jj = 1, npovma
            POVM( jj + (kk -1)*npovma, :, :) = Kronecker_product( POVM_A(jj,:,:), POVM_B(kk,:,:) )
        enddo
    enddo
    ! check if they sum to identity
    if(all( (real(POVM(1,:,:)+POVM(2,:,:)+POVM(3,:,:)+POVM(4,:,:) - id_t)) >= 1e-4)) stop 'AB POVM is not a POVM'
    ! check hermiticity and positivity
    do jj = 1, npovm
        call checkpoint(is_positive(POVM(jj,:,:),dtot),text= "POVM not POSITIVE" ,var=jj) !positivity
        call checkpoint(is_hermitian(POVM(jj,:,:),dtot),text="POVM not HERMITIAN",var=jj) !hermiticity
    enddo

! 5* Kraus representation of announcements
!    Alice KRAUS:
    K_A(1,:,:) = Kronecker_product(sqrtM(POVM_A(1,:,:), da) ,&  ! A
                    & Kronecker_product( column(Z(1,:)),&       ! tilde{A}
                    & column(Z(1,:)))) &                        ! bar{A}
                &+&
                & Kronecker_product(sqrtM(POVM_A(2,:,:), da) ,& ! A
                    & Kronecker_product( column(Z(1,:)), &      ! tilde{A}
                    & column(Z(2,:))))                          ! bar{A}
    K_A(2,:,:) = Kronecker_product(sqrtM(POVM_A(3,:,:), da) ,&  ! A
                    & Kronecker_product( column(Z(2,:)),&       ! tilde{A}
                    & column(Z(1,:)))) &                        ! bar{A}
                &+&
                & Kronecker_product(sqrtM(POVM_A(4,:,:), da) ,& ! A
                    & Kronecker_product( column(Z(2,:)), &      ! tilde{A}
                    & column(Z(2,:))))                          ! bar{A}

!    Bob KRAUS:    
    K_B(1,:,:) = Kronecker_product(sqrtM(POVM_B(1,:,:), db) ,&  ! A
                    & Kronecker_product( column(Z(1,:)),&       ! tilde{A}
                    & column(Z(1,:)))) &                        ! bar{A}
                &+&
                & Kronecker_product(sqrtM(POVM_B(2,:,:), db) ,& ! A
                    & Kronecker_product( column(Z(1,:)), &      ! tilde{A}
                    & column(Z(2,:))))                          ! bar{A}
    K_B(2,:,:) = Kronecker_product(sqrtM(POVM_B(3,:,:), db) ,&  ! A
                    & Kronecker_product( column(Z(2,:)),&       ! tilde{A}
                    & column(Z(1,:)))) &                        ! bar{A}
                &+&
                & Kronecker_product(sqrtM(POVM_B(4,:,:), db) ,& ! A
                    & Kronecker_product( column(Z(2,:)), &      ! tilde{A}
                    & column(Z(2,:))))                          ! bar{A}
!    TOTAL KRAUS Alice tensor product Bob gives the sifting operator in Kraus representation
    do kk = 1, nka
        do jj = 1, nkb
            Kraus(jj+(kk-1)*nka,:,:)=Kronecker_product(K_A(kk,:,:),K_B(jj,:,:))
        enddo
    enddo

! 6* calculate the projector of the sifting-phase
    sifting = Kronecker_product( id_a , &                                        ! A
            & Kronecker_product( sigma_00 , &   ! tilde{A}
            & Kronecker_product( id_b, &                                      ! bar{A}
            & Kronecker_product( id, &                                      ! B 
            & Kronecker_product( sigma_00, &      ! tilde{B}
            & id))))) &                                                     ! bar{B}
        & + &                                                               !+                        
            & Kronecker_product( id_a , &                                     ! A
            & Kronecker_product( sigma_11 , &     ! tilde{A}
            & Kronecker_product( id_b, &                                      ! bar{A}
            & Kronecker_product( id, &                                      ! B
            & Kronecker_product( sigma_11, &      ! tilde{B}
            & id)))))                                                       ! bar{B}

! 7* consider the isometry V such that Alice stores 0 (1) in her key when she obtains outcome P1 or P3 (P2 or P4)
    !  V = |0>_R .o. |0><0|_t{A} .o. |0><0|_b{A} .o. |0><0|_t{B} +
    !    + |1>_R .o. |0><0|_t{A} .o. |1><1|_b{A} .o. |0><0|_t{B}
    isometry = Kronecker_product( column(Z(1,:)) , &                            ! R
        & Kronecker_product( id_a, &                                      ! A
        & Kronecker_product( id , &   ! tilde{A}
        & Kronecker_product( sigma_00 , &   ! bar{A}
        & Kronecker_product( id_b, &                                      ! B
        & Kronecker_product( id, &    ! tilde{B}
        & id &                                                          ! bar{B}
    &))))))&
    & + &
        & Kronecker_product( column(Z(2,:)), &                          ! R
        & Kronecker_product( id_a, &                                      ! A
        & Kronecker_product( id , &   ! tilde{A}
        & Kronecker_product( sigma_11 , &   ! bar{A}
        & Kronecker_product( id_b, &                                      ! B
        & Kronecker_product( id, &    ! tilde{B}
        & id &                                                          ! bar{B}
    &))))))

! 8* key_map channel Z: decohere R in his basis, which turns R into a classical register denoted Z^R
    key_map = cmplx(0.,0.)
    key_map(1,:,:) =  Kronecker_product(sigma_00, identity(int8(nsi)))
    key_map(2,:,:) =  Kronecker_product(sigma_11, identity(int8(nsi)))

! 9* define Gamma for constraints
    Gamma_i(:npovm,:,:) = POVM
    Gamma_i(npovm + 1:,:,:) = Theta

! 10* define the set of Hermitian operators for tomography
    sz = 1
    allocate(Gamma_tilde(sz,dtot,dtot))
    Gamma_tilde(1,:,:) = POVM(1,:,:)
    do ii = 1, npovm
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
    call extend_basism(Theta, Gamma_tilde, ntheta, dtot)
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
            do kk = 1, size(pauli, 1)
                allocate(rho_temp(dtot,dtot))
                rho_temp = 0.5/sqrt(2.) * Kronecker_product( pauli(kk,:,:), Kronecker_product( pauli(ll,:,:), pauli(mm,:,:) ))
                rho_temp = rho_temp / abs(mat_trace(matmul(rho_temp,rho_temp)))
                !check if is a Tr[G_mu G_mu]==1
                if(abs(mat_trace(matmul(rho_temp,rho_temp))- 1.)>=1e-5)print*,"Tr[G_mu G_mu] != 1",&
                & real(mat_trace(matmul(rho_temp,rho_temp)))
                Omega(mm+(ll-1)*(size(pauli, 1)),:,:) = rho_temp
                DEALLOCATE(rho_temp)
            enddo
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
        psi_aa = psi_aa + P_A(jj)*row(kronecker_product( column(basis(jj,:)), column(signal(jj,:)) ))
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
    do ii = 1, N_steps
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
        depo(1,:,:) = sqrt(1-3./4.*uu)*kronecker_product(id_a, pauli(1,:,:))
        depo(2,:,:) =      sqrt(uu/4.)*kronecker_product(id_a, pauli(2,:,:))
        depo(3,:,:) =      sqrt(uu/4.)*kronecker_product(id_a, pauli(3,:,:))
        depo(4,:,:) =      sqrt(uu/4.)*kronecker_product(id_a, pauli(4,:,:))

! 4* action on the state rho
        sz = size(rho_AB, 1)
        allocate(rho_temp(sz, sz)); rho_temp = cmplx(0.,0.)

        ! apply the depolarization to the state rho_AB
        do jj = 1, size(depo, 1)
            rho_temp = rho_temp + matmul(depo(jj,:,:), matmul(rho_AB, conjg(transpose( depo(jj,:,:))) ) )
        enddo
        ! check if the state after the depolarization is still physical
        rho_AB = rho_temp / mat_trace(rho_temp)
        call checkpoint(real(mat_trace(rho_AB))-1 <= 1e-8, text="Tr(rho_AB)/=1", var=mat_trace(rho_AB))
        call checkpoint(is_hermitian(rho_AB,sz), text="rho_AB is not hermitian")
        call checkpoint(is_positive(rho_AB,sz), text="rho_AB is not positive")
        deallocate(rho_temp)
        print*, "purity after quantum channel", real(mat_trace(matmul(rho_ab,rho_ab)))

! 6* Calculate the constraints      
        ! create rho_0
        allocate(rho_0(dtot, dtot))
        rho_0 = cmplx(0., 0.)
        do jj = 1, size(Gamma_tilde_k, 1)
            rho_0 = rho_0 + mat_trace(matmul(Gamma_tilde_k(jj,:,:),rho_ab))*Gamma_tilde_k(jj,:,:)
        enddo
        do jj = 1, size(Omega_j, 1)
            rho_0 = rho_0 + mat_trace(matmul(Omega_j(jj,:,:),rho_ab))*Omega_j(jj,:,:)
        enddo
        rho_0 = rho_0 / mat_trace(rho_0) ! renormalization
        rho_0 = rho_ab
        ! check if rho_0 is physical
        sz = size(rho_0, 1)
        call checkpoint(real(mat_trace(rho_0))-1 <= 1e-8, text="Tr(rho_0)/=1",var=mat_trace(rho_0))
        call checkpoint(is_hermitian(rho_0,sz), text="rho_0 is not hermitian", mat=rho_0)
        call checkpoint(is_positive(rho_0,sz), text="rho_0 is not positive")

! ALGORITHM 
! 7* STEP 1
        ALLOCATE(grad_f(size(rho_0, 1), size(rho_0,2)))
        ! call compute_primal(rho_0, f_rho, grad_f, id_1_4, id_4, id_4, key_map_test, Omega_j, Maxit, finesse, epsilon)
        call compute_primal(rho_0, f_rho, grad_f, Kraus, sifting, isometry, key_map, Omega_j, Maxit, finesse, epsilon)

! 8* STEP 2
        sp = f_rho - hp
        if(trn<1E-10)trn=0.
        write(10,'(es20.10," ",es20.10," ",es20.10," ",es20.10)')uu, sp, th, f_rho
        deallocate(rho_0, grad_f, stat=ios)
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
    DEALLOCATE(Gamma_tilde, Omega_j, stat=ios)
    call checkpoint(ios==0, text="final deallocation failed")
    
! end computation
    call cpu_time(finish)
    WRITE(*, '(f10.3, A)') finish-start, " s"
    call execute_command_line(("gnuplot -p "//trim(ps_file)))

endprogram main