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
include "src/DUALSIMP_MOD.f90"

program QKDsimulation
    use debugging
    use matrices
    use QKD
    use DUALSIMP_MOD
    use minimization
    implicit none

! Declarations
    complex(8), dimension(4)                    :: psi_aa
    complex(8), dimension(2,2)                  :: Z, X, id, null_matrix, sigma_00, sigma_11,&
                                                & sigma_pp, sigma_mm
    complex(8), dimension(4,4)                  :: rho_ab, id_4
    complex(8), dimension(:,:,:), ALLOCATABLE   :: POVM_A, POVM_B, POVM, K_A, K_B, Kraus, &
                                                & Gamma_i, Gamma_tilde, pinching, Theta, &
                                                & tttemp, Omega_j, Gamma_tilde_k
    complex(8), dimension(2,4,4)                :: key_map
    complex(8), dimension(:,:), ALLOCATABLE     :: rho_temp, proj, rho_4, rho_5, &
                                                & isometry, rho_0, delta_rho, grad_f, logrho
    complex(8), dimension(:), ALLOCATABLE       :: p_i, p_tilde, theta_j
    complex(8)                                  :: signal(4,2), Omega(16,4,4), pauli(4,2,2), depo(4,4,4)
    real(8), ALLOCATABLE                        :: C(:), AT(:,:), X_sol(:), Y_sol(:), upper_bounds(:)
    real(8)                                     :: f_rho, f_1, f_2, epsilon, tt, uu
    real                                        :: Pz, Pz_A(4), Pz_B(4), Ppass, trn
    real                                        :: start, finish
    real                                        :: PI
    INTEGER                                     :: ios, ii, jj, kk, ll, mm, sz, ierr, N_steps, Maxit,&
                                                & counter, finesse
    INTEGER(8)                                  :: KeyLength
    INTEGER, ALLOCATABLE                        :: IBASIS(:)
    CHARACTER(200)                              :: ps_file

    external trace_fn

    ! start computing
    call cpu_time(start)

! VALUES TO BE SET:
    Pz = 0.5        ! Z-basis probability
    N_steps = 10    ! # of steps 
    Maxit = 5       ! maximum # of iterations
    epsilon = 1E-8  ! convergence tightness
    finesse = 10    ! finesse need for searching tt

! constants
    PI = 4.0*ATAN(1.0)
    KeyLength = 10
    null_matrix(1,:) = [cmplx(0.,0.), cmplx(0.,0.)]
    null_matrix(2,:) = [cmplx(0.,0.), cmplx(0.,0.)]
    id = identity(2_8)+cmplx(0.,0.); id_4 = identity(4_8)+cmplx(0.,0.)!; id_64 = identity(64_8)+cmplx(0.,0.)
    pauli(1,:,:) = id; pauli(2,:,:) = sigma("x"); pauli(3,:,:) = sigma("y"); pauli(4,:,:) = sigma("z")

! probabilities
    Pz_A = [Pz/2, Pz/2, (1-Pz)/2, (1-Pz)/2]
    Pz_B = Pz_A

! encoding qubits
    Z(1,:) = [1,0]; Z(2,:) = [0, 1]; X(1,:) = [1, 1]/sqrt(2.0); X(2,:) = [1, -1]/sqrt(2.0)
    signal(1,:) = [1,0]; signal(2,:) = [0, 1]; signal(3,:) = [1, 1]/sqrt(2.0); signal(4,:) = [1, -1]/sqrt(2.0) 

! alice measurments
    sigma_00 = Kronecker_product(signal(1,:), conjg(signal(1,:)))
    sigma_11 = Kronecker_product(signal(2,:), conjg(signal(2,:)))
    sigma_pp = Kronecker_product(signal(3,:), conjg(signal(3,:)))
    sigma_mm = Kronecker_product(signal(4,:), conjg(signal(4,:)))

! output file
    ps_file = "analysis/QKD.ps"
    open(unit=10, file=ps_file)
    write(10, '("$data<< EOD")')

!---------------------------------------------------------------------
!   PART 1): definition of the operators
!---------------------------------------------------------------------
    key_map = 

!---------------------------------------------------------------------
!   PART 2): algorithm
!---------------------------------------------------------------------
! 1* construct the state between A and B
!    create Psi_AA'
    psi_aa = cmplx(0., 0.)
    do jj = 1, 4
        psi_aa = psi_aa + Pz_A(jj)*row(kronecker_product( column(signal(jj,:)), column(signal(jj,:)) ))
    enddo
!    creation of the density matrix rho_{AB} given by (eqn.1)
    rho_ab = Outer_product( psi_aa, conjg(psi_aa) )
    !!!! PAY ATTENTION THIS PASS NORMALIZE IN ANY CASE !!!!
    rho_AB = rho_AB / mat_trace(rho_AB)
    ! check if the density operator is physical
    sz = size(rho_AB,1)
    call checkpoint(real(mat_trace(rho_AB))-1 <= 1e-8, text="Tr(rho_AB)/=1",var=mat_trace(rho_AB))
    call checkpoint(is_hermitian(rho_AB,sz), text="rho_AB is not hermitian")
    call checkpoint(is_positive(rho_AB,sz), text="rho_AB is not positive")
    print*, "rho_ab purity", real(mat_trace(matmul( rho_ab, rho_ab )))

! 2* starting point of the algorithm
    do ii = 0, N_steps
        uu = dble(ii)/dble(N_steps)
        write(*,*) " - - - "
        write(*,*) " depo prob.", uu

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
            p_i(jj) = mat_trace( matmul( Gamma_i(jj,:,:), rho_ab) )
        enddo
        do jj = 1, size(p_tilde)
            p_tilde(jj) = mat_trace( matmul( Gamma_tilde_k(jj,:,:), rho_ab) )
        enddo

! ALGORITHM 
! 7* STEP 1
!      7.1* set counter to zero
        counter = 1
!      7.2* create rho_0
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
!      7.3* while counter <= Maxit
        do while( counter <= Maxit )
!      7.4* Calculate f(rho)
        ! apply the CP map G_e
            allocate(rho_4(size(isometry, 1), size(isometry, 1)), rho_5(size(isometry, 1), size(isometry, 1)))
            call CP_map(rho_0, Kraus, proj, isometry, rho_4, Ppass)
        ! apply pinching channel
            rho_5 = cmplx(0.,0.)
            do jj = 1, size(pinching, 1)
                rho_5 = rho_5 + matmul( matmul( pinching(jj,:,:), rho_4 ), pinching(jj,:,:) )
            enddo
            ! check if the density operator is physical
            sz = size(rho_5,1)
            call checkpoint(real(mat_trace(rho_5))-1 <= 1e-6, text="Tr(rho_5)/=1",var=mat_trace(rho_5))
            call checkpoint(is_hermitian(rho_5,sz),text="rho_5 is not hermitian")
            call checkpoint(is_positive(rho_5,sz),text="rho_5 is not positive")
            !if(all( abs(rho_4 - rho_5) <= 1e-4)) stop "rho_t_4 == rho_t_5"
            f_rho = RelativeEntropy(rho_4, rho_5, sz)
            write(*,*) " f(rho) =", f_rho

!      7.5* define gradient [grad_f(rho)]^T = G**+(log(G(rho))) - G**+(logZ(G(rho)))
            sz = size(rho_0, 1) ! 4
            ios = size(rho_4, 1)! 128
            allocate(grad_f(sz, sz)); grad_f = cmplx(0.,0.)
            allocate(logrho(ios, ios))
            allocate(rho_temp(sz, sz))
            logrho = logM(rho_4, ios)
            call CP_map_inverse( logrho, Kraus, proj, isometry, rho_temp)
            grad_f = rho_temp
            deallocate(rho_temp, logrho)
            allocate(logrho(ios,ios))
            allocate(rho_temp(sz, sz))
            logrho = logM(rho_5, ios)
            call CP_map_inverse( logrho, Kraus, proj, isometry, rho_temp)
            grad_f = transpose( grad_f - rho_temp )
            deallocate(logrho, rho_temp)

            !      7.6* minimize c^T X = sum_jw_j Tr[O_j^T graf_f]
            !          subject to AT <= B
            !   ATTENTION !!! VERIFY THAT THE MATRICES ARE REAL !!! OTHERWISE, TWICE THE DIMENSION OF THE PROBLEM  !!!
            sz = size(Omega_j, 1)
            ios = size(Omega_j, 2)
            allocate(AT(ios, ios), C(sz))
            ! find c^T = Tr[ O^T_j * grad_f ]
            C = 0.
            do jj = 1, sz
                C(jj) = real(mat_trace( matmul( transpose(  Omega_j(jj,:,:) ), grad_f) ), kind=8)! real
            enddo
            ! find AT constraints
            AT = 0.
            do jj = 1, sz
                AT = AT - real( Omega_j(jj,:,:) , kind=8)! the minus convert it in \geq
            enddo
            AT = transpose( real(AT + rho_0, kind=8) )! real
        ! define the output vector X solutions and Y dual solutions
            allocate(X_sol(sz), Y_sol(sz), IBASIS(sz), upper_bounds(sz))
            call FEASIBLEBASIS (sz, sz, AT, C, IBASIS, IERR)
            upper_bounds = 0.
        ! IBASIS = [(jj, jj=1,sz)]
            call DUALSIMPLEX(sz, sz, AT, upper_bounds, C, IBASIS, X_sol, Y_sol, IERR)
            ! call mat_dump(X_sol)
            ! call mat_dump(Y_sol)
        ! define Delta_rho = sum_j w_j O_j
            allocate(delta_rho(size(rho_0, 1), size(rho_0, 2)))
            delta_rho = cmplx(0., 0.)
            do jj = 1, size(Omega, 1)
                delta_rho = delta_rho + transpose( X_sol(jj)*Omega(jj, :, :) )
            enddo

!      7.8* convergence check
            print*, "WWWWWWW:", real(mat_trace(matmul(transpose(delta_rho), grad_f)))
            if(real(mat_trace(matmul(transpose(delta_rho), grad_f))) <= epsilon) then
                write(*,*) "algorithm exited at: ", counter, "step"
                DEALLOCATE(delta_rho, X_sol, Y_sol, IBASIS, upper_bounds, C, AT, grad_f, rho_4, rho_5, stat=ios)
                exit
            else if(counter == maxit) then
                write(*,*) "algorithm reached Maxit: ", Maxit
                DEALLOCATE(delta_rho, X_sol, Y_sol, IBASIS, upper_bounds, C, AT, grad_f, rho_4, rho_5, stat=ios)
                exit
            endif

!      7.9* find tt \in [0, 1]
            tt = 0.
            f_1 = f_rho
            sz = size(rho_4, 1)
            do jj = 0, finesse
                print*,jj
                rho_4 = cmplx(0., 0.); rho_5 = cmplx(0., 0.)    ! set var to zero
                call CP_map(rho_0 + dble(jj)/dble(finesse) * delta_rho, Kraus, proj, isometry, rho_4, Ppass)
                do kk = 1, size(pinching, 1)
                    rho_5 = rho_5 + matmul( matmul( pinching(kk,:,:), rho_4 ), pinching(kk,:,:) )
                enddo
                f_2 = RelativeEntropy(rho_4, rho_5, sz)
                if(f_2 < f_1) then
                    tt = dble(jj)/dble(finesse)
                    f_1 = f_2
                endif
            enddo
            print*, "fin"

!       7.10* assign new rho_0
            rho_0 = rho_0 + delta_rho * tt
!       7.11* increment counter
            counter = counter + 1
!       deallocation
            DEALLOCATE(delta_rho, X_sol, Y_sol, IBASIS, upper_bounds, C, AT, grad_f, rho_4, rho_5, stat=ios)
            call checkpoint(ios==0, text="while loop deallocation failed")
        enddo
! 8* STEP 2
        write(10,'(es20.10," ",es20.10)')uu, f_rho
        deallocate(rho_0, p_i, p_tilde, stat=ios)
        call checkpoint(ios==0, text="rho_0 deallocation failed")
    enddo
    write(10,'("EOD")')
    write(10,'(a)')"set terminal wxt 0"        
    write(10,'("set xlabel ''depolarization probability''")')
    write(10,'("set ylabel ''secret key rate''")')
    write(10,'("set title ''reliable lower bound for BB84 protocol''")')
    write(10,'("set grid")')
    write(10,'("unset key")')    
    write(10,'("#set xr [:]")')    
    write(10,'("plot $data u 1:2 w lp")')
    close(10)!close file
    call execute_command_line(("gnuplot -p "//trim(ps_file)))
!   END PROCEDURE
!---------------------------------------------------------------------------------------------------------------------

! deallocation
    DEALLOCATE(POVM_A, POVM_B, POVM, Gamma_i, Gamma_tilde, K_A, K_B, Kraus, isometry,&
                & proj, pinching, Omega_j, stat=ios)
    call checkpoint(ios==0, text="final deallocation failed")
    
! end computation
    call cpu_time(finish)
    WRITE(*, '(f10.3, A)') finish-start, " s"

endprogram QKDsimulation