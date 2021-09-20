!    MODULE: QKD 
!        *) subroutine: random_key
!        *) subroutine: depolarizing_channel
!        *) subroutine: binary_entropy
!        *) subroutine: complex_to_realm
!        *) subroutine: SDPA_write_problem
!        *) subroutine: SDPA_read_solution
!        *) subroutine: compute_primal
!
module QKD
    use debugging
    use matrices
    implicit none
    
    contains
    !random_key(
    !           strin(char(*)) (inout) = 'random string'
    !)
    subroutine random_key(string)
        implicit none
        !
        CHARACTER(*), INTENT(INOUT) :: string
        INTEGER(8)                  :: ii
        real                        :: uu
        INTEGER                     :: rnd
        !
        do ii = 1, len(string), 1
            call random_number(uu)
            rnd = floor(2*uu)
            write(string(ii:ii), '(i1)') rnd
        enddo
        return
    endsubroutine random_key

    !depolarizing_channel(
    !                    A(double complex) = 'qubit which is depolarized'
    !                    prob(real*4) = 'probability of being depolarized'
    !                   )
    subroutine depolarizing_channel(A, prob)
        implicit none
        !
        complex(8), DIMENSION(:,:), INTENT(INOUT)   :: A
        real(4)                                     :: prob
        A = prob/2 *identity(int8(size(A,1))) + (1-prob)*A
        return
    endsubroutine depolarizing_channel
    !angle_displacement(
    !                    A(double complex) = 'qubit which is depolarized'
    !                    angle(real*4) = 'angle theta which displace qubit polarization'
    !                   )
    subroutine angle_displacement(A, angle)
        implicit none
        !
        complex(8), DIMENSION(2,2), INTENT(INOUT)   :: A
        complex(8), DIMENSION(2,2)                  :: phase
        real(4)                                     :: angle
        !
        phase(1,:) = (/exp(cmplx(0.,1.)*angle), cmplx(0.,0.)/)
        phase(2,:) = (/cmplx(0.,0.), exp(-cmplx(0.,1.)*angle)/)
        A = matmul( phase, A)
        return
    endsubroutine angle_displacement

    subroutine CP_map(rho_input, Kraus_operators, projector, isometry, G_rho, Pprobability_pass)
        implicit none
        complex(8), INTENT(IN)                  :: rho_input(:,:), kraus_operators(:,:,:), projector(:,:), isometry(:,:)
        complex(8), INTENT(INOUT)               :: G_rho(:,:)
        real(4), INTENT(INOUT)                  :: Pprobability_pass
        complex(8), ALLOCATABLE                 :: isometry_H(:,:), rho_tilde(:,:)
        integer                                 :: siz, jj, ios

        allocate(rho_tilde( size(Kraus_operators,2), size(Kraus_operators,2) ))
        rho_tilde = cmplx(0.,0.)
        do jj = 1, size(Kraus_operators,1)
            rho_tilde = rho_tilde + matmul(Kraus_operators(jj,:,:), matmul(rho_input, conjg(transpose(Kraus_operators(jj,:,:)))))
        enddo
        ! rho_tilde = rho_tilde / mat_trace(rho_tilde)
        ! check if the density operator is physical
        siz = size(rho_tilde,1)
        ! call checkpoint(real(mat_trace(rho_tilde))-1 <= 1e-6, text="CP_map::1 Tr(rho_tilde)/=1",var=mat_trace(rho_tilde))
        ! call checkpoint(is_hermitian(rho_tilde,siz),text="CP_map::1 rho_tilde is not hermitian",flag=1)
        ! call checkpoint(is_positive(rho_tilde,siz),text="CP_map::1 rho_tilde is not positive",flag=1)
        ! apply the projectorector
        Pprobability_pass = real(mat_trace( matmul( rho_tilde, projector) ), kind=4)
        rho_tilde = matmul( projector, matmul(rho_tilde, projector))/Pprobability_pass
        ! rho_tilde = rho_tilde / mat_trace(rho_tilde)
        ! check if the density operator is physical
        siz = size(rho_tilde,1)
        ! call checkpoint(real(mat_trace(rho_tilde))-1 <= 1e-6, text="CP_map::2 Tr(rho_tilde)/=1",var=mat_trace(rho_tilde))
        ! call checkpoint(is_hermitian(rho_tilde,siz),text="CP_map::2 rho_tilde is not hermitian",flag=1)
        ! call checkpoint(is_positive(rho_tilde,siz),text="CP_map::2 rho_tilde is not positive",flag=1)
        ! isometry
        allocate(isometry_H(size(isometry,2), size(isometry,1) )); isometry_H = cmplx(0.,0.)
        isometry_H = conjg(transpose(isometry))

        ! apply the isometry V
        G_rho = matmul(matmul( isometry, rho_tilde ) , isometry_H)
        ! G_rho = G_rho / mat_trace(G_rho)
        ! ! check if the density operator is physical
        ! siz = size(G_rho,1)
        ! call checkpoint(real(mat_trace(G_rho))-1 <= 1e-6, text="CP_map:: Tr(G_rho)/=1",var=mat_trace(G_rho),flag=1)
        ! call checkpoint(is_hermitian(G_rho,siz),text="CP_map:: G_rho is not hermitian",flag=1)
        ! call checkpoint(is_positive(G_rho,siz,1e-2),text="CP_map:: G_rho is not positive",flag=1)
        DEALLOCATE(rho_tilde, isometry_H, stat=ios)
        call checkpoint(ios == 0, text = "CP_map:: deallocation failed. ")
        return
    endsubroutine CP_map

    subroutine CP_map_inverse(rho_input, Kraus_operators, projector, isometry, G_rho_inverse)
        implicit none
        complex(8), INTENT(IN)      :: rho_input(:,:), kraus_operators(:,:,:), projector(:,:), isometry(:,:)
        complex(8), INTENT(INOUT)   :: G_rho_inverse(:,:)
        real(4)                     :: Pprobability_pass
        complex(8), ALLOCATABLE     :: isometry_H(:,:), rho_tilde(:,:)
        integer                     :: jj, ios!, siz

        ! isometry
        allocate(isometry_H(size(isometry,2), size(isometry,1) )); isometry_H = cmplx(0.,0.)
        isometry_H = conjg(transpose(isometry))
        ! apply the isometry V
        allocate(rho_tilde(size(isometry_H,1),size(isometry_H,1)))
        ! call mat_dump( matmul(matmul( isometry_H, rho_input ) , isometry) )
        ! stop
        rho_tilde = matmul(matmul( isometry_H, rho_input ) , isometry)
        ! check if the density operator is physical
        ! siz = size(rho_tilde,1)
        ! call checkpoint(real(mat_trace(rho_tilde))-1 <= 1e-6, text="CP_map_inverse:: Tr(rho_tilde)/=1",var=mat_trace(rho_tilde))
        ! call checkpoint(is_hermitian(rho_tilde,siz),text="CP_map_inverse::1 rho_tilde is not hermitian",flag=1)
        ! call checkpoint(is_positive(rho_tilde,siz),text="CP_map_inverse::1 rho_tilde is not positive",flag=1)

        ! apply the projectorector
        Pprobability_pass = real(mat_trace( matmul( rho_tilde, projector) ), kind=4)
        rho_tilde = matmul( projector, matmul(rho_tilde, projector))/Pprobability_pass
        ! check if the density operator is physical
        ! siz = size(rho_tilde,1)
        ! call checkpoint(real(mat_trace(rho_tilde))-1 <= 1e-6, text="CP_map_inverse:: Tr(rho_tilde)/=1",var=mat_trace(rho_tilde))
        ! call checkpoint(is_hermitian(rho_tilde,siz),text="CP_map_inverse::2 rho_tilde is not hermitian",flag=1)
        ! call checkpoint(is_positive(rho_tilde,siz),text="CP_map_inverse::2 rho_tilde is not positive",flag=1)

        G_rho_inverse = cmplx(0.,0.)
        do jj = 1, size(Kraus_operators,1)
            G_rho_inverse = G_rho_inverse + matmul( conjg(transpose(Kraus_operators(jj,:,:))),&
                                                    & matmul(rho_tilde, Kraus_operators(jj,:,:)))
        enddo
        ! G_rho_inverse = G_rho_inverse / mat_trace(G_rho_inverse)

        ! ! check if the density operator is physical
        ! siz = size(G_rho_inverse,1)
        ! call checkpoint(real(mat_trace(G_rho_inverse))-1 <= 1e-6, text="CP_map_inverse:: Tr(G_rho_inverse)/=1",&
                        ! &var=mat_trace(G_rho_inverse))
        ! call checkpoint(is_hermitian(G_rho_inverse,siz),text="CP_map_inverse:: G_rho_inverse is not hermitian",flag=1)
        ! call checkpoint(is_positive(G_rho_inverse,siz),text="CP_map_inverse:: G_rho_inverse is not positive",flag=1)
        DEALLOCATE(rho_tilde, isometry_H, stat=ios)
        call checkpoint(ios == 0, text = "CP_map_inverse:: deallocation failed. ")
        return
    endsubroutine CP_map_inverse

    function VonNeumannEntropy(matrix, n, tol)
        implicit none
        integer, INTENT(IN)         :: n
        real, INTENT(IN), optional  :: tol
        complex(8), intent(in)      :: matrix(n,n)    ! input matrix 
        real(8), ALLOCATABLE        :: eigval(:)
        complex(8), ALLOCATABLE     :: mat(:,:), eigvec(:,:) 
        real(8)                     :: VonNeumannEntropy, tolerance, avl
        integer                     :: ii

        allocate(mat(n,n))
        mat = matrix
        call eigensolver(mat,eigval,eigvec,n)

        VonNeumannEntropy = 0.
        if(present(tol).eqv..true.) then; tolerance = tol; else; tolerance = 1e-10; endif
        
        do ii = 1, n
            avl = eigval(ii)
            if(avl <= 0.0)then
                if( abs(avl) > tolerance )then
                    WRITE(*,'(A,i0,A,f10.6)')"VonNeumannEntropy:: eigenvalue ",ii," is negative ",avl
                    stop
                endif
            else       
                VonNeumannEntropy = - avl*DLOG(avl) + VonNeumannEntropy
            endif
        enddo
        DEALLOCATE(eigval, mat, eigvec)
    endfunction VonNeumannEntropy

    function RelativeEntropy(rho, sigmap, n)
        implicit none
        integer, INTENT(IN)         :: n
        complex(8), intent(in)      :: rho(n,n), sigmap(n,n)    ! input rho
        real(8), ALLOCATABLE        :: p(:), q(:)
        complex(8), ALLOCATABLE     :: rho_avt(:,:), sigma_avt(:,:)
        real(8)                     :: RelativeEntropy, tol
        INTEGER                     ii

        ! print*,"re1"
        ! RelativeEntropy = VonNeumannEntropy(rho,n,1e-2)
        ! print*,"re2"
        ! RelativeEntropy = ( RelativeEntropy - real(mat_trace( matmul(rho, logM(sigmap,n) ) ), kind=8) )/log(2.)

        ! print*,"egv1"
        call eigensolver(rho, p, rho_avt, n)
        ! print*,"egv2"
        call eigensolver(sigmap, q, sigma_avt, n)

        tol = 1e-10
        RelativeEntropy = 0.
        do ii = 1, n
            if(p(ii)>tol.and.q(ii)>tol)then
                RelativeEntropy = RelativeEntropy + p(ii)*log(p(ii)/q(ii))/log(2.)
            endif
        enddo
    endfunction RelativeEntropy

    function binary_entropy(p) result(be)
    ! computes the binary entropy
    !       -p*log_2(p) - (1-p)*log_2(1-p)
    ! ---------------------------------------------------------------
    ! Arguemnts
    !   p = probability \in [0,1]
        implicit none
        real(8), INTENT(inout) :: p
        real(8)                :: be
        
        be = 0.
        if(abs(p) >= 1e-10 .and. abs(1.- p) >= 1e-10 ) then
            be = -p*log(p)/log(2.) - (1-p)*log(1-p)/log(2.)
        endif
    endfunction binary_entropy

    ! subroutine solve_sdp(r0, dr, gf, gb)
    !     ! r0 = rho_0
    !     ! dr = delta_rho
    !     ! gr = grad_f
    !     ! gb = gamma_basis
    !     implicit none
    
    !     complex(8), INTENT(IN)                      :: gb(:,:,:), gf(:,:), r0(:,:)
    !     complex(8), INTENT(INOUT)                   :: dr(:,:)
    !     real(8), ALLOCATABLE                        :: AT(:,:), C(:), sol_x(:), sol_y(:),&
    !     & up_bounds(:)
    !     integer, ALLOCATABLE                        ::  INTBASIS(:)
    !     integer nsize, ierror, jj
    
    !     !          subject to AT <= B
    !     !   ATTENTION !!! VERIFY THAT THE MATRICES ARE REAL !!! OTHERWISE, TWICE THE DIMENSION OF THE PROBLEM  !!!
    !     nsize = size(gb, 2)
    !     allocate(AT(nsize, nsize))
    !     nsize = size(gb, 1)
    !     allocate(C(nsize))
    !     ! find c^T = Tr[ O^T_j * grad_f ]
    !     C = 0.
    !     do jj = 1, nsize
    !         C(jj) = real(mat_trace( matmul( transpose( gb(jj,:,:) ), gf) ), kind=8)! real
    !     enddo
    !     ! find AT constraints
    !     AT = 0.
    !     do jj = 1, nsize
    !         AT = AT - real( gb(jj,:,:) , kind=8 )! the minus convert it in \geq
    !     enddo
    !     AT = transpose( real(AT + r0, kind=8) )! real
    !     ! define the output vector X solutions and Y dual solutions
    !     allocate(sol_x(nsize), sol_y(nsize), INTBASIS(nsize), up_bounds(nsize))
    !     call FEASIBLEBASIS (nsize, nsize, AT, C, INTBASIS, ierror)
    !     up_bounds = 0.
    !     ! INTBASIS = [(jj, jj=1,nsize)]
    !     call DUALSIMPLEX(nsize, nsize, AT, up_bounds, C, INTBASIS, sol_x, sol_y, ierror)
    !     call mat_dump(sol_x)
    !     call mat_dump(sol_y)
    !     ! define dr = sum_j w_j O_j
    !     dr = cmplx(0., 0.)
    !     do jj = 1, size(gb, 1)
    !         dr = dr +  sol_x(jj)*gb(jj, :, :) 
    !     enddo
    !     DEALLOCATE(at, c, sol_x, sol_y, up_bounds, INTBASIS, stat=ierror)
    ! endsubroutine solve_sdp

    subroutine complex_to_realm(n, cm, rm)
    ! Given an complex matrix cm \in C^{nxn}, this subroutine returns
    ! the real matrix rm \in R^{2nx2n}.
    ! ---------------------------------------------------------------
    ! Arguments:
    !   n = # of columns of om
    !   om(n,n) = complex input matrix
    !   nm(2n,2n) = real output matrix
        implicit none
        INTEGER, INTENT(IN)     :: n
        real(8), INTENT(INOUT)  :: rm(2*n,2*n)
        complex(8), INTENT(IN)  :: cm(n,n)
    
        rm = 0.
        rm(:n,:n) = real(cm)
        rm(:n,n+1:) = - aimag(cm)
        rm(n+1:,:n) = aimag(cm)
        rm(n+1:,n+1:) = real(cm)
        return
    endsubroutine complex_to_realm

    subroutine SDPA_write_problem(m, n, j, c, F0, Fi)
    ! Writes the SDP problem given by
    !   minimize  : sum_i( x_i*c_i )
    !   subject to: sum_i( x_i*F_i - F0 ) >= 0
    ! onto a text file sdp.dat and launches SDPA.exe in order to so-
    ! lve it.
    ! ---------------------------------------------------------------
    ! Arguments:
    !   m = number of c components
    !   n = matrix dimension of Gi
    !   j = length of Gi
    !   c = costant vector
    !   F0 = constant matrix 
    !   Fi = coefficient matrix 
    !   Gi = set of constraints operators generating constraints
    !   pi = set of constraints
        implicit none

        integer,intent(in) :: m            ! number of x_i
        integer,intent(in) :: n            ! dimension of F0
        integer,intent(in) :: j            ! dimension of F0
        real(8),INTENT(IN) :: c(m)         ! costant vector of c_i
        real(8),INTENT(IN) :: F0(n,n)       
        real(8),INTENT(IN) :: Fi(m,n,n)    
        ! real(8),INTENT(IN) :: Gi(j,n,n)    
        ! real(8),INTENT(IN) :: pi(j)    
        INTEGER            ii, jj, kk, sz

        ! open file
        open(unit=15, file="sdp.dat")

        ! write titleand comment
        write(15,'(A,A,A)') '"',"sdp problem",'"'
        write(15,'(i0,A)') m," = mDIM"          ! m
        ! write(15,'(i0,A)') j+2," = nBLOCK"       
        write(15,'(i0,A)') 1," = nBLOCK"       
        write(15,'(A,i0,A)',advance="no") "(",2*n,","
        do jj = 1, j
            write(15,'("1,")',advance="no")
        enddo
        write(15,'(A)',advance="no") "1) = bLOCKsTRUCT"   ! n

        ! write C
        sz = size(c)
        write(15,*)''
        write(15,'("{")',advance="no")
        do jj = 1, sz
            write(15, '(f7.4)' ,advance="no") c(jj)
            if(jj<sz)then
                write(15,'(",")',advance="no")
            endif
        enddo
        write(15,'("}")')

        ! write F0
        sz = size(F0,1)
        ! open
        write(15,'(" {")',advance="no")
        ! positivity
        do jj = 1, sz
            write(15,'(" {")',advance="no")
            do kk = 1, sz
                write(15, '(f7.4)' ,advance="no") F0(jj, kk)
                if(kk<sz)then
                    write(15,'(",")',advance="no")
                endif
            enddo
            write(15,'(" }")',advance="no")
            if(jj<sz)then
                write(15,'(",")',advance="no")
            endif
        enddo
        ! close 
        write(15,'("}")')

        ! write Fi
        do ii = 1, m
            ! open
            write(15,'("{")',advance="no")
            ! positivity
            do jj = 1, sz
                write(15,'(" {")',advance="no")
                do kk = 1, sz
                    write(15,'(es20.10)',advance="no") Fi(ii, jj, kk)
                    if(kk<sz)then
                        write(15,'(",")',advance="no")
                    endif
                enddo
                write(15,'("}")',advance="no")
                if(jj<sz)then
                    write(15,'(",")',advance="no")
                endif
            enddo
            ! close
            write(15,'("}")')
        enddo

        ! close sdp.out
        close(15)
        ! execute command to perform SDPA solving
        call execute_command_line("sdpa sdp.dat sdp.out", wait=.true.)
        ! stop
        return
    endsubroutine SDPA_write_problem
   
    subroutine SDPA_read_solution(m, arr)
    ! Reads the SDPA solution written in the result file sdp.out
    ! ---------------------------------------------------------------
    ! Arguments:
    !   m = number of c components
    !   arr = array where to store the solution
        implicit none

        integer, INTENT(IN)           :: m    ! dimension of x_i
        real(8), INTENT(INOUT)        :: arr(m)  ! solution
        CHARACTER(len=:), allocatable :: vals(:)
        CHARACTER(200)                line
        CHARACTER(:), ALLOCATABLE     :: stripped_line
        INTEGER                       ii, jj, kk, NUM_LINES

        ii = 0; kk = 0
        open(unit=20, file="sdp.out")
        ! num lines
        do while (kk == 0)
            NUM_LINES = NUM_LINES + 1
            read(20,'(A)', iostat=kk) line
            if(ii==1)then
                call strip_spaces(line,stripped_line)
                stripped_line = stripped_line(2:len(stripped_line)-1)
                call split(stripped_line, vals, ",")
                do jj = 1, size(vals)
                    read(vals(jj),*) arr(jj)
                enddo
                ! call mat_dump(arr)
                close(20); return
            endif
            if(line == "xVec = ")then; ii=1; else; ii=0; endif
        end do
        close(20)
        return
    endsubroutine SDPA_read_solution

    subroutine compute_primal(r0, fr, gf, kr, sf, is, km, Oj, mi, fi, ep)
    ! Compute_primal: computes the primal SDP problem
    ! ---------------------------------------------------------------
    ! Arguments:
    !   r0 = rho_0 the initial state
    !   fr = inout value f_rho = D(G(rho)||Z(G(rho)))
    !   gr = inout gradient of f_rho
    !   kr = kraus representation of enlarging register operation
    !   sf = sifting phase projector
    !   is = isometry V
    !   km = key map POVM
    !   Gi = Gamma_i
    !   pi = Tr[rho_0 @ Gamma_i] = p_i
    !   Oj = set of orthonormal Hermitian operators that complement the Gamma_i
    !   mi = maximum iteration
    !   fi = finesse of the tt search
    !   ep = tolerance of the algorithm    
        implicit none
    
        complex(8), INTENT(IN)                     :: r0(:,:)
        real(8), INTENT(INOUT)                     :: fr
        complex(8), INTENT(INOUT)                  :: gf(:,:), is(:,:), sf(:,:), kr(:,:,:),&
        &Oj(:,:,:), km(:,:,:)
        INTEGER, INTENT(IN), OPTIONAL              :: mi, fi
        REAL(8), INTENT(IN), OPTIONAL              :: ep
        complex(8), dimension(:,:), ALLOCATABLE    :: rho_i, rho_4, rho_5, delta_rho, logrho, rho_temp
        real(8), ALLOCATABLE                       :: oj_dbl(:,:,:), c_i(:), F0_real(:,:), x_i(:)
        INTEGER                                    :: mit, finesse, counter, m, j, siz, iostat, jj, tt, kk
        real                                       :: Ppass
        real(8)                                    :: f_1, f_2, epsilon
    
        if(present(mi))then
            mit = mi
        else 
            mit = 20
        endif
    
        if(present(fi))then
            finesse = fi
        else 
            finesse = 20
        endif
    
        if(present(ep))then
            epsilon = ep
        else 
            epsilon = 1E-10
        endif
        
        ! find the real matrices starting from the complex oj set
        m = size(oj,1)
        siz = size(oj,2)
        allocate(oj_dbl(m, 2*siz, 2*siz))
        do jj = 1, m
            call complex_to_realm(siz, oj(jj,:,:), oj_dbl(jj,:,:))
        enddo

        ! while counter <= mit
        counter = 1
        rho_i = r0
        do while( counter <= mit )
        ! Calculate f(rho)
        ! apply the CP map G_e
            allocate(rho_4(size(is, 1), size(is, 1)), rho_5(size(is, 1), size(is, 1)))
            call CP_map(rho_i, kr, sf, is, rho_4, Ppass)
        ! apply key map
            rho_5 = cmplx(0., 0.)
            do jj = 1, size(km, 1)
                rho_5 = rho_5 + matmul( matmul( km(jj,:,:), rho_4 ), km(jj,:,:) )
            enddo
            ! check if the density operator is physical
            siz = size(rho_5,1)
            call checkpoint(real(mat_trace(rho_5))-1 <= 1e-6, text="Tr(rho_5)/=1",var=mat_trace(rho_5))
            call checkpoint(is_hermitian(rho_5,siz),text="rho_5 is not hermitian")
            call checkpoint(is_positive(rho_5,siz),text="rho_5 is not positive")

        ! f_rho
            fr = RelativeEntropy(rho_4, rho_5, siz)
            write(*,*) " f(rho) =", fr

        ! define gradient [grad_f(rho)]^T = G**+(log[G(rho)]) - G**+(log[Z(G(rho))])
            gf = cmplx(0.,0.)
            siz = size(rho_4, 1)! 128
            allocate(logrho(siz,siz))
            siz = size(rho_i, 1) ! 4
            allocate(rho_temp(siz,siz))
            siz = size(rho_4, 1)
            logrho = logM(rho_4,siz)
            ! inverse function
            call CP_map_inverse(logrho, kr, sf, is, rho_temp)
            gf = rho_temp
            deallocate(rho_temp,logrho)
            siz = size(rho_4, 1)
            allocate(logrho(siz,siz))
            siz = size(rho_i, 1) ! 4
            allocate(rho_temp(siz,siz))
            siz = size(rho_4, 1)
            logrho = logM(rho_5,siz)
            ! inverse function
            call CP_map_inverse(logrho, kr, sf, is, rho_temp)
            gf = transpose(gf - rho_temp)/log(2.)
            deallocate(logrho, rho_temp)

        ! solve SDP
            siz = size(rho_i, 1)
            ALLOCATE(c_i(m))
            do jj = 1, m
                c_i(jj) = real(mat_trace(matmul(oj(jj,:,:), gf)),4)
            enddo
            allocate(F0_real(2*siz,2*siz))
            call complex_to_realm(siz, -rho_i, F0_real)
            call SDPA_write_problem(m, 2*siz, j, c_i, F0_real, oj_dbl)
            allocate(x_i(m))
            call SDPA_read_solution(m, x_i)

            ! find delta_rho
            siz = size(rho_i, 1)
            allocate(delta_rho(siz, siz))
            delta_rho = cmplx(0.,0.)
            do jj = 1, m
                delta_rho = delta_rho + x_i(jj)*Oj(jj,:,:)
            enddo
            ! call solve_sdp(rho_i, delta_rho, gf, Oj)

        ! convergence check
            if(abs(real(mat_trace(matmul(transpose(delta_rho),gf)))).le.epsilon) then
                write(*,*) "algorithm exited at: ", counter, "step (precision=", &
                &real(mat_trace(matmul(transpose(delta_rho), gf))),")"
                DEALLOCATE(delta_rho,rho_i,rho_4,rho_5,x_i,oj_dbl,c_i,F0_real,stat=iostat)
                exit
            else if(counter==mit) then
                write(*,*) "algorithm reached maxit: ", mit
                DEALLOCATE(delta_rho,rho_i,rho_4,rho_5,x_i,oj_dbl,c_i,F0_real,stat=iostat)
                exit
            endif

        ! find tt \in [0, 1]
            tt = 0.
            f_1 = fr
            siz = size(rho_4, 1)
            do jj = 0, finesse
                rho_4 = cmplx(0., 0.); rho_5 = cmplx(0., 0.)    ! set var to zero
                call CP_map((rho_i + dble(jj)/dble(finesse)*delta_rho), kr, sf, is, rho_4, Ppass)
                do kk = 1, size(km, 1)
                    rho_5 = rho_5 + matmul( matmul( km(kk,:,:), rho_4 ), km(kk,:,:) )
                enddo
                f_2 = RelativeEntropy(rho_4, rho_5, siz)
                if(f_2<=f_1) then
                    tt = int(dble(jj)/dble(finesse))
                    f_1 = f_2
                endif
            enddo

        ! if f_1 == fr EXIT
            if(abs(f_1-fr) <= 1e-10) exit
        ! assign new rho_i
            rho_i = rho_i + delta_rho * tt
        ! increment counter
            counter = counter + 1
        ! deallocation
            DEALLOCATE(delta_rho, rho_4, rho_5, x_i, c_i, F0_real, stat=iostat)
            call checkpoint(iostat==0, text="while loop deallocation failed")
        enddo
        DEALLOCATE(oj_dbl, rho_i, stat=iostat)
        return
    endsubroutine compute_primal

endmodule QKD




! subroutine SDPA_write_problem(m, n, c, F0, Fi, Gi, pi)
!     ! Writes the SDP problem given by
!     !   minimize  : sum_i( x_i*c_i )
!     !   subject to: sum_i( x_i*F_i - F0 ) >= 0
!     ! onto a text file sdp.dat and launches SDPA.exe in order to so-
!     ! lve it.
!     ! ---------------------------------------------------------------
!     ! Arguments:
!     !   m = number of c components
!     !   n = dimension of Fi
!     !   c = costant vector
!     !   F0 = constant matrix 
!     !   Fi = coefficient matrix 
!     !   Gi = set of constraints operators generating constraints
!     !   pi = set of constraints
!         implicit none

!         integer,intent(in)    :: m            ! number of x_i
!         integer,intent(in)    :: n            ! dimension of F0
!         real(8),INTENT(IN)    :: c(m)         ! costant vector of c_i
!         real(8),INTENT(IN)    :: F0(n,n)       
!         real(8),INTENT(IN)    :: Fi(n,n,n)    
!         complex(8),INTENT(IN) :: Gi(:,:,:)    
!         real(8),INTENT(IN)    :: pi(:)    
!         INTEGER               ii, jj, kk, sz

!         ! open file
!         open(unit=15, file="sdp.dat")

!         ! write titleand comment
!         write(15,'(A)') '"sdp problem"'
!         write(15,'(i0,A)') m," = mDIM"          ! m
!         write(15,'(A)') "1 = nBLOCK"            ! 
!         write(15,'(i0,A)') n," = bLOCKsTRUCT"   ! n

!         ! write C
!         sz = size(c)
!         write(15,'("{")',advance="no")
!         do jj = 1, sz
!             write(15, '(f7.4)' ,advance="no") c(jj)
!             if(jj<sz)then
!                 write(15,'(",")',advance="no")
!             endif
!         enddo
!         write(15,'("}")')

!         ! write F0
!         sz = size(F0,1)
!         write(15,'("{")',advance="no")
!         do jj = 1, sz
!             write(15,'("{")',advance="no")
!             do kk = 1, sz
!                 write(15, '(f7.4)' ,advance="no") F0(jj, kk)
!                 if(kk<sz)then
!                     write(15,'(",")',advance="no")
!                 endif
!             enddo
!             write(15,'("}")',advance="no")
!             if(jj<sz)then
!                 write(15,'(",")',advance="no")
!             endif
!         enddo
!         write(15,'("}")',advance="no")

!         ! write Fi
!         do ii = 1, m
!             write(15,*)''
!             write(15,'("{")',advance="no")
!             do jj = 1, sz
!                 write(15,'("{")',advance="no")
!                 do kk = 1, sz
!                     write(15, '(f7.4)' ,advance="no") Fi(ii, jj, kk)
!                     if(kk<sz)then
!                         write(15,'(",")',advance="no")
!                     endif
!                 enddo
!                 write(15,'("}")',advance="no")
!                 if(jj<sz)then
!                     write(15,'(",")',advance="no")
!                 endif
!             enddo
!             write(15,'("}")',advance="no")
!             if(jj<sz)then
!                 write(15,'(",")',advance="no")
!             endif
!         enddo

!         ! close sdp.out
!         close(15)

!         ! execute command to perform SDPA solving
!         call execute_command_line("sdpa sdp.dat sdp.out", wait=.true.)
!         return
!     endsubroutine SDPA_write_problem