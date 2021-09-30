!    MODULE: QKD 
!        *) subroutine: random_key
!        *) subroutine: depolarizing_channel
!        *) subroutine: VonNeumannEntropy
!        *) subroutine: RelativeEntropy
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

    subroutine random_key(string)
    !random_key(
    !           strin(char(*)) (inout) = 'random string'
    !)
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

    subroutine depolarizing_channel(A, prob)
    !depolarizing_channel(
    !                    A(double complex) = 'qubit which is depolarized'
    !                    prob(real*4) = 'probability of being depolarized'
    !                   )
        implicit none
        !
        complex(8), DIMENSION(:,:), INTENT(INOUT)   :: A
        real(4)                                     :: prob
        A = prob/2 *identity(int8(size(A,1))) + (1-prob)*A
        return
    endsubroutine depolarizing_channel
    
    subroutine angle_displacement(A, angle)
    !angle_displacement(
    !                    A(double complex) = 'qubit which is depolarized'
    !                    angle(real*4) = 'angle theta which displace qubit polarization'
    !                   )
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

    subroutine CP_map(ri, kr, sf, is, ro, pp)
    ! G map
    ! Args:
    !   ri = rho_input
    !   kr = kraus operators
    !   sf = sifting phase projector
    !   is = isometry V
    !   ro = rho_output
    !   pp = probability of passing
        implicit none
        complex(8),INTENT(IN)       :: ri(:,:), kr(:,:,:), sf(:,:), is(:,:)
        complex(8),INTENT(INOUT)    :: ro(:,:)
        real(4),INTENT(INOUT)       :: pp
        complex(8),ALLOCATABLE      :: rho_tilde(:,:)
        integer                     :: jj, ios

        ! apply enlarging
        allocate(rho_tilde( size(kr,2), size(kr,2) ))
        rho_tilde = cmplx(0.,0.)
        do jj = 1, size(kr,1)
            rho_tilde=rho_tilde+matmul(matmul(kr(jj,:,:),ri),conjg(transpose(kr(jj,:,:))))
        enddo
        
        ! apply the sifting phase
        rho_tilde = matmul( rho_tilde, sf)
        pp = real(mat_trace( rho_tilde ), kind=4)
        rho_tilde = matmul( sf, rho_tilde )/pp

        ! apply the is V
        ro = matmul(matmul( is, rho_tilde ), conjg(transpose(is)))
        DEALLOCATE(rho_tilde, stat=ios)
        call checkpoint(ios == 0, text = "CP_map:: deallocation failed. ")

        return
    endsubroutine CP_map

    subroutine CP_map_inverse(rho_input, Kraus_operators, projector, isometry, G_rho_inverse)
    ! inverse G map
        implicit none
        complex(8), INTENT(IN)      :: rho_input(:,:), kraus_operators(:,:,:), projector(:,:), isometry(:,:)
        complex(8), INTENT(INOUT)   :: G_rho_inverse(:,:)
        !real(4)                     :: Pprobability_pass
        complex(8), ALLOCATABLE     :: rho_tilde(:,:)
        integer                     :: jj, ios!, siz

        ! apply the isometry V
        rho_tilde = matmul(matmul( conjg(transpose(isometry)), rho_input ) , isometry)

        ! apply the projectorector
        !Pprobability_pass = real(mat_trace( matmul( rho_tilde, projector) ), kind=4)
        rho_tilde = matmul( projector, matmul(rho_tilde, projector))!/Pprobability_pass

        ! apply inverse enlarging
        G_rho_inverse = cmplx(0.,0.)
        do jj = 1, size(Kraus_operators,1)
            G_rho_inverse = G_rho_inverse + matmul( conjg(transpose(Kraus_operators(jj,:,:))),&
                            & matmul(rho_tilde, Kraus_operators(jj,:,:)))
        enddo

        DEALLOCATE(rho_tilde, stat=ios)
        call checkpoint(ios == 0, text = "CP_map_inverse:: deallocation failed. ")
        return
    endsubroutine CP_map_inverse

    function VonNeumannEntropy(rho, n)
    ! Computes the von Neumann entropy of a given quantum state.
    !----------------------------------------------------------------
    ! args
    !   n = # columns rho
    !   rho = density matrix
        use ieee_arithmetic
        implicit none
        integer, INTENT(IN)         :: n
        complex(8), intent(in)      :: rho(n,n)    ! input rho 
        complex(8)                  :: new_rho(n,n)
        real(8)                     :: VonNeumannEntropy, fudge
        fudge = 1e-15
    
        new_rho = (1-fudge) * rho + fudge * identity(int8(n))

        new_rho = logmz(new_rho,n)

        ! print*,"vv logmz",mat_trace(new_rho)
        ! print*,"vv      ",mat_trace(matmul(rho,new_rho))
        if (isnan(real(mat_trace(rho),4))) stop " logmz nan"
        VonNeumannEntropy=-1 *real(mat_trace(matmul(rho,new_rho)))
    endfunction VonNeumannEntropy

    function RelativeEntropy(rho1, rho2, n)
    ! Computes the quantum relative entropy D(rho1 || rho2)
    ! ---------------------------------------------------------------
    ! args:
    !   n = dimension of rho1, rho2
    !   rho1
    !   rho2
        implicit none
        integer, INTENT(IN)         :: n
        complex(8), intent(in)      :: rho1(n,n), rho2(n,n)    ! input rho
        complex(8)                  :: new_rho2(n,n)
        real(8)                     :: RelativeEntropy, fudge
        
        fudge = 1e-16
        if(size(rho1) /= size(rho2)) stop "Need two matrices to be the same shape"
    
        new_rho2 = (1-fudge) * rho2 + fudge * identity(int8(n))
    
        RelativeEntropy = -1*VonNeumannEntropy(rho1,n)
        ! print*,"re2",RelativeEntropy
        RelativeEntropy = real( RelativeEntropy - mat_trace(matmul(rho1, logmz(new_rho2,n)))) /log(2.)
        ! print*,"re3",RelativeEntropy

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


    subroutine SDPA_write_problem1(m, j, n, c, F0, Fi, Gi, pi)
    ! Writes the SDP problem given by
    !   minimize  : sum_i( x_i*c_i )
    !   subject to: sum_i( x_i*F_i - F0 ) >= 0
    !               Tr[(sum_j x_j*F_j - F0)@Gi] = pi 
    ! onto a text file sdp.dat and launches SDPA.exe in order to so-
    ! lve it.
    ! ---------------------------------------------------------------
    ! Arguments:
    !   m = number of c components
    !   j = number of pi components
    !   n = dimension of Fi
    !   c = costant vector
    !   F0 = constant matrix 
    !   Fi = coefficient matrix 
    !   Gi = set of constraints operators generating constraints
    !   pi = set of constraints
        implicit none

        integer,intent(in)    :: m            ! number of x_i
        integer,intent(in)    :: j            ! number of x_i
        integer,intent(in)    :: n            ! dimension of F0
        real(8),INTENT(IN)    :: c(m)         ! costant vector of c_i
        real(8),INTENT(IN)    :: F0(n,n)      ! it is rho_0 from complex to real
        real(8),INTENT(IN)    :: Fi(m,n,n)    
        real(8),INTENT(IN)    :: Gi(j,n,n)    
        real(8),INTENT(IN)    :: pi(j)    
        INTEGER               ii, jj, kk
        INTEGER,PARAMETER     :: uf=15      !file unit

        ! open file
        open(unit=uf, file="sdp.dat", action="write")

        ! write titleand comment
        write(uf,'(A)') '"sdp problem"'
        write(uf,'(i0,A)') m," = mDIM"          ! m
        write(uf,'(i0,A)') 1+j," = nBLOCK"  ! j-gamma constraints, 1 trace constr, 1 positivity constr, 1 hermiticity constr    
        write(uf,'(A,i0,A)',advance="no") "(",n,","   ! positivity (-n because it is antisymmetric)
        do ii = 1, j-1
            write(uf,'("1,")',advance="no")            ! gamma_i
        enddo
        write(uf,'("1) = bLOCKsTRUCT")')               ! tr = 1

        ! write C
        write(uf,'("{")',advance="no")
        do jj = 1, m
            write(uf, '("",es20.10,"")' ,advance="no") c(jj)
            if(jj<m)then
                write(uf,'(",")',advance="no")
            endif
        enddo
        write(uf,'("}")')

        ! write F0
        ! open
        write(uf,'("{")')
        ! positivty
        write(uf,'(" {")',advance="no")
        do jj = 1, n
            write(uf,'(" {")',advance="no")
            do kk = 1, n
                write(uf, '(" {",es20.10,"}")',advance="no") F0(kk,jj)
                if(kk<n)then
                    write(uf,'(",")',advance="no")
                endif
            enddo
            if(jj<n)then
                write(uf,'(" },")')
            else
                write(uf,'(" }")',advance="no")
            endif
        enddo
        write(uf,'(" }")')
        ! tr[rho_i Gj]
        do ii = 1, j!'(" {",f10.7,"}")'
            write(uf,'(" ",es20.10,"")')(-real(pi(ii))+real(mat_trace(matmul(F0,Gi(ii,:,:)))))
        enddo
        ! tr = 1
        write(uf,'(" ",es20.10,"")') (-1.+real(mat_trace(F0)))
        ! close
        write(uf,'("}")')

        ! write Fi
        do ii = 1, m
            ! open
            write(uf,'("{")')
            ! Omega_i
            write(uf,'(" {")',advance="no")
            do jj = 1, n
                write(uf,'(" {")',advance="no")
                do kk = 1, n
                    write(uf, '(" {",es20.10,"}")' ,advance="no") Fi(ii,kk,jj)
                    if(kk<n)then
                        write(uf,'(",")',advance="no")
                    endif
                enddo
                if(jj<n)then
                    write(uf,'(" },")')
                else
                    write(uf,'(" }")',advance="no")
                endif
            enddo
            write(uf,'(" }")')
            ! x_k Tr[Omega_k Gj]
            do jj = 1, j
                write(uf,'(" ",es20.10,"")') real(mat_trace(matmul(Fi(ii,:,:),Gi(jj,:,:))))
            enddo
            ! tr = 1
            write(uf,'(" ",es20.10,"")') real(mat_trace(Fi(ii,:,:)))
            ! close
            write(uf,'("}")')
        enddo

        ! close sdp.out
        close(uf)

        ! execute command to perform SDPA solving
        call execute_command_line("sdpa sdp.dat sdp.out param.sdpa", wait=.true.)
        return
    endsubroutine SDPA_write_problem1
   
    subroutine SDPA_write_problem(m, n, c, F0, Fi)
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
        implicit none

        integer,intent(in) :: m            ! number of x_i
        integer,intent(in) :: n            ! dimension of F0
        real(8),INTENT(IN) :: c(m)         ! costant vector of c_i
        real(8),INTENT(IN) :: F0(n,n)       
        real(8),INTENT(IN) :: Fi(m,n,n)    
        INTEGER            ii, jj, kk, sz, uf
        parameter          (uf=15)

        ! open file
        open(unit=uf, file="sdp.dat", action="write")

        ! write titleand comment
        write(uf,'(A,A,A)') '"',"sdp problem",'"'
        write(uf,'(i0,A)') m," = mDIM"          
        write(uf,'(i0,A)') 1," = nBLOCK"       
        write(uf,'(i0," = bLOCKsTRUCT")') n

        ! write C
        sz = size(c)
        write(uf,'("{")',advance="no")
        do jj = 1, sz
            write(uf, '(f10.7)' ,advance="no") c(jj)
            if(jj<sz)then
                write(uf,'(",")',advance="no")
            endif
        enddo
        write(uf,'("}")')

        ! write F0
        sz = size(F0,1)
        ! open
        write(uf,'("{")',advance="no")
        ! positivity
        do jj = 1, sz
            write(uf,'(" {")',advance="no")
            do kk = 1, sz
                write(uf,'(f10.7)',advance="no") F0(kk,jj)
                if(kk<sz)then
                    write(uf,'(",")',advance="no")
                endif
            enddo
            write(uf,'(" }")',advance="no")
            if(jj<sz)then
                write(uf,'(",")',advance="no")
            endif
        enddo
        ! close 
        write(uf,'("}")')

        ! write Fi
        do ii = 1, m
            ! open
            write(uf,'("{")',advance="no")
            ! positivity
            do jj = 1, sz
                write(uf,'(" {")',advance="no")
                do kk = 1, sz
                    write(uf,'(f10.7)',advance="no")Fi(ii,kk,jj)
                    if(kk<sz)then
                        write(uf,'(",")',advance="no")
                    endif
                enddo
                write(uf,'("}")',advance="no")
                if(jj<sz)then
                    write(uf,'(",")',advance="no")
                endif
            enddo
            ! close
            write(uf,'("}")')
        enddo

        ! close sdp.out
        close(uf)
        ! execute command to perform SDPA solving
        call execute_command_line("sdpa sdp.dat sdp.out",wait=.true.)
        ! stop
        return
    endsubroutine SDPA_write_problem 

    subroutine solve_sdpa(m, j, n, c, F0, Fi, Gi, pi)
    ! Writes the SDP problem given by
    !   minimize  : sum_i( x_i*c_i )
    !   subject to: sum_i( x_i*F_i - F0 ) >= 0
    !               Tr[(sum_j x_j*F_j - F0)@Gi] = pi 
    ! onto a text file sdp.dat and launches SDPA.exe in order to so-
    ! lve it.
    ! ---------------------------------------------------------------
    ! Arguments:
    !   m = number of c components
    !   j = number of pi components
    !   n = dimension of Fi
    !   c = costant vector
    !   F0 = constant matrix 
    !   Fi = coefficient matrix 
    !   Gi = set of constraints operators generating constraints
    !   pi = set of constraints
        implicit none

        integer,intent(in)  :: m            ! number of x_i
        integer,intent(in)  :: j            ! number of x_i
        integer,intent(in)  :: n            ! dimension of F0
        real(8),INTENT(IN)  :: c(m)         ! costant vector of c_i
        real(8),INTENT(IN)  :: F0(n,n)      ! it is rho_0 from complex to real
        real(8),INTENT(IN)  :: Fi(m,n,n)    
        real(8),INTENT(IN)  :: Gi(j,n,n)    
        real(8),INTENT(IN)  :: pi(j)    
        INTEGER             ii, jj
        INTEGER,PARAMETER   :: uf=15      !file unit

        ! open file
        open(unit=uf, file="sdp.dat", action="write")

        ! write titleand comment
        write(uf,'(A)') '"sdp problem"'
        write(uf,'(i0,A)') m," = mDIM"          ! m
        write(uf,'(i0,A)') j+1," = nBLOCK"  ! j-gamma constraints, 1 trace constr, 1 positivity constr, 1 hermiticity constr    
        write(uf,'(A)',advance="no") "("   ! positivity (-n because it is antisymmetric)
        do ii = 1, j
            write(uf,'("-1,")',advance="no")            ! gamma_i
        enddo
        write(uf,'("-1) = bLOCKsTRUCT")')               ! tr = 1

        ! write C
        write(uf,'("{")',advance="no")
        do jj = 1, m
            write(uf, '("",es20.10,"")' ,advance="no") c(jj)
            if(jj<m)then
                write(uf,'(",")',advance="no")
            endif
        enddo
        write(uf,'("}")')

        ! write F0
        ! open
        write(uf,'("{")')
        ! tr[rho_i Gj]
        do ii = 1, j!'(" {",f10.7,"}")'
            write(uf,'(" ",es20.10,"")')(-real(pi(ii))+real(mat_trace(matmul(F0,Gi(ii,:,:)))))
        enddo
        ! tr = 1
        write(uf,'(" ",es20.10,"")') (-1.+real(mat_trace(F0)))
        ! close
        write(uf,'("}")')

        ! write Fi
        do ii = 1, m
            ! open
            write(uf,'("{")')
            ! x_k Tr[Omega_k Gj]
            do jj = 1, j
                write(uf,'(" ",es20.10,"")') real(mat_trace(matmul(Fi(ii,:,:),Gi(jj,:,:))))
            enddo
            ! tr = 1
            write(uf,'(" ",es20.10,"")') real(mat_trace(Fi(ii,:,:)))
            ! close
            write(uf,'("}")')
        enddo

        ! close sdp.out
        close(uf)

        ! execute command to perform SDPA solving
        call execute_command_line("sdpa sdp.dat sdp.out",wait=.true.)
        return
    endsubroutine solve_sdpa

    subroutine solve_sparse_sdpa(m, j, n, c, F0, Fi, Gi, pi)
    ! Writes the SDP problem given by
    !   minimize  : sum_i( x_i*c_i )
    !   subject to: sum_i( x_i*F_i - F0 ) >= 0
    !               Tr[(sum_j x_j*F_j - F0)@Gi] = pi 
    ! onto a text file sdp.dat and launches SDPA.exe in order to so-
    ! lve it.
    ! ---------------------------------------------------------------
    ! Arguments:
    !   m = number of c components
    !   j = number of pi components
    !   n = dimension of Fi
    !   c = costant vector
    !   F0 = constant matrix 
    !   Fi = coefficient matrix 
    !   Gi = set of constraints operators generating constraints
    !   pi = set of constraints
        implicit none

        integer,intent(in)  :: m            ! number of x_i
        integer,intent(in)  :: j            ! number of x_i
        integer,intent(in)  :: n            ! dimension of F0
        real(8),INTENT(IN)  :: c(m)         ! costant vector of c_i
        real(8),INTENT(IN)  :: F0(n,n)      ! it is rho_0 from complex to real
        real(8),INTENT(IN)  :: Fi(m,n,n)    
        real(8),INTENT(IN)  :: Gi(j,n,n)    
        real(8),INTENT(IN)  :: pi(j)    
        INTEGER             ii, jj
        INTEGER,PARAMETER   :: uf=15      !file unit

        ! open file
        open(unit=uf, file="sdp.dat-s", action="write")

        ! write titleand comment
        write(uf,'(A)') '"sdp problem"'
        write(uf,'(i0,A)') m," = mDIM"         ! m
        write(uf,'("1 = nBLOCK")')             ! j-gamma constraints, 1 trace constr, 1 positivity constr, 1 hermiticity constr    
        write(uf,'(i0,A)') -(m+1)," = bLOCKsTRUCT"

        ! write C
        write(uf,'("{")',advance="no")
        do jj = 1, m
            write(uf, '("",es20.10,"")' ,advance="no") c(jj)
            if(jj<m)then
                write(uf,'(",")',advance="no")
            endif
        enddo
        write(uf,'("}")')
        
        ! write F0
        ! tr[rho_i Gj]
        do ii = 1, j!'(" {",f10.7,"}")'
            write(uf,'("0   1   ",i0,"   ",i0,"   ",f20.10)')ii,ii,&
            &(-real(pi(ii))+real(mat_trace(matmul(F0,Gi(ii,:,:)))))
        enddo

        ! write Fi
        do ii = 1, m
            ! open
            ! x_k Tr[Omega_k Gj]
            do jj = 1, j
                write(uf,'(i0,"   1   ",i0,"   ",i0,"   ",f20.10)')ii,jj,jj,&
                &real(mat_trace(matmul(Fi(ii,:,:),Gi(jj,:,:))))
            enddo
            ! close
        enddo

        ! close sdp.out
        close(uf)

        ! execute command to perform SDPA solving
        call execute_command_line("sdpa sdp.dat-s sdp.out",wait=.true.)
        return
    endsubroutine solve_sparse_sdpa

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
            open(unit=20, file="sdp.out", action='READ')
    
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

    subroutine compute_primal(r0, fr, gf, kr, sf, is, km, Oj, Gi, pi, mi, fi, ep)
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
    
        complex(8),INTENT(IN)                  :: r0(:,:)
        real(8),INTENT(OUT)                    :: fr
        complex(8),INTENT(OUT)                 :: gf(:,:)
        complex(8),intent(In)                  :: is(:,:),sf(:,:),kr(:,:,:),Oj(:,:,:),km(:,:,:),Gi(:,:,:)
        real(8),INTENT(IN)                     :: pi(:)
        INTEGER,INTENT(IN),OPTIONAL            :: mi, fi
        REAL(8),INTENT(IN),OPTIONAL            :: ep
        complex(8),dimension(:,:),ALLOCATABLE  :: rho_i, rho_4, rho_5, delta_rho, logrho, rho_temp
        real(8),ALLOCATABLE                    :: oj_real(:,:,:),c_i(:),F0_real(:,:),x_i(:),Gi_real(:,:,:)
        INTEGER                                :: mit, finesse, counter, m, j, siz, iostat, jj, kk
        real(8)                                :: f_1, f_2, uu, tt, epsilon
        real                                   :: Ppass
    
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
        allocate(oj_real(m, 2*siz, 2*siz))
        do jj = 1, m
            call complex_to_realm(siz, oj(jj,:,:), oj_real(jj,:,:))
        enddo

        ! find the real matrices of Gamma_i
        j = size(gi,1)
        siz = size(gi,2)
        allocate(gi_real(j, 2*siz, 2*siz))
        do jj = 1, j
            call complex_to_realm(siz, gi(jj,:,:), gi_real(jj,:,:))
        enddo

        ! while counter <= mit
        counter = 1
        rho_i = r0
        do while( counter <= mit )
        ! Calculate f(rho)
        ! apply the CP map G_e
            allocate(rho_4(size(is, 1), size(is, 1)), rho_5(size(is, 1), size(is, 1)))
            call CP_map(rho_i, kr, sf, is, rho_4, Ppass)
            rho_4 = rho_4/mat_trace(rho_4)
        ! apply key map
            rho_5 = cmplx(0., 0.)
            do jj = 1, size(km, 1)
                rho_5 = rho_5 + matmul(matmul(km(jj,:,:),rho_4),km(jj,:,:))
            enddo
            rho_5 = rho_5/ real(mat_trace(rho_5))
            ! check if the density operator is physical
            siz = size(rho_5,1)
            ! call checkpoint(real(mat_trace(rho_5))-1 <= 1e-6, text="compute_primal::Tr(rho_5)/=1",var=mat_trace(rho_5))
            ! call checkpoint(is_hermitian(rho_5,siz),text="compute_primal::rho_5 is not hermitian")
            ! call checkpoint(is_positive(rho_5,siz,1e-5),text="compute_primal::rho_5 is not positive")

        ! f_rho
            fr = RelativeEntropy(rho_4, rho_5, siz)!22636!*ppass

        ! define gradient [grad_f(rho)]^T = G**+(log[G(rho)]) - G**+(log[Z(G(rho))])
            gf = cmplx(0.,0.)
            siz = size(rho_4, 1)! 128
            allocate(logrho(siz,siz))
            siz = size(rho_i, 1) ! 4
            allocate(rho_temp(siz,siz))
            siz = size(rho_4, 1)
            logrho = logmz(rho_4,siz)
            ! inverse function
            call CP_map_inverse(logrho, kr, sf, is, rho_temp)
            gf = rho_temp
            deallocate(rho_temp,logrho)
            siz = size(rho_4, 1)
            allocate(logrho(siz,siz))
            siz = size(rho_i, 1) ! 4
            allocate(rho_temp(siz,siz))
            siz = size(rho_5, 1)
            logrho = logmz(rho_5,siz)
            ! inverse function
            call CP_map_inverse(logrho, kr, sf, is, rho_temp)
            gf = (gf - rho_temp)/log(2.)
            deallocate(logrho, rho_temp)

        ! solve SDP
            siz = size(rho_i, 1)
            ALLOCATE(c_i(m))
            do jj = 1, m
                c_i(jj) = real(mat_trace(matmul(oj(jj,:,:),gf)))
            enddo
            allocate(F0_real(2*siz,2*siz))
            call complex_to_realm(siz,-rho_i,F0_real)
            ! call SDPA_write_problem(m,2*siz,c_i,F0_real,Oj_real)
            ! call SDPA_write_problem1(m,j,2*siz,c_i,F0_real,Oj_real,Gi_real,pi)
            call solve_sparse_SDPA(m,j,2*siz,c_i,F0_real,Oj_real,Gi_real,pi)!solve_SDPA
            allocate(x_i(m))
            call SDPA_read_solution(m,x_i)
            write(*,*) " the value at iteration ", counter
            write(*,*) " f(rho) =", fr

        ! find delta_rho
            siz = size(rho_i,1)
            allocate(delta_rho(siz, siz))
            delta_rho = cmplx(0.,0.)
            do jj = 1, m
                delta_rho = delta_rho + x_i(jj)*Oj(jj,:,:)
            enddo
            ! normalize
            ! delta_rho = delta_rho/norm2(x_i)
            print*, "norm xi ",norm2(x_i)
            print*, "sum pi ",sum(pi)

        ! convergence check
            if(abs(real(mat_trace(matmul(delta_rho,gf)))).le.epsilon) then
                write(*,*)"algorithm exited at: ",counter,"step (precision=",&
                &real(mat_trace(matmul(delta_rho,gf))),")"
                DEALLOCATE(delta_rho,rho_i,rho_4,rho_5,x_i,oj_real,c_i,F0_real,stat=iostat)
                return
            elseif(counter==mit) then
                write(*,*) "algorithm reached maxit: ", mit
                DEALLOCATE(delta_rho,rho_i,rho_4,rho_5,x_i,oj_real,c_i,F0_real,stat=iostat)
                return
            endif

        ! find tt \in [0, 1]
            tt = 0.
            f_1 = fr
            siz = size(rho_4, 1)
            do jj = 1, finesse
                uu = dble(jj)/dble(finesse)
                rho_4 = cmplx(0.,0.); rho_5 = cmplx(0.,0.)    ! set var to zero
                call CP_map((rho_i + uu*delta_rho), kr, sf, is, rho_4, Ppass)
                do kk = 1, size(km, 1)
                    rho_5 = rho_5 + matmul( matmul( km(kk,:,:), rho_4 ), km(kk,:,:) )
                enddo
                rho_5 = rho_5 / real(mat_trace(rho_5))
                f_2 = RelativeEntropy(rho_4, rho_5, siz)!22636
                if(f_2.le.f_1) then
                    tt = uu
                    f_1 = f_2
                endif
            enddo

        ! if f_1 == fr EXIT
            if(abs(f_1-fr) <= 1e-8) return
        ! assign new rho_i
            rho_i = rho_i + delta_rho * tt
        ! increment counter
            counter = counter + 1
        ! deallocation
            DEALLOCATE(delta_rho,rho_4,rho_5,x_i,c_i,F0_real,stat=iostat)
            call checkpoint(iostat==0,text="compute_primal::while loop deallocation failed")
        enddo
        DEALLOCATE(oj_real,rho_i,c_i,stat=iostat)
        return
    endsubroutine compute_primal

endmodule QKD





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