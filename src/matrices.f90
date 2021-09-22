module matrices
    use debugging
    implicit none
    
    interface mat_trace
        module procedure mat_trace_complex8, mat_trace_real8
    endinterface

    interface det
        module procedure det_real, det_complex
    endinterface

    contains

    real(8) function DET_real(aa)
        real(8) aa(:,:)
        real(8) tmp,c(size(aa,dim=1),size(aa,dim=2))
        real(8) max
        integer i,k,l,m,num(size(aa,dim=1)),n
        n=size(aa,dim=1)
        DET_real=1.    
        do k=1,n
            max=aa(k,k);num(k)=k;
            do i=k+1,n 
                if(abs(max)<abs(aa(i,k))) then
                    max=aa(i,k)
                    num(k)=i
                endif
            enddo
            if (num(k)/=k) then
                do l=k,n 
                    tmp=aa(k,l)
                    aa(k,l)=aa(num(k),l)
                    aa(num(k),l)=tmp
                enddo
                DET_real=-1.*DET_real
            endif
            do m=k+1,n
                c(m,k)=aa(m,k)/aa(k,k)
                do l=k,n 
                    aa(m,l)=aa(m,l)-c(m,k)*aa(k,l)
                enddo
            enddo !There we made matrix triangular!    
        enddo

        do i=1,n
            DET_real=DET_real*aa(i,i)
        enddo
        return
    endfunction DET_real

    real(8) function DET_complex(cc)
        complex(8) cc(:,:)
        real(8) rr(2*size(cc,dim=1),2*size(cc,dim=2))
        
        call complex_to_realm(size(cc,dim=1),cc,rr)

        DET_complex = 0.
        DET_complex = DET_real(rr)
        DET_complex = sqrt(DET_complex)
        return
    endfunction DET_complex

    function Kronecker_product(A,B)
    !Kronecker_product(
    !                    A,B(double complex)='input matrices'
    !                   )
    !returns the Kronecker product AxB
        implicit none
        !
        double complex, intent(in)        :: A(:,:), B(:,:)
        double complex, allocatable       :: Kronecker_product(:,:)
        integer*8                         :: ii, jj, kk, ll
        !
        allocate(Kronecker_product( size(A,1)*size(B,1) , size(A,2)*size(B,2) ))
        !
        do ii=1,size(A,1)
            do jj=1,size(A,2)
                do kk=1,size(B,1)
                    do ll=1,size(B,2)
                        Kronecker_product(kk + (ii-1)*size(B,1) , ll + (jj-1)*size(B,2)) = A(ii,jj) * B(kk,ll)
                    enddo
                enddo
            enddo
        enddo
        return
    endfunction Kronecker_product
    
    function Outer_product(A,B)
    !Outer_product(
    !                A, B(complex(8)) = 'array inputs'
    !              )
    !return the outer product of to arrays
        implicit none
        double complex, intent(in)  :: A(:), B(:)
        double complex, allocatable :: Outer_product(:,:)
        integer*8                   :: ii, jj
        allocate(Outer_product(size(A),size(B)))
        do ii=1,size(B)
            do jj=1,size(A)
                Outer_product(jj,ii) = A(jj) * B(ii)
            enddo
        enddo
        return
    endfunction Outer_product

    function Direct_sum(A,B) result(C)
    !Direct_sum(
    !                A, B(complex(8)) = 'array inputs'
    !              )
    !return the direct sum of to matrices
        complex(8), dimension (:,:), intent(in)     :: A, B
        complex(8), dimension (:,:), allocatable    :: C
        integer                                     :: p = 0, q = 0
        allocate(C(size(A,1)+size(B,1),size(A,2)+size(B,2)))
        C = 0
        p = size(A,1) + size(B,1) 
        q = size(A,2) + size(B,2) 
        C(1:size(A,1),1:size(A,2)) = A
        C(size(A,1)+1:p,size(A,2)+1:q) = B
        return
    end function Direct_sum

    subroutine eigensolver(mat,EIGVAL,EIGVEC,N)
    !eigensolver(
    !              A(double complex) = 'matrice hermitiana input che deve già essere allocata'
    !              W(double precision) = 'autovalori deve essere già allocato' 
    !              NV(char) = 'type of jonz, V or N'
    !              info(int) = 'info'
    !            )
    !    returns -> A (matrix of the eigenvectors), W (array with eigenvalues)
        use debugging
        implicit none
        !
        INTEGER, INTENT(IN)         :: N
        complex(8), INTENT(IN)      :: mat(N,N)
        real(8), INTENT(INOUT), ALLOCATABLE      :: eigval(:)
        complex(8), INTENT(INOUT), ALLOCATABLE   :: eigvec(:,:)
        INTEGER                     :: LDA
        INTEGER                     LWMAX
        INTEGER                     :: INFO, LWORK
        DOUBLE PRECISION,ALLOCATABLE:: RWORK( : )
        COMPLEX(8), ALLOCATABLE     :: WORK( : )
        EXTERNAL                    :: ZHEEV
        EXTERNAL                    :: PRINT_MATRIX, PRINT_RMATRIX
        INTRINSIC                   :: INT, MIN
        !
        LWMAX = N*1000
        !allochiamo la memoria per i vettori e le matrici
        !RWORK dimension should be at least MAX(1,3*N-2)
        allocate(RWORK( 3*N-2 ))
        allocate(WORK( LWMAX ))

        allocate(eigval(n), eigvec(n,n))
        eigvec = mat
        LDA = N
        LWORK = -1
        CALL ZHEEV( 'V', 'L', N, eigvec, LDA, eigval, WORK, LWORK, RWORK, INFO )
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
        !
        CALL ZHEEV( 'V', 'L', N, eigvec, LDA, eigval, WORK, LWORK, RWORK, INFO )
        !
        call checkpoint( INFO == 0, text= 'eigensolver:: (i < 0):  if INFO = -i, the i-th argument had an illegal value&
                                            &(i > 0):  if INFO = i, the algorithm failed to converge; i&
                                            &off-diagonal elements of an intermediate tridiagonal&
                                            &form did not converge to zero.',var = info)
        DEALLOCATE(RWORK, WORK)
    end subroutine eigensolver

    subroutine reduced_density_matrix(rho,rho_red,ind_k,sd)
    !reduced_density_matrix(
    !                        rho(dcplx) = 'density matrix in input'
    !                        rho_red(dcmplx) = 'reduced desity matrix in output'
    !                        ind_k(int) = 'the index k in rho_k=Tr _{not k} rho'
    !                        subdim(int) = 'local dimension, dimension of the exited reduced matrix'
    !                        )
    !changes the value of rho_red 
        implicit none
        double complex, intent(in)        :: rho(:,:)                ! dens matrix
        double complex, intent(inout)    :: rho_red(:,:)            ! reduce dens. matrix
        integer*8, intent(in)            :: ind_k, sd            ! rho_k
        integer*8                        :: ii, jj, kk, rrow, col
        do ii=1,size(rho_red,1)
            do jj=1,size(rho_red,2)
                rho_red(ii,jj) = dcmplx(0d0,0d0)
                do kk=1,sd
                    ! ρA[ii;jj]=∑dBiρAB[(ii−1)sd+kk;(jj−1)sd+kk]
                    rrow=1+(mod((ii-1),(sd)**(ind_k-1)))+ (kk-1)*(sd*(ind_k-1)) +&
                    & sd*(ii-1-(mod((ii-1),(sd)**(ind_k-1)))) + (kk-1)*(2-ind_k)
                    col=1+(mod((jj-1),(sd)**(ind_k-1)))+ (kk-1)*(sd*(ind_k-1)) +&
                    & sd*(jj-1-(mod((jj-1),(sd)**(ind_k-1)))) + (kk-1)*(2-ind_k)
                    !write(*,'("[i=",i0",j=",i0,"]","[r=",i0",c=",i0,"]")')ii,jj,rrow,col
                    rho_red(ii,jj)=rho_red(ii,jj)+rho(rrow,col)
                enddo
            enddo
        enddo
        ! tr=1
        call checkpoint(abs(abs(mat_trace(rho_red))-1.d0)<1E-4,&
                        text='reduced density matrix TRACE is NOT 1! ',&
                        var = abs(mat_trace(rho_red)))
        ! is hermitian?
        call checkpoint(is_hermitian(rho_red,size(rho_red,1)),text='reduced matrix is NOT hermitian')
        return
    end subroutine reduced_density_matrix

    function column(array)
    !column (
    !           array(complex(8)) = 'from horizontal array of length n'
    !       )
    !returns the same array written as a column of a nx1 matrix
        implicit none
        complex(8), INTENT(IN):: array(:)
        complex(8), ALLOCATABLE:: column(:,:)
        INTEGER ii
        allocate(column(size(array),1))
        do ii = 1, size(array)
            column(ii, 1) = array(ii)
        enddo
        return 
    endfunction column

    function row(mat)
    !row (
    !           mat(complex(8)) = 'vertical nx1 matrix, i.e. column vector'
    !       )
    !returns the same array written as a row, .i.e. rank-1 array
        implicit none
        complex(8), INTENT(IN)  :: mat(:,:)
        complex(8), allocatable :: row(:)
        INTEGER ii
        allocate(row(size(mat,1)))
        do ii = 1, size(mat,1)
            row(ii) = mat(ii,1)
        enddo
        return 
    endfunction row
    
    function invert_hermitian(m, n)
    !invert_hermitian
        implicit none
        integer, intent(in) :: n
        complex(8), intent(in) :: m(n,n)
        complex(8)  :: invert_hermitian(n,n)
        complex(8) :: work(n)
        integer :: ipiv(n), i, j, ierr
        
        invert_hermitian = m

        call zhetrf('U',n,invert_hermitian,n,ipiv,work,n,ierr)
        if (ierr.ne.0) stop 'Upper triangular factorization failed'
    
        call zhetri('U',n,invert_hermitian,n,ipiv,work,ierr)
        if (ierr.ne.0) stop 'Invert failed'
    
        do i=1,n
          do j=1,i-1
            invert_hermitian(i,j) = dconjg(invert_hermitian(j,i))
          enddo
        enddo
    endfunction invert_hermitian

    function sqrtM(matrix,n) result(a)
    !sqrtM (
    !           matrix(complex(8)) = 'input matrix'
    !           n(integer) = '# rows of matrix'
    !       )
    !returns the sqrt of the matrix
    !calculate the sqrt of a matrix using the Schur decomposition
        implicit none
        integer, INTENT(IN)     :: n
        complex(8), INTENT(IN)  :: matrix(n,n)
        Integer                 :: info, lda, ldc, ldd, ldvs, lwork, sdim, nb, ii
        complex(8), Allocatable :: a(:, :), c(:, :), d(:, :), vs(:, :), &
                                    w(:), work(:), sqr(:,:), vrinv(:,:), ipiv(:)
        complex(8)              :: wdum(1)
        Real(8), Allocatable    :: rwork(:)
        Logical                 :: dummy(1), select
        Intrinsic               :: cmplx, epsilon, max, nint, real
        nb = 64; lda = n; ldc = n; ldd = n; ldvs = n
        Allocate (a(lda,n), vs(ldvs,n), c(ldc,n), d(ldd,n), w(n), rwork(n))
    
        ! Use routine workspace query to get optimal workspace.
        lwork = -1
        Call zgees('V', 'N', select, n, a, lda, sdim, w, vs, ldvs, wdum, lwork, rwork, dummy, info)
    
        ! Make sure that there is enough workspace for block size nb.
        lwork = max((nb+1)*n, nint(real(wdum(1))))
        Allocate (work(lwork))
    
        ! Read in the matrix A
        a = matrix

        ! Find the Schur factorization of A
        Call zgees('V', 'N', select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, dummy, info)       
        !Check for convergence.
        call checkpoint(&
        &INFO==0,&
        &text="logM, zgees::&
        &(i<0): the i-th argument had an illegal value (program is terminated);\n&
        &(i=1 to n): the QR algorithm failed to compute the i eigenvalue;\n&
        &(i=n+1): The eigenvalues could not be reordered;\n&
        &(i=n+2): After reordering, roundoff changed values. [n,i]:",&
        &vec = [n, info]&
        &)

        allocate(sqr(n,n))
        sqr = cmplx(0.,0.)
        do ii = 1,n
            if(abs(w(ii)) > 1e-8 ) then
                sqr(ii, ii) = sqrt(W(ii))
            endif
        enddo

        ALLOCATE(vrinv(size(vs,1),size(vs,2)), ipiv(n))
        vrinv = vs
    
        call zgetrf( n, n, VRinv, n, ipiv, info )
        !Check for convergence.
        call checkpoint( INFO == 0, text= 'logM, zgetrf:: (i < 0):  if INFO = -i, the i-th argument had an illegal value&
                                            &(i > 0):  if INFO = i, the algorithm failed to converge; i&
                                            &off-diagonal elements of an intermediate tridiagonal&
                                            &form did not converge to zero.',var = info)    
        call zgetri( n, VRinv, n, ipiv, WORK, LWORK, INFO )
        !Check for convergence.
        call checkpoint( INFO == 0, text= 'logM, zgetri:: (i < 0):  if INFO = -i, the i-th argument had an illegal value&
                                            &(i > 0):  if INFO = i, the algorithm failed to converge; i&
                                            &off-diagonal elements of an intermediate tridiagonal&
                                            &form did not converge to zero.',var = info)        
        a = matmul(matmul(vs,sqr),VRinv)  
        deallocate(sqr, vs, c, vrinv, w, work, rwork, ipiv)
    endfunction sqrtM

    function logM(matrix,n,prec) result(A)
    !logM (
    !           mat(complex(8)) = 'input matrix'
    !       )
    !returns the log of the matrix
        implicit none
        integer, INTENT(IN)     :: n
        complex(8), INTENT(IN)  :: matrix(n,n)
        real(4), INTENT(IN), optional  :: prec
        Integer                 :: info, lda, ldc, ldd, ldvs, lwork, sdim, nb, ii
        complex(8), Allocatable :: a(:, :), c(:, :), d(:, :), vs(:, :), &
                                    w(:), work(:), sqr(:,:), vrinv(:,:), ipiv(:)
        complex(8)              :: wdum(1)
        Real(8), Allocatable    :: rwork(:)
        real(8)                 :: fudge
        Logical                 :: dummy(1), select
        Intrinsic               :: cmplx, epsilon, max, nint, real
        nb = 64; lda = n; ldc = n; ldd = n; ldvs = n
        Allocate (a(lda,n), vs(ldvs,n), c(ldc,n), d(ldd,n), w(n), rwork(n))

        ! precision
        if(present(prec))then
            fudge = prec
        else
            fudge = 1e-10
        endif

        ! Read in the matrix A
        a = (1-fudge)*matrix+fudge*identity(int8(n))

        ! Use routine workspace query to get optimal workspace.
        lwork = -1
        Call zgees('V', 'N', select, n, a, lda, sdim, w, vs, ldvs, wdum, lwork, rwork, dummy, info)
    
        ! Make sure that there is enough workspace for block size nb.
        lwork = max((nb+1)*n, nint(real(wdum(1))))
        Allocate (work(lwork))
    
        info = 1        
        do while(info/=0)
            ! Read in the matrix A
            a = (1.-fudge)*matrix + fudge*identity(int8(n))


            ! Find the Schur factorization of A
            Call zgees('V', 'N', select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, dummy, info)       

            if (info/=0)then
                fudge = fudge*1e5
                write(*,'(A,es20.10)') "logM:: precision set to ",fudge
            endif
        enddo

        allocate(sqr(n,n))
        sqr = cmplx(0.,0.)
        do ii = 1,n
            if(abs(w(ii)) > 1e-8 ) then
                sqr(ii, ii) = zlog(W(ii))
            endif
        enddo

        ALLOCATE(vrinv(size(vs,1),size(vs,2)), ipiv(n))
        vrinv = vs
    
        call zgetrf( n, n, VRinv, n, ipiv, info )
        !Check for convergence.
        call checkpoint( INFO == 0, text= 'logM, zgetrf:: (i < 0):  if INFO = -i, the i-th argument had an illegal value&
                                            &(i > 0):  if INFO = i, the algorithm failed to converge; i&
                                            &off-diagonal elements of an intermediate tridiagonal&
                                            &form did not converge to zero.',var = info)    
        call zgetri( n, VRinv, n, ipiv, WORK, LWORK, INFO )
        !Check for convergence.
        call checkpoint( INFO == 0, text= 'logM, zgetri:: (i < 0):  if INFO = -i, the i-th argument had an illegal value&
                                            &(i > 0):  if INFO = i, the algorithm failed to converge; i&
                                            &off-diagonal elements of an intermediate tridiagonal&
                                            &form did not converge to zero.',var = info)        
        a = matmul(matmul(vs,sqr),VRinv)  
        deallocate(sqr, vs, c, vrinv, w, work, rwork, ipiv)
    endfunction logM

    function logMV(matrix,n,prec) result(A)
    !logMV (
    !           mat(complex(8)) = 'input matrix'
    !       )
    !returns the log of the matrix
        implicit none
        integer,INTENT(IN)     :: n
        complex(8),INTENT(IN)  :: matrix(n,n)
        real(4), INTENT(IN), optional::prec
        Integer                :: info, lwork, nb, ii
        parameter               (nb=64)
        complex(8),Allocatable :: a(:, :), w(:), work(:), diagm(:,:),&
        & vrinv(:,:), ipiv(:), vl(:,:),vr(:,:)
        Real(8),Allocatable    :: rwork(:)
        real(8)                 fudge

        Allocate (a(n,n), w(n), rwork(2*n))
        Allocate (vl(n,n), vr(n,n))
        lwork = (nb+1)*n!, nint(real(wdum(1))))
        Allocate (work(lwork))

        ! precision
        if(present(prec))then
            fudge=prec
        else
            fudge=1e-16
        endif

        info = 1        
        do while(info/=0)
            ! Read in the matrix A
            a = (1.-fudge)*matrix + fudge*identity(int8(n))

            call zgeev('N','V',n,a,n,w,VL,n,VR,n,WORK,LWORK,rwork,INFO)

            if (info/=0)then
                fudge = fudge*10
                write(*,'(A,es20.10)') "logMV:: precision set to ",fudge
            endif
        enddo

        allocate(diagm(n,n))
        diagm = cmplx(0.,0.)
        do ii = 1,n
            if(abs(w(ii)) > 1e-10 ) then
                diagm(ii, ii) = zlog(W(ii))
            endif
        enddo

        ALLOCATE(vrinv(n,n), ipiv(n))
        vrinv = VR

        call zgetrf( n, n, VRinv, n, ipiv, info )
        !Check for convergence.
        call checkpoint(&
        &INFO==0,text='logMV, zgetrf:: (i < 0):  if INFO = -i, the i-th argument had an illegal value&
        &(i > 0):  if INFO = i, the algorithm failed to converge; i&
        &off-diagonal elements of an intermediate tridiagonal&
        &form did not converge to zero.',var=info)

        call zgetri( n, VRinv, n, ipiv, WORK, LWORK, INFO )
        !Check for convergence.
        call checkpoint(&
        &INFO==0,text='logMV, zgetri:: (i < 0):  if INFO = -i, the i-th argument had an illegal value&
        &(i > 0):  if INFO = i, the algorithm failed to converge; i&
        &off-diagonal elements of an intermediate tridiagonal&
        &form did not converge to zero.',var=info)

        a = matmul(matmul(VR,diagm),VRinv) 
        deallocate(diagm, vrinv, w, work, rwork, ipiv, vl, vr)
    endfunction logMV

    function logMZ(matrix,n) result(A)
        !logMZ (
        !           mat(complex(8)) = 'input matrix'
        !       )
        !returns the log of the matrix
            implicit none
            integer,INTENT(IN)     :: n
            complex(8),INTENT(IN)  :: matrix(n,n)
            Integer                :: info, lwork, nb, ii
            complex(8),Allocatable :: a(:,:), w(:), work(:), diagm(:,:),&
            & vrinv(:,:), ipiv(:), vl(:,:),vr(:,:)
            Real(8),Allocatable    :: rwork(:)
    
            Allocate (a(n,n), w(n), rwork(2*n))
            Allocate (vl(n,n), vr(n,n))
    
            ! Read in the matrix A
            vr = matrix
            ! Make sure that there is enough workspace for block size nb.
            nb = 64
            lwork = (nb+1)*n!, nint(real(wdum(1))))
            Allocate (work(lwork))
            call zheev('v','u',n,vr,n,w,WORK,LWORK,rwork,INFO)
    
            !Check for convergence.
            call checkpoint(INFO==0,text="logMZ, ZHEEV::&
            &(i<0): the i-th argument had an illegal value (program is terminated);\n&
            &(i=1 to n): the QR algorithm failed to compute the i eigenvalue;\n&
            &(i=n+1): The eigenvalues could not be reordered;\n&
            &(i=n+2): After reordering, roundoff changed values. [n,i]:",&
            &vec = [n, info]&
            &)
    
            allocate(diagm(n,n))
            diagm = cmplx(0.,0.)
            do ii = 1,n
                if(abs(w(ii)) > 1e-10 ) then
                    diagm(ii, ii) = zlog(W(ii))
                endif
            enddo
    
            ALLOCATE(vrinv(n,n), ipiv(n))
            vrinv = VR
    
            call zgetrf( n, n, VRinv, n, ipiv, info )
            !Check for convergence.
            call checkpoint(&
            &INFO==0,text='logMZ, zgetrf:: (i < 0):  if INFO = -i, the i-th argument had an illegal value&
            &(i > 0):  if INFO = i, the algorithm failed to converge; i&
            &off-diagonal elements of an intermediate tridiagonal&
            &form did not converge to zero.',var=info)
    
            call zgetri( n, VRinv, n, ipiv, WORK, LWORK, INFO )
            !Check for convergence.
            call checkpoint(&
            &INFO==0,text='logMZ, zgetri:: (i < 0):  if INFO = -i, the i-th argument had an illegal value&
            &(i > 0):  if INFO = i, the algorithm failed to converge; i&
            &off-diagonal elements of an intermediate tridiagonal&
            &form did not converge to zero.',var=info)
    
            a = matmul(matmul(VR,diagm),VRinv) 
            deallocate(diagm, vrinv, w, work, rwork, ipiv, vl, vr)
        endfunction logMZ

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

    subroutine Gram_Schmidt_decomposition(A,Q,R,m,n)
    ! Gram_Schmidt_decomposition 
        implicit none
        integer, INTENT(IN)         :: m,n
        complex(8), INTENT(IN)      :: A(m,n)
        complex(8), INTENT(INOUT)   :: Q(m,m), R(m,n)
        INTEGER                     :: LDA, LWORK, info, ii, jj
        complex(8), ALLOCATABLE     :: mat(:,:), TAU(:), WORK(:), id(:,:), v(:,:)
        allocate(mat(m,n), TAU(n), WORK(1))
        LWORK = -1; LDA = max(1,m)
        call zgeqrf(m, n, mat, LDA, TAU, WORK, LWORK, info)
  
        ! Make sure that there is enough workspace for block size nb.
        LWORK = MAX(1, int(work(1),kind=4))
        deallocate(work)
        Allocate (work(LWORK))
        call zgeqrf(m, n, mat, LDA, TAU, WORK, LWORK, info)

        !R
        DO ii = 1, n
            do jj = 1, m
                if (jj <= ii)then
                    R(jj,ii) = A(jj,ii)
                endif
            enddo
        ENDDO

        !Q: Q = H(1) H(2) . . . H(k), where k = min(m,n).
        allocate(id(m,m))
        id = identity(int(m,kind=8))
        allocate(v(m,1))
        do jj = 1, min(m,n)
            v = cmplx(0.,0.); v(jj,1) = cmplx(1.,0.); v(jj+1:m,1) = mat(jj+1:m,jj) 
            Q = Q * (id - TAU(jj)* matmul(v,transpose(v)))
        ENDDO
        DEALLOCATE(mat,tau, work, id, v)
    endsubroutine Gram_Schmidt_decomposition

    subroutine gram_schmidtm(A,n,k,l)
    !gram_schmidt for matrices(
    !               A(complex*8, intent(in)) = 'input matrix'
    !               Q(complex*8, intent(in)) = 'output decomposed matrix'
    !               R(complex*8, intent(in)) = 'output upper triangular matrix'
    !){second algorithm}
        implicit none
        INTEGER, INTENT(IN)             :: n, k, l !dimensions of A
        complex(8), INTENT(INOUT)       :: A(n,k,l)
        complex(8), ALLOCATABLE         :: R(:,:,:)
        complex(8)                      norm
        INTEGER                         jj, ii
        allocate(R(n,k,l))
        R = cmplx(0.,0.)
        do jj = 1, n
            R(jj,:,:) = A(jj,:,:)
            do ii = 1, (jj - 1)
                R(jj,:,:) = R(jj,:,:) - mat_trace( matmul(conjg(transpose( R(ii,:,:) )), R(jj,:,:) ) )/&
                                       &mat_trace( matmul(conjg(transpose( R(ii,:,:) )), R(ii,:,:) ) ) * R(ii,:,:)
            enddo
            norm = real(mat_trace( matmul(conjg(transpose(R(jj,:,:))), R(jj,:,:)) ))
            if ( abs(norm) >= 1e-8) then
                R(jj,:,:) =  R(jj,:,:) / norm
            else
                !delete the jj element: this means that the set of matrices is NOT orthogonal
                R = R(:jj-1,:,:)
                write(*,*) "gram_schmidtm:: set of matrices is NOT orthogonal. Exited at", jj-1
                exit
            endif
        enddo
        A = R
        return
    endsubroutine gram_schmidtm

    subroutine extend_basism(V, B, j, n)
    !Extend a set of j.leq.k orthogonal matrices and k.leq.m orthonormal matrices to basis of the space H^m.
    ! ---------------------------------------------
    ! Keyword arguments:
    ! V -- the new set of matrices (size (j, n, n) )
    ! B -- the orthonormal set to be extended (size (k, n, n))
    ! return C = (B,V^{orth})
        integer, INTENT(IN)                     :: j, n ! j=len(V), k=len(B), n=dim of the space 
        complex(8), intent(in)                  :: V(j, n, n)
        complex(8), intent(inout), ALLOCATABLE  :: B(:,:,:)
        complex(8), ALLOCATABLE                 :: U(:,:), C(:,:,:)
        INTEGER                                 ii, jj, kk, sz
        sz = size(B, 1)
        allocate(U(n,n))
        do ii = 1, j ! remember the Hilbert-schmidt norm of two operators A, B is NORM = Tr[ A @ B ]
            U = V(ii, :, :)
            do jj = 1, sz
                U = U - (mat_trace( matmul(transpose(conjg(B(jj,:,:))), U)) / &
                & mat_trace( matmul(transpose(conjg(B(jj,:,:))), B(jj,:,:) )) ) * B(jj,:,:)
            enddo
            if( abs(mat_trace( matmul( transpose(conjg(U)) , U) )) >= 1e-8 ) then
                allocate(C(sz,size(B,2),size(B,3)))
                C = B
                DEALLOCATE(B)
                sz = sz + 1
                ALLOCATE(B(sz,size(C,2),size(C,3)))
                do kk = 1, sz -1         
                    B(kk,:,:) = C(kk,:,:)
                enddo
                B(sz,:,:) = U
                DEALLOCATE(C)
            endif
        ENDDO
        DEALLOCATE(U)
        return
    endsubroutine extend_basism

    function identity(N)
    !identity(
    !            N(int)='number of elements on the diagonal'
    !          )
    !returns a Identity NxN matrix
        implicit none
            integer*8, intent(in)    :: N
            integer*8                    :: ii
            double complex            :: identity(N,N)
            identity=dcmplx(0.d0,0.d0)
            do ii=1,N
                identity(ii,ii)=dcmplx(1.d0,0.d0)
            enddo
        return
    endfunction identity

    function sigma(N)
    !sigma(
    !            N(char)='X,Y or Z index of pauli 2x2 matrix'
    !          )
    !returns a sigma_n 2x2 matrix    
        implicit none
        character(*), intent(in)    :: N
        double complex                :: sigma_x(2,2), sigma_y(2,2), sigma_z(2,2), sigma(2,2)
        sigma_x(1,:) = [dcmplx(0.d0,0.d0),dcmplx(1.d0,0.d0)]
        sigma_x(2,:) = [dcmplx(1.d0,0.d0),dcmplx(0.d0,0.d0)]
        sigma_y(1,:) = [dcmplx(0.d0,0.d0),dcmplx(0.d0,-1.d0)]            
        sigma_y(2,:) = [dcmplx(0.d0,1.d0),dcmplx(0.d0,0.d0)]
        sigma_z(1,:) = [dcmplx(1.d0,0.d0),dcmplx(0.d0,0.d0)]
        sigma_z(2,:) = [dcmplx(0.d0,0.d0),dcmplx(-1.d0,0.d0)]
        select case(N)
            case("x")
                sigma = sigma_x
            case("y")
                sigma = sigma_y
            case("z")
                sigma = sigma_z
            case default
                write(*,'("ERROR uncorrect input N ",A," must be X,Y, or Z")')N
                stop
        endselect
        return
    endfunction

    function is_hermitian(m,n)
    ! is_hermitian(
    !                mat(double complex) = 'input matrix'
    !                dump(int,optional) = '=1 print, else NOT'
    !                tol(real,optional) = ''
    !              )
    !return true if matrix is hermitian or false if not.
        implicit none
        integer, intent(in) :: n
        complex(8), intent(in) :: m(n,n)
        logical :: is_hermitian
        integer(8) ii,jj
    
        is_hermitian=.True.
        !is_hermitian = all(abs(mt-m) .le. 1e-8)
        do ii=1,n
            do jj=1,n
                if(abs(m(ii,jj)-conjg(m(jj,ii)))>1e-8)then

                    is_hermitian=.False.
                   
                    exit
                endif
            enddo
        enddo
    endfunction is_hermitian

    function is_positive(m,n,tol) result(CC)
    ! is_positive(
    !                mat(double complex) = 'input matrix'
    !                dump(int,optional) = '=1 print, else NOT'
    !                tol(real,optional) = ''
    !              )
    !return true if matrix is hermitian or false if not.
        implicit none
        integer, intent(in)         :: n
        real, intent(in), optional:: tol
        complex(8), intent(in)  :: m(n,n)
        complex(8), ALLOCATABLE :: mp(:,:)
        complex(8), ALLOCATABLE :: egv(:,:)
        real(8), ALLOCATABLE    :: egl(:)
        real                    tolerance
        integer                 :: ii
        logical                 :: CC

        ! Duplicate
        mp = m
        if(present(tol).eqv..true.)then;tolerance=tol;else;tolerance=1e-8;endif
        call eigensolver(mp,egl,egv,n)
        CC = .true.
        do ii = 1, n
            if (egl(ii) < - tolerance)then
                CC = .false.
                write(*,'(A,i0,A,f10.5)')"is_positive:: ",ii," eigenvalue negative ",egl(ii)
                RETURN
            endif
        enddo
        DEALLOCATE(mp, egl,egv)
        return
    endfunction is_positive

    function mat_trace_complex8(mat)
    !mat_trace_complex8(
    !            mat='array 1-dimensionale che contiene la lista degli elementi della matrice',
    !              )
    !la traccia è definita come la somma di tutti gli elementi di matrice presenti sulla diagonale della matrice mat
        implicit none
        !
        double complex, intent(in)    :: mat(:,:)
         double complex                :: mat_trace_complex8
         integer*4                    :: ii
        !
        mat_trace_complex8=cmplx(0.0,0.0)
        do ii=1,size(mat,1)
            mat_trace_complex8=mat_trace_complex8+mat(ii,ii)
        enddo
        !
        return
    endfunction mat_trace_complex8

    function mat_trace_real8(mat)
    !mat_trace_real8(
    !            mat='array 1-dimensionale che contiene la lista degli elementi della matrice',
    !              )
    !la traccia è definita come la somma di tutti gli elementi di matrice presenti sulla diagonale della matrice mat
            implicit none
            !
            real(8), intent(in)    :: mat(:,:)
            real(8)                :: mat_trace_real8
            integer*4              :: ii
            !
            mat_trace_real8=0_8
            do ii=1,size(mat,1)
                mat_trace_real8=mat_trace_real8+mat(ii,ii)
            enddo
            return
    endfunction mat_trace_real8
    
    function cnorm(psi, dx)
    !cnorm(
    !             psi(double complex) = 'array 1 dimensionale',
    !             dx(integer,opt) = 'larghezza step, default = 1'
    !           )
    ! return real norm for a complex array
        implicit none
        !
        double complex,  INTENT(INOUT)            :: psi(:)
        real*8, INTENT(IN),    optional            :: dx
        real*8                                    :: cnorm, step=1.d0
        integer*4                                :: ii
        !
        if(present(dx))step=dx
        !
        cnorm=0.d0
        do ii=1,size(psi)
           cnorm=cnorm+(real(psi(ii))**(2)+aimag(psi(ii))**(2))*step
        enddo
        cnorm=sqrt(cnorm)
        !
        return
    endfunction cnorm

    real function infinity()
        implicit none
        real :: varx
        varx = huge(1.)
        infinity = varx + varx
    end function infinity
endmodule matrices
