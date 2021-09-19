MODULE DUALSIMP_MOD
   ! Fortran 90 module for solving the asymetric dual of a problem
   !
   ! max C^T X
   ! s.t. A X \leq B
   !
   ! where A is dense and the dual solution is unique.
   !
   ! As far as the author knows, this code has tested robust in all situations
   ! where Y is unique. However, this code exists primarily for learning
   ! purposes, and has not been thoroughly tested for robustness. Neither does
   ! it represent the most efficient method for solving LPs of the above form
   ! in many situations.
   !
   ! Two subroutines are provided:
   !
   ! DUALSIMPLEX for solving an LP when the initial basis is known (Phase II of
   !    the simplex algorithm),
   !
   ! and
   !
   ! FEASIBLEBASIS for finding an initial dual feasible basis when none is
   !    known (Phase I of the simplex algorithm) by solving the auxiliary
   !    problem using DUALSIMPLEX.
   !
   ! Author: Tyler Chang
   ! Last Update: July, 2019
   
   ! Includes the R8 data type for approximately 64-bit precision.
   INTEGER, PARAMETER:: R8=SELECTED_REAL_KIND(13)
   
   CONTAINS
   
   SUBROUTINE DUALSIMPLEX (N, M, AT, B, C, IBASIS, X, Y, IERR, &
                         & EPS, IBUDGET, OBASIS)
   ! Solves a primal problem of the form
   !
   ! maximize     C^T X
   ! such that  A X \leq B
   !
   ! where A \in R^{M \times N}, C,X \in R^N, and B \in R^M.
   !
   ! Assuming M > N, this is done most efficiently by applying the simplex
   ! method on the asymmetric dual problem
   !
   ! minimize     B^T Y
   ! such that  A^T Y = C
   ! and        Y \geq 0
   !
   ! where Y \in R^M.
   !
   ! To solve the dual problem, the revised simplex algorithm is applied. In
   ! each iteration of the dual simplex algorithm, a complete LU factorization
   ! is performed via DGETRF. Dantzig's rule is used to pivot, until either a
   ! solution is found, or a pivot does not improve the objective. When a pivot
   ! fails to improve the objective, Bland's rule is used until improvement
   ! resumes, at which point Dantzig's pivoting strategy is resumed.
   !
   ! This strategy is most effective when A is dense and when M > N. For
   ! efficient memory access patterns, the constraints are specified by
   ! inputting A^T instead of A. For efficient linear algebra, LAPACK is
   ! used for computing LU factorizations and performing triangular solve
   ! operations, and BLAS is used for computing dot products and matrix-vector
   ! multiplication.
   !
   ! While a complete LU factorization in each iteration produces maximum
   ! numerical stability, a Woodbury update could alternatively be performed
   ! if iteration speed is of greater concern. When M >> N, the cost of the
   ! LU factorization is negligible compared to other computational costs,
   ! and the difference in speed would be negligible.
   !
   ! On input:
   !
   ! N is the integer number of variables in the primal problem.
   !
   ! M is the integer number of constraints in the primal problem.
   !
   ! AT(N,M) is the transpose of the real valued constraint matrix A.
   !
   ! B(M) is the real valued vector of upper bounds for the constraints.
   !
   ! C(N) is the real valued cost vector for the objective function.
   !
   ! IBASIS(N) is an integer valued vector containing the indices
   ! (from AT) of an intitial basis that is dual feasible.
   !
   ! On output:
   !
   ! X(N) is a real valued vector, which contains the primal solution.
   !
   ! Y(M) is a real valued vector, which contains the dual solution.
   !
   ! IERR is an integer valued error flag. The error codes are listed below:
   !
   ! Tens-digit is 0:
   !
   ! These codes indicate expected termination conditions.
   !
   !  0 : C^T X has been successfully maximized.
   !  1 : The dual problem is unbounded, therefore the primal must be infeasible.
   !
   ! Tens-digit is 1:
   !
   ! These codes indicate that the problem dimensions do not agree.
   !
   ! 10 : Illegal problem dimensions: N < 1.
   ! 11 : Illegal problem dimensions: M < N. If you wish to solve a problem with
   !      more variables than constraints, consider using a primal method.
   ! 12 : N does not match the first dimension of the constraint matrix AT.
   ! 13 : M does not match the second dimension of the constraint matrix AT.
   ! 14 : M does not match the length of the upper bounds B.
   ! 15 : N does not match the length of the cost vector C.
   ! 16 : N does not match the length of the initial basis IBASIS.
   ! 17 : N does not match the length of the primal solution vector X.
   ! 18 : M does not match the length of the dual solution vector Y.
   !
   ! Tens-digit is 2:
   !
   ! These codes indicate that the optional arguments contain illegal values
   ! or dimensions.
   !
   ! 20 : The optional argument EPS must be strictly positive.
   ! 21 : The optional argument IBUDGET must be nonnegative.
   ! 22 : The optional argument OBASIS must be length N.
   !
   ! Tens-digit is 3:
   !
   ! These codes indicate that the initial basis (IBASIS) was not feasible.
   !
   ! 30 : The provided initial basis IBASIS for AT contains indices that are
   !      out of the bounds of AT (either greater than M or less than 1).
   ! 31 : The provided initial basis IBASIS for AT contains duplicate indices,
   !      making it rank-deficient.
   ! 32 : The provided initial basis IBASIS for AT, while not redundant,
   !      produced a singularity
   ! 33 : The provided initial basis IBASIS for AT is not feasible.
   !
   ! Tens-digit is 4:
   !
   ! These codes indicate
   !
   ! 40 : The pivot budget (IBUDGET) was exceeded before a solution could be
   !      found. Consider increasing IBUDGET, or increasing the working
   !      precision EPS.
   ! 41 : A pivot has produced a singular basis. Consider increasing the
   !      working precision (EPS).
   ! 
   ! Tens-digit is 5:
   !
   ! These codes indicate a LAPACK error. If one of these errors occurs,
   ! it most likely indicates a system or compiler failure of some kind.
   !
   ! 50 : The subroutine DGETRF reported an illegal value (rare).
   ! 51 : The subroutine DGETRS reported an illegal value (rare).
   !
   ! Optional arguments:
   !
   ! EPS contains the working precision for the problem. EPS must be a
   ! strictly positive real number, and by default EPS is the square-root of
   ! the unit roundoff.
   !
   ! IBUDGET contains the integer budget for the maximum number of pivots
   ! allowed. By default, IBUDGET=50,000.
   !
   ! When present, OBASIS(:) returns the integer indices of the final basis
   ! as listed in AT.
   
   IMPLICIT NONE
   ! Parameter list.
   INTEGER, INTENT(IN) :: N
   INTEGER, INTENT(IN) :: M
   REAL(KIND=R8), INTENT(IN) :: AT(:,:)
   REAL(KIND=R8), INTENT(IN) :: B(:)
   REAL(KIND=R8), INTENT(IN) :: C(:)
   INTEGER, OPTIONAL, INTENT(INOUT) :: IBASIS(:)
   REAL(KIND=R8), INTENT(OUT) :: X(:)
   REAL(KIND=R8), INTENT(OUT) :: Y(:)
   INTEGER, INTENT(OUT) :: IERR
   REAL(KIND=R8), OPTIONAL, INTENT(IN) :: EPS
   INTEGER, OPTIONAL, INTENT(IN) :: IBUDGET
   INTEGER, OPTIONAL, INTENT(OUT) :: OBASIS(:)
   ! Local variables.
   INTEGER :: I, J, IBUDGETL
   REAL(KIND=R8) :: EPSL, NEWSOL, OLDSOL
   INTEGER :: IPIV(N), JPIV(M)
   REAL(KIND=R8) :: APIV(N,M), BPIV(M), S(M-N), LU(N,N)
   ! External procedures.
   REAL(KIND=R8), EXTERNAL :: DDOT  ! Inner product (BLAS).
   EXTERNAL :: DGEMV  ! General matrix vector multiply (BLAS)
   EXTERNAL :: DGETRF ! Perform a LU factorization with partial pivoting (LAPACK).
   EXTERNAL :: DGETRS ! Use the output of DGETRF to solve a linear system (LAPACK).
   
   ! Check inputs for errors.
   IF (N < 1) THEN
      IERR = 10; RETURN
   ELSE IF (M < N) THEN
      IERR = 11; RETURN
   ELSE IF (SIZE(AT,1) .NE. N) THEN
      IERR = 12; RETURN
   ELSE IF (SIZE(AT,2) .NE. M) THEN
      IERR = 13; RETURN
   ELSE IF (SIZE(B,1) .NE. M) THEN
      IERR = 14; RETURN
   ELSE IF (SIZE(C,1) .NE. N) THEN
      IERR = 15; RETURN
   ELSE IF (SIZE(IBASIS,1) .NE. N) THEN
      IERR = 16; RETURN
   ELSE IF (SIZE(X,1) .NE. N) THEN
      IERR = 17; RETURN
   ELSE IF (SIZE(Y,1) .NE. M) THEN
      IERR = 18; RETURN
   END IF
   ! Check for optionals.
   IF (PRESENT(EPS)) THEN
      IF(EPS .LE. 0.0_R8) THEN ! Must be strictly positive.
         IERR = 20; RETURN; END IF
      EPSL = EPS
   ELSE ! Set the default value.
      EPSL = EPSILON(0.0_R8)
   END IF
   IF (PRESENT(IBUDGET)) THEN
      IF(IBUDGET < 0) THEN ! Must be nonnegative.
         IERR = 21; RETURN; END IF
      IBUDGETL = IBUDGET
   ELSE ! Set the default value.
      IBUDGETL = 50000
   END IF
   IF (PRESENT(OBASIS)) THEN
      IF(SIZE(OBASIS, 1) .NE. N) THEN ! Must match the size of IBASIS.
         IERR = 22; RETURN; END IF
      OBASIS = 0 ! Initialize to zeros.
   END IF
   
   ! Check for issues with the initial basis, IBASIS.
   IF (ANY(IBASIS < 1) .OR. ANY(IBASIS > M)) THEN ! Check for illegal values.
      IERR = 30; RETURN; END IF
   DO I = 1, N ! Check for illegal bases.
      IF (ANY(IBASIS(I+1:N) .EQ. IBASIS(I))) THEN
         IERR = 31; RETURN; END IF
   END DO
   
   ! Initialize JPIV.
   FORALL ( I = 1 : M ) JPIV(I) = I
   ! Initilaize APIV and BPIV.
   APIV(:,:) = AT(:,:)
   BPIV(:) = B(:)
   ! Pivot the indices of APIV and BPIV to match IBASIS.
   DO I = 1, N
      ! Locate the current index of the basis element (after swapping).
      DO J = 1, M
         IF (JPIV(J) .EQ. IBASIS(I)) EXIT
      END DO
      ! Pivot APIV and BPIV to match the initial basis specified in IBASIS.
      CALL DSWAP(N, APIV(:,I), 1, APIV(:,J), 1)
      OLDSOL = BPIV(I)
      BPIV(I) = BPIV(J)
      BPIV(J) = OLDSOL
      ! Track the changes in JPIV.
      IPIV(1) = JPIV(I)
      JPIV(I) = JPIV(J)
      JPIV(J) = IPIV(1)
   END DO
   
   ! Get a solution using the LU factorization.
   LU(:,:) = APIV(:,1:N)
   CALL DGETRF(N, N, LU, N, IPIV, IERR)
   IF (IERR > 0) THEN ! LU is exactly singular.
      IERR = 32; RETURN
   ELSE IF (IERR < 0) THEN
      IERR = 50; RETURN; END IF
   ! Use the LU factorization to get the first N elements of the dual solution.
   Y(1:N) = C(:)
   CALL DGETRS('N', N, 1, LU, N, IPIV, Y, M, IERR)
   IF (IERR .NE. 0) THEN
      IERR = 51; RETURN; END IF
   IF (ANY(Y(1:N) .LT. -EPSL)) THEN
      IERR = 33; RETURN; END IF
   Y(N+1:M) = 0.0_R8 ! The last N elements are zeros.
   ! Given S(1:N)=0, use the LU factorization to get the primal solution.
   X(:) = BPIV(1:N)
   CALL DGETRS('T', N, 1, LU, N, IPIV, X, N, IERR)
   IF (IERR .NE. 0) THEN
      IERR = 51; RETURN; END IF
   ! Get the rest of the slack variables by solving B - A*X for the kernel.
   S(:) = BPIV(N+1:M)
   CALL DGEMV('T', N, M-N, -1.0_R8, APIV(:,N+1:M), N, X, 1, 1.0_R8, S, 1)
   ! Check if the KKT conditions have been satisfied.
   IF (ALL(S .GE. -EPSL) ) THEN
      IERR = 0
      ! Undo the pivots in Y, so that Y can be output.
      CALL DLAPMT( .FALSE., 1, M, Y, 1, JPIV )
      ! Store the final basis in OBASIS.
      IF(PRESENT(OBASIS)) OBASIS(:) = JPIV(1:N)
      RETURN
   END IF
   
   ! If not, compute the current solution and begin the iteration.
   NEWSOL = DDOT(N, BPIV, 1, Y, 1)
   OLDSOL = NEWSOL + 1.0_R8
   ! Loop until a solution is found.
   DO I = 1, IBUDGETL
      ! Choose the pivot rule based on the improvement.
      IF (OLDSOL - NEWSOL > EPSL) THEN
         CALL PIVOT_DANTZIG() ! Use Dantzig's rule when improvement is made.
         IF (IERR .NE. 0) RETURN
      ELSE
         CALL PIVOT_BLAND() ! Use Bland's rule when stalled.
         IF (IERR .NE. 0) RETURN
      END IF
      ! Get a new solution using the LU factorization.
      LU(:,:) = APIV(:,1:N)
      CALL DGETRF(N, N, LU, N, IPIV, IERR)
      IF (IERR > 0) THEN ! LU is exactly singular.
         IERR = 41; RETURN
      ELSE IF (IERR < 0) THEN
         IERR = 50; RETURN; END IF
      ! Use the LU factorization to get the first N elements of the dual solution.
      Y(1:N) = C(:)
      CALL DGETRS('N', N, 1, LU, N, IPIV, Y, M, IERR)
      IF (IERR .NE. 0) THEN
         IERR = 51; RETURN; END IF
      ! Given S(1:N)=0, use the LU factorization to get the primal solution.
      X(:) = BPIV(1:N)
      CALL DGETRS('T', N, 1, LU, N, IPIV, X, N, IERR)
      IF (IERR .NE. 0) THEN
         IERR = 51; RETURN; END IF
      ! Get the rest of the slack variables by solving B - A*X for the kernel.
      S(:) = BPIV(N+1:M)
      CALL DGEMV('T', N, M-N, -1.0_R8, APIV(:,N+1:M), N, X, 1, 1.0_R8, S, 1)
      ! Check if the KKT conditions are satisfied.
      IF (ALL(S .GE. -EPSL) ) THEN
         ! Undo the pivots in Y, so that Y can be output.
         CALL DLAPMT( .FALSE., 1, M, Y, 1, JPIV )
         ! Store the final basis in OBASIS.
         IF(PRESENT(OBASIS)) OBASIS(:) = JPIV(1:N)
         RETURN
      END IF
      ! Save the current solution and update the new solution.
      OLDSOL = NEWSOL
      NEWSOL = DDOT(N, BPIV, 1, Y, 1)
   END DO
   ! Budget expired.
   IERR = 40
   RETURN
   
   CONTAINS
   
   SUBROUTINE PIVOT_DANTZIG()
   ! Pivot using Dantzig's minimum ratio method for fast convergence.
   !
   ! On input, assume that APIV(:,:) contains the basis and kernel of AT,
   ! in that order. Also, assume that LU contains the LU factorization
   ! of the basis and IPIV contains the corresponding pivot indices.
   ! Also assume that Y contains a dual feasible solution and S contains
   ! the non-basic slack variables.
   !
   ! Given the above, compute the entering and exiting indices (IENTER and IEXIT)
   ! and update APIV(:,:) and BPIV(:) accordingly, tracking the pivots in JPIV.
   INTEGER :: IENTER, IEXIT
   REAL(KIND=R8) :: W(N), CURRMIN
   ! Compute the entering index.
   IENTER = MINLOC(S, 1) + N
   ! Build a weight vector for the entering vertex using the LU factorization.
   W(:) = APIV(:,IENTER)
   CALL DGETRS('N', N, 1, LU, N, IPIV, W, N, IERR)
   IF (IERR .NE. 0) THEN
      IERR = 51; RETURN; END IF
   ! Compute the weight ratios and choose the exiting index.
   CURRMIN = HUGE(0.0_R8)
   IEXIT = 0
   DO J = 1, N
      IF (W(J) .LT. EPSL) CYCLE
      W(J) = Y(J) / W(J)
      IF (W(J) < CURRMIN) THEN
         CURRMIN = W(J)
         IEXIT = J
      END IF
   END DO
   ! Check that an exiting index was found.
   IF (IEXIT .EQ. 0) THEN
      IERR = 1
      RETURN
   END IF
   ! Perform the pivot operation on both AT and B.
   CALL DSWAP(N, APIV(:,IEXIT), 1, APIV(:,IENTER), 1) ! Pivot AT using DSWAP.
   W(1) = BPIV(IEXIT) ! Pivot B using W(1) as a temp variable.
   BPIV(IEXIT) = BPIV(IENTER)
   BPIV(IENTER) = W(1)
   ! Record the pivot in JPIV(:).
   J = JPIV(IENTER)
   JPIV(IENTER) = JPIV(IEXIT)
   JPIV(IEXIT) = J
   RETURN
   END SUBROUTINE PIVOT_DANTZIG
   
   SUBROUTINE PIVOT_BLAND()
   ! Pivot using Bland's anticycling rule to guarantee convergence on degenerate
   ! point sets.
   !
   ! On input, assume that APIV(:,:) contains the basis and kernel of AT,
   ! in that order. Also, assume that LU contains the LU factorization
   ! of the basis and IPIV contains the corresponding pivot indices.
   ! Also assume that Y contains a dual feasible solution and S contains
   ! the non-basic slack variables.
   !
   ! Given the above, compute the entering and exiting indices (IENTER and IEXIT)
   ! and update APIV(:,:) and BPIV(:) accordingly, tracking the pivots in JPIV.
   INTEGER :: IENTER, IEXIT
   REAL(KIND=R8) :: W(N), CURRMIN
   ! Compute the entering index. It is the first negative entry in the kernel.
   DO J = 1, M-N
      IF(S(J) < -EPSL) THEN
         IENTER = J+N; EXIT; END IF
   END DO
   ! Build a weight vector for the entering vertex using the LU factorization.
   W(:) = APIV(:,IENTER)
   CALL DGETRS('N', N, 1, LU, N, IPIV, W, N, IERR)
   IF (IERR .NE. 0) THEN
      IERR = 51; RETURN; END IF
   ! Compute the weight ratios and choose the exiting index.
   CURRMIN = HUGE(0.0_R8)
   IEXIT = 0
   DO J = 1, N
      IF (W(J) .LT. EPSL) CYCLE
      W(J) = Y(J) / W(J)
      IF (W(J) - CURRMIN < -EPSL) THEN
         CURRMIN = W(J)
         IEXIT = J
      END IF
   END DO
   ! Check that an exiting index was found.
   IF (IEXIT .EQ. 0) THEN
      IERR = 1
      RETURN
   END IF
   ! Perform the pivot operation on both AT and B.
   CALL DSWAP(N, APIV(:,IEXIT), 1, APIV(:,IENTER), 1) ! Pivot AT using DSWAP.
   W(1) = BPIV(IEXIT) ! Pivot B using W(1) as a temp variable.
   BPIV(IEXIT) = BPIV(IENTER)
   BPIV(IENTER) = W(1)
   ! Record the pivot in JPIV(:).
   J = JPIV(IENTER)
   JPIV(IENTER) = JPIV(IEXIT)
   JPIV(IEXIT) = J
   RETURN
   END SUBROUTINE PIVOT_BLAND
   
   END SUBROUTINE DUALSIMPLEX
   
   SUBROUTINE FEASIBLEBASIS (N, M, AT, C, BASIS, IERR, EPS, IBUDGET)
   ! Implement the simplex method to find a dual feasible basis for a primal
   ! problem of the form
   !
   ! maximize     C^T X
   ! such that  A X \leq B
   !
   ! where A \in R^{M \times N}, X \in R^N, and B \in R^M.
   !
   ! The asymmetric dual problem is then of the form
   !
   ! minimize     B^T Y
   ! such that  A^T Y = C
   ! and        Y \geq 0
   !
   ! where Y \in R^M.
   !
   ! A basis V \in R^{N \times N} for A^T is dual feasible if V Y = C is
   ! solvable with Y \geq 0.
   !
   ! Find the indices of A corresponding to V by solving the auxiliary problem
   !
   ! minimize     SUM(Z)
   ! such that  A_AUX [Z^T Y^T]^T = C
   ! and        Z \geq 0, Y \geq 0.
   !
   ! where Z \in R^N and A_AUX = [ SIGN_N | A^T ] \in R^{N \times N+M}.
   ! Here, SIGN_N denotes the N-by-N identity matrix, except that when
   ! C(I) < 0 then SIGN_N(I,I) = -1.
   !
   ! Trivially, columns 1, ..., N of A_AUX provide an initial basis for this
   ! problem, with solution Z = ABS(B) and Y = 0. When SUM(Z) has been minimized,
   ! if SUM(Z) > 0, the problem is infeasible. Otherwise, if SUM(Z) = 0, then
   ! a feasible basis is given by the columns of A^T that are in the basis of
   ! A_AUX.
   !
   ! Uses DUALSIMPLEX to solve the auxiliary problem.
   !
   ! On input:
   !
   ! N is the integer number of variables in the primal problem.
   !
   ! M is the integer number of constraints in the primal problem.
   !
   ! AT(N,M) is the transpose of the real valued constraint matrix A.
   !
   ! C(N) is the real valued cost vector for the objective function.
   !
   ! On output:
   !
   ! BASIS(N) is an integer valued vector, which contains the indices of a
   ! dual feasible basis for AT.
   !
   ! IERR is an integer valued error flag. The error codes are listed below:
   !
   ! Tens-digit is 0:
   !
   ! These codes indicate expected termination conditions.
   !
   !  0 : C^T X has been successfully maximized.
   !  1 : The dual problem is unbounded, therefore the primal must be infeasible.
   !  2 : The dual problem is infeasible, the primal may be unbounded or
   !      infeasible.
   !
   ! Tens-digit is 1:
   !
   ! These codes indicate that the problem dimensions do not agree.
   !
   ! 10 : Illegal problem dimensions: N < 1.
   ! 11 : Illegal problem dimensions: M < N. If you wish to solve a problem with
   !      more variables than constraints, consider using a primal method.
   ! 12 : N does not match the first dimension of the constraint matrix AT.
   ! 13 : M does not match the second dimension of the constraint matrix AT.
   ! 15 : N does not match the length of the cost vector C.
   ! 16 : N does not match the length of the output BASIS.
   !
   ! Tens-digit is 2:
   !
   ! These codes indicate that the optional arguments contain illegal values
   ! or dimensions.
   !
   ! 20 : The optional argument EPS must be strictly positive.
   ! 21 : The optional argument IBUDGET must be nonnegative.
   !
   ! Tens-digit is 4:
   !
   ! These codes indicate
   !
   ! 40 : The pivot budget (IBUDGET) was exceeded before a solution could be
   !      found. Consider increasing IBUDGET, or increasing the working
   !      precision EPS.
   ! 41 : A pivot has produced a singular basis. Consider increasing the
   !      working precision (EPS).
   ! 
   ! Tens-digit is 5:
   !
   ! These codes indicate a LAPACK error. If one of these errors occurs,
   ! it most likely indicates a system or compiler failure of some kind.
   !
   ! 50 : The subroutine DGETRF reported an illegal value (rare).
   ! 51 : The subroutine DGETRS reported an illegal value (rare).
   !
   ! Optional arguments:
   !
   ! EPS contains the working precision for the problem. EPS must be a
   ! strictly positive real number, and by default EPS is the square-root of
   ! the unit roundoff.
   !
   ! IBUDGET contains the integer budget for the maximum number of pivots
   ! allowed. By default, IBUDGET=50,000.
   
   IMPLICIT NONE
   ! Parameter list.
   INTEGER, INTENT(IN) :: N
   INTEGER, INTENT(IN) :: M
   REAL(KIND=R8), INTENT(IN) :: AT(:,:)
   REAL(KIND=R8), INTENT(IN) :: C(:)
   INTEGER, INTENT(OUT) :: BASIS(:)
   INTEGER, INTENT(OUT) :: IERR
   REAL(KIND=R8), OPTIONAL, INTENT(IN) :: EPS
   INTEGER, OPTIONAL, INTENT(IN) :: IBUDGET
   ! Local variables.
   INTEGER :: I, IDX, IBUDGETL, IBASIS(N)
   REAL(KIND=R8) :: EPSL
   REAL(KIND=R8) :: A_AUX(N,N+M)
   REAL(KIND=R8) :: B_AUX(N+M)
   REAL(KIND=R8) :: Y_AUX(M+N)
   REAL(KIND=R8) :: X(N)
   ! External (BLAS) function DDOT for inner products.
   REAL(KIND=R8), EXTERNAL :: DDOT
   
   ! Check inputs for errors.
   IF (N < 1) THEN
      IERR = 10; RETURN
   ELSE IF (M < N) THEN
      IERR = 11; RETURN
   ELSE IF (SIZE(AT,1) .NE. N) THEN
      IERR = 12; RETURN
   ELSE IF (SIZE(AT,2) .NE. M) THEN
      IERR = 13; RETURN
   ELSE IF (SIZE(C,1) .NE. N) THEN
      IERR = 15; RETURN
   ELSE IF (SIZE(BASIS,1) .NE. N) THEN
      IERR = 16; RETURN
   END IF
   ! Check for optionals.
   IF (PRESENT(EPS)) THEN
      IF(EPS .LE. 0.0_R8) THEN ! Must be strictly positive.
         IERR = 20; RETURN; END IF
      EPSL = EPS
   ELSE ! Set the default value.
      EPSL = EPSILON(0.0_R8)
   END IF
   IF (PRESENT(IBUDGET)) THEN
      IF(IBUDGET < 0) THEN ! Must be nonnegative.
         IERR = 21; RETURN; END IF
      IBUDGETL = IBUDGET
   ELSE ! Set the default value.
      IBUDGETL = 50000
   END IF
   
   ! Copy AT into A_AUX to solve Phase 1 problem and zero B_AUX.
   A_AUX(:,N+1:N+M) = AT(:,:)
   B_AUX(1:N) = 1.0_R8
   B_AUX(N+1:N+M) = 0.0_R8
   ! Create artificial variables to solve the problem.
   A_AUX(:,1:N) = 0.0_R8
   FORALL ( I = 1 : N ) A_AUX(I,I) = SIGN(1.0_R8, C(I))
   ! Set IBASIS.
   FORALL (I = 1 : N) IBASIS(I) = I
   ! Get solution.
   CALL DUALSIMPLEX(N, N+M, A_AUX, B_AUX, C, IBASIS, X, Y_AUX, IERR, &
      & EPS=EPSL, IBUDGET=IBUDGETL, OBASIS=BASIS)
   IF (IERR .NE. 0) RETURN
   ! Check for infeasible dual solution.
   IF (DDOT(N, Y_AUX, 1, B_AUX, 1) > EPSL) THEN; IERR = 2; RETURN; END IF
   BASIS(:) = BASIS(:) - N
   ! Check that all basis elements are legal (in case of degeneracies).
   DO I = 1, N
      IF (BASIS(I) < 1) THEN
         IDX = 1 ! Find an available legal basis element.
         DO WHILE (ANY(BASIS(:) .EQ. IDX)); IDX = IDX + 1; END DO
         BASIS(I) = IDX
      END IF
   END DO
   RETURN
   
   END SUBROUTINE FEASIBLEBASIS
   

   subroutine nelmin ( fn, n, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
      icount, numres, ifault )
    
    !*****************************************************************************80
    !
    !! NELMIN minimizes a function using the Nelder-Mead algorithm.
    !
    !  Discussion:
    !
    !    This routine seeks the minimum value of a user-specified function.
    !
    !    Simplex function minimisation procedure due to Nelder and Mead (1965),
    !    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
    !    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
    !    25, 97) and Hill(1978, 27, 380-2)
    !
    !    The function to be minimized must be defined by a function of
    !    the form
    !
    !      function fn ( x, f )
    !      real ( kind = 8 ) fn
    !      real ( kind = 8 ) x(*)
    !
    !    and the name of this subroutine must be declared EXTERNAL in the
    !    calling routine and passed as the argument FN.
    !
    !    This routine does not include a termination test using the
    !    fitting of a quadratic surface.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    27 February 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by R ONeill.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    John Nelder, Roger Mead,
    !    A simplex method for function minimization,
    !    Computer Journal,
    !    Volume 7, 1965, pages 308-313.
    !
    !    R ONeill,
    !    Algorithm AS 47:
    !    Function Minimization Using a Simplex Procedure,
    !    Applied Statistics,
    !    Volume 20, Number 3, 1971, pages 338-345.
    !
    !  Parameters:
    !
    !    Input, external FN, the name of the function which evaluates
    !    the function to be minimized.
    !
    !    Input, integer ( kind = 4 ) N, the number of variables.
    !    0 < N is required.
    !
    !    Input/output, real ( kind = 8 ) START(N).  On input, a starting point
    !    for the iteration.  On output, this data may have been overwritten.
    !
    !    Output, real ( kind = 8 ) XMIN(N), the coordinates of the point which
    !    is estimated to minimize the function.
    !
    !    Output, real ( kind = 8 ) YNEWLO, the minimum value of the function.
    !
    !    Input, real ( kind = 8 ) REQMIN, the terminating limit for the variance
    !    of the function values.  0 < REQMIN is required.
    !
    !    Input, real ( kind = 8 ) STEP(N), determines the size and shape of the
    !    initial simplex.  The relative magnitudes of its elements should reflect
    !    the units of the variables.
    !
    !    Input, integer ( kind = 4 ) KONVGE, the convergence check is carried out
    !    every KONVGE iterations. 0 < KONVGE is required.
    !
    !    Input, integer ( kind = 4 ) KCOUNT, the maximum number of function
    !    evaluations.
    !
    !    Output, integer ( kind = 4 ) ICOUNT, the number of function evaluations
    !    used.
    !
    !    Output, integer ( kind = 4 ) NUMRES, the number of restarts.
    !
    !    Output, integer ( kind = 4 ) IFAULT, error indicator.
    !    0, no errors detected.
    !    1, REQMIN, N, or KONVGE has an illegal value.
    !    2, iteration terminated because KCOUNT was exceeded without convergence.
    !
      implicit none
    
      integer ( kind = 4 ) n
    
      real ( kind = 8 ), parameter :: ccoeff = 0.5D+00
      real ( kind = 8 ) del
      real ( kind = 8 ), parameter :: ecoeff = 2.0D+00
      real ( kind = 8 ), parameter :: eps = 0.001D+00
      real ( kind = 8 ), external :: fn
      integer ( kind = 4 ) i
      integer ( kind = 4 ) icount
      integer ( kind = 4 ) ifault
      integer ( kind = 4 ) ihi
      integer ( kind = 4 ) ilo
      integer ( kind = 4 ) j
      integer ( kind = 4 ) jcount
      integer ( kind = 4 ) kcount
      integer ( kind = 4 ) konvge
      integer ( kind = 4 ) l
      integer ( kind = 4 ) numres
      real ( kind = 8 ) p(n,n+1)
      real ( kind = 8 ) p2star(n)
      real ( kind = 8 ) pbar(n)
      real ( kind = 8 ) pstar(n)
      real ( kind = 8 ), parameter :: rcoeff = 1.0D+00
      real ( kind = 8 ) reqmin
      real ( kind = 8 ) rq
      real ( kind = 8 ) start(n)
      real ( kind = 8 ) step(n)
      real ( kind = 8 ) x
      real ( kind = 8 ) xmin(n)
      real ( kind = 8 ) y(n+1)
      real ( kind = 8 ) y2star
      real ( kind = 8 ) ylo
      real ( kind = 8 ) ynewlo
      real ( kind = 8 ) ystar
      real ( kind = 8 ) z
    !
    !  Check the input parameters.
    !
      if ( reqmin <= 0.0D+00 ) then
        ifault = 1
        return
      end if
    
      if ( n < 1 ) then
        ifault = 1
        return
      end if
    
      if ( konvge < 1 ) then
        ifault = 1
        return
      end if
    !
    !  Initialization.
    !
      icount = 0
      numres = 0
      jcount = konvge
      del = 1.0D+00
      rq = reqmin * real ( n, kind = 8 )
    !
    !  Initial or restarted loop.
    !
      do
    
        p(1:n,n+1) = start(1:n)
        y(n+1) = fn ( start )
        icount = icount + 1
    !
    !  Define the initial simplex.
    !
        do j = 1, n
          x = start(j)
          start(j) = start(j) + step(j) * del
          p(1:n,j) = start(1:n)
          y(j) = fn ( start )
          icount = icount + 1
          start(j) = x
        end do
    !
    !  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
    !  the vertex of the simplex to be replaced.
    !
        ilo = minloc ( y(1:n+1), 1 )
        ylo = y(ilo)
    !
    !  Inner loop.
    !
        do while ( icount < kcount )
    !
    !  YNEWLO is, of course, the HIGHEST value???
    !
          ihi = maxloc ( y(1:n+1), 1 )
          ynewlo = y(ihi)
    !
    !  Calculate PBAR, the centroid of the simplex vertices
    !  excepting the vertex with Y value YNEWLO.
    !
          do i = 1, n
            pbar(i) = ( sum ( p(i,1:n+1) ) - p(i,ihi) ) / real ( n, kind = 8 )
          end do
    !
    !  Reflection through the centroid.
    !
          pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
          ystar = fn ( pstar )
          icount = icount + 1
    !
    !  Successful reflection, so extension.
    !
          if ( ystar < ylo ) then
    
            p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
            y2star = fn ( p2star )
            icount = icount + 1
    !
    !  Retain extension or contraction.
    !
            if ( ystar < y2star ) then
              p(1:n,ihi) = pstar(1:n)
              y(ihi) = ystar
            else
              p(1:n,ihi) = p2star(1:n)
              y(ihi) = y2star
            end if
    !
    !  No extension.
    !
          else
    
            l = 0
            do i = 1, n + 1
              if ( ystar < y(i) ) then
                l = l + 1
              end if
            end do
    
            if ( 1 < l ) then
    
              p(1:n,ihi) = pstar(1:n)
              y(ihi) = ystar
    !
    !  Contraction on the Y(IHI) side of the centroid.
    !
            else if ( l == 0 ) then
    
              p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
              y2star = fn ( p2star )
              icount = icount + 1
    !
    !  Contract the whole simplex.
    !
              if ( y(ihi) < y2star ) then
    
                do j = 1, n + 1
                  p(1:n,j) = ( p(1:n,j) + p(1:n,ilo) ) * 0.5D+00
                  xmin(1:n) = p(1:n,j)
                  y(j) = fn ( xmin )
                  icount = icount + 1
                end do
    
                ilo = minloc ( y(1:n+1), 1 )
                ylo = y(ilo)
    
                cycle
    !
    !  Retain contraction.
    !
              else
                p(1:n,ihi) = p2star(1:n)
                y(ihi) = y2star
              end if
    !
    !  Contraction on the reflection side of the centroid.
    !
            else if ( l == 1 ) then
    
              p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
              y2star = fn ( p2star )
              icount = icount + 1
    !
    !  Retain reflection?
    !
              if ( y2star <= ystar ) then
                p(1:n,ihi) = p2star(1:n)
                y(ihi) = y2star
              else
                p(1:n,ihi) = pstar(1:n)
                y(ihi) = ystar
              end if
    
            end if
    
          end if
    !
    !  Check if YLO improved.
    !
          if ( y(ihi) < ylo ) then
            ylo = y(ihi)
            ilo = ihi
          end if
    
          jcount = jcount - 1
    
          if ( 0 < jcount ) then
            cycle
          end if
    !
    !  Check to see if minimum reached.
    !
          if ( icount <= kcount ) then
    
            jcount = konvge
    
            x = sum ( y(1:n+1) ) / real ( n + 1, kind = 8 )
            z = sum ( ( y(1:n+1) - x )**2 )
    
            if ( z <= rq ) then
              exit
            end if
    
          end if
    
        end do
    !
    !  Factorial tests to check that YNEWLO is a local minimum.
    !
        xmin(1:n) = p(1:n,ilo)
        ynewlo = y(ilo)
    
        if ( kcount < icount ) then
          ifault = 2
          exit
        end if
    
        ifault = 0
    
        do i = 1, n
          del = step(i) * eps
          xmin(i) = xmin(i) + del
          z = fn ( xmin )
          icount = icount + 1
          if ( z < ynewlo ) then
            ifault = 2
            exit
          end if
          xmin(i) = xmin(i) - del - del
          z = fn ( xmin )
          icount = icount + 1
          if ( z < ynewlo ) then
            ifault = 2
            exit
          end if
          xmin(i) = xmin(i) + del
        end do
    
        if ( ifault == 0 ) then
          exit
        end if
    !
    !  Restart the procedure.
    !
        start(1:n) = xmin(1:n)
        del = eps
        numres = numres + 1
    
      end do
    
      return
    end
    subroutine timestamp ( )
    
    !*****************************************************************************80
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    31 May 2001   9:45:54.872 AM
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 May 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
      implicit none
    
      character ( len = 8 ) ampm
      integer ( kind = 4 ) d
      integer ( kind = 4 ) h
      integer ( kind = 4 ) m
      integer ( kind = 4 ) mm
      character ( len = 9 ), parameter, dimension(12) :: month = (/ &
        'January  ', 'February ', 'March    ', 'April    ', &
        'May      ', 'June     ', 'July     ', 'August   ', &
        'September', 'October  ', 'November ', 'December ' /)
      integer ( kind = 4 ) n
      integer ( kind = 4 ) s
      integer ( kind = 4 ) values(8)
      integer ( kind = 4 ) y
    
      call date_and_time ( values = values )
    
      y = values(1)
      m = values(2)
      d = values(3)
      h = values(5)
      n = values(6)
      s = values(7)
      mm = values(8)
    
      if ( h < 12 ) then
        ampm = 'AM'
      else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h < 12 ) then
          ampm = 'PM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if
    
      write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
        d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
    
      return
    end
   END MODULE DUALSIMP_MOD