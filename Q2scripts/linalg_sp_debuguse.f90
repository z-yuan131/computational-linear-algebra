! ============================================================================ !
!                          *** linalg_sp.f90 ***                               !
! ============================================================================ !
! A library of FORTRAN subroutines for linear algebra in SINGLE PRECISION.
!
! INSTRUCTIONS
! -------------
! Complete this library by implementing the two missing subroutines:
!
! *  SP_LU_PP: LU factorization with partial pivoting in single precision
! *  SP_LU_CP: LU factorization with complete pivoting in single precision
!
! Include suitable comments, explaining CLEARLY how the subroutines should be
! used, their inputs and their outputs. The list of contents below gives a brief
! description of the three subroutines that this library should contain.
!
! CONTENTS
! --------
! subroutine SP_LU(n, A)
!       Subroutine to compute the LU factorization of a matrix A "in place".
!       Store the strictly lower triangular part of L in place of the strictly
!       lower triangular part of A, and the upper triangular factor U in place
!       of the upper triangular part of A. No pivoting is implemented.
!       *** THE CODE FOR THIS SUBROUTINE IS GIVEN AND HAS BEEN TESTED ***
!
! ! subroutine SP_LU_PP(n, A)
!       Subroutine to compute the LU factorization of a matrix A "in place"
!       With PARTIAL pivoting. This subroutine should store the strictly lower
!       triangular part of L in place of the strictly lower triangular part of
!       A, and the upper triangular factor U in place of the upper triangular
!       part of A. Also return a permutation vector p such that A(p,:) = L*U
!       *** YOU NEED TO WRITE THE CODE FOR THIS SUBROUTINE ***
!
! subroutine SP_LU_CP(n, A)
!       Subroutine to compute the LU factorization of a matrix A "in place"
!       With COMPLETE pivoting. This subroutine should store the strictly lower
!       triangular part of L in place of the strictly lower triangular part of
!       A, and the upper triangular factor U in place of the upper triangular
!       part of A. Also return two permutation vectors, p and q, such that
!       A(p,q) = L*U
!       *** YOU NEED TO WRITE THE CODE FOR THIS SUBROUTINE ***
!
! subroutine SP_UNIT_LSOLVE()
!       Subroutine to solve a lower triangular system Lx=b with L unit lower
!       triangular. Only the strictly lower triangular part of L is used since
!       the diagonal entries are known to be 1. To be used in conjunction with
!       the subroutine LU() to solve a general linear system of equations.
!       IMPORTANT: All real inputs and outputs should be SINGLE PRECISION arrays
!       *** THE CODE FOR THIS SUBROUTINE IS GIVEN AND HAS BEEN TESTED ***
!
! subroutine SP_USOLVE()
!       Subroutine to solve an upper triangular system Ux=b with upper
!       triangular coefficient matrix U. Only the upper triangular part of U is
!       used. To be used in conjunction with the subroutine LU() to solve a
!       general linear system of equations.
!       IMPORTANT: All real inputs and outputs should be SINGLE PRECISION arrays
!       *** THE CODE FOR THIS SUBROUTINE IS GIVEN AND HAS BEEN TESTED ***
!
! ============================================================================ !

program LU_zhenyangyuan
!
  implicit none 
  integer,parameter :: n=80
  real(4) :: A(n,n), p(1,n) , L(n,n) , U(n,n) , b(n,1) , b_temp(n,1) , q(1,n) , residual , temp(n,1) , b_res(n,1) , A_res(n,n)
  integer :: k, j , i
  

	
    do i = 1, n
		p(1 , i) = i;
		q(1 , i) = i;
	    b(i , 1) = sqrt(i+0.0);      !!!SQRT(X) where X must be real. i is an integer. So i + 0.0 is real.
	enddo
	b_res = b;
	
	do i = 1, n
	A(i , i) = 1;
	A(: , n) = 1;
		if (i > 1) then
			A(i , 1 : i - 1) = -1;
			A(i - 1 , i : n - 1) = 0;
		end if
	enddo
	A_res = A;
	
	
	!print *,'Before pivoting, A = '
	!do i = 1,n
	!		print *, A(i,:);
	!enddo
	!print *, A;
	!print *, p;
	
	!call SP_LU(n , A)
	call SP_LU_PP(n , p , A)
	!call SP_LU_CP(n , p , q , A)
	
	!do i = 1 , n
	!	do j = 1 , n
	!		if (i > j) then
	!			L(i , j) = A(i , j);
	!			U(i , j) = 0;
	!		else if (j >= i) then
	!			L(i , j) = 0;
	!			U(i , j) = A(i , j);
	!		end if
	!	enddo
	!	L(i , i) = 1;
	!enddo;
	!print *,'L = '
	!do i = 1,n
	!		print *, L(i,:);
	!enddo
	!print *,'U = '
	!do i = 1,n
	!		print *, U(i,:);
	!enddo
	
	
	
	b_temp = b;
	do i = 1 , n
		 b(p(1,i) , 1) = b_temp(i , 1);
	enddo
	
	
	!print *,'After pivoting, A = '
	!do i = 1,n
	!		print *, A(i,:);
	!enddo
	!print *,'premutation vector for column, p = '
	!print *, p;
	!print *, L;
	!print *, U;
	!print *, b;
	!print *,'premutation vector for row, q = '
	!print *, q;
	
	
	call SP_UNIT_LSOLVE(n , A , b)
	call SP_USOLVE(n , A , b)
	
	b_temp = b;
	do i = 1 , n
		 b(q(1,i) , 1) = b_temp(i , 1);
	enddo
	print *,'b = '
	print *, b
	
	
	!!!! residual norm
	!do i = 1 , n
			temp = abs(MATMUL(A_res,b) - b_res);
	!enddo
	residual = MAXVAL(temp(:,1))/MAXVAL(abs(b_res(:,1)));
	print *,'residual = '
	print *, residual;
	
end program LU_zhenyangyuan
				





subroutine SP_LU(n, A)
    ! "In place" LU factorization of an n-by-n single precision matrix A. No
    ! pivoting is implemented.
    !
    ! Inputs:
    ! -------
    ! integer :: n - the size of A
    ! real    :: A - the matrix A, datatype real (single precision)
    !
    ! Outputs:
    ! --------
    ! real :: A - a real n-by-n array, overwritten onto the input matrix A,
    !             containing the LU factors of A. Precisely:
    !             (1) The strictly lower triangular part of L is stored in the
    !                 strictly lower triangular part of A
    !             (2) The upper triangular factor U is stored in the upper
    !                 triangular part of A.

    ! Declare all variables, no implicits
    implicit none
    integer, intent(in) :: n
    real(4) :: A(n,n)
    integer :: k, j

    ! Perform Gaussian elimination "in place"
    do k = 1, n-1
        if (A(k,k).NE.0.0) then
            ! Compute all Gauss weights
            A(k+1:n,k) = A(k+1:n,k)/A(k,k)
            ! Perform gaussian elimination. Proceed by columns using array
            ! operations supported by FORTRAN for simplicity.
            do j = k+1, n
                A(k+1:n,j) = A(k+1:n,j) - A(k+1:n,k)*A(k,j)
            end do
        else
            ! A does not admit an LU factorization without pivoting
            ! Warn the user and exit
            print *, "A does not admit an LU factorization without pivoting"
            exit
        end if
    end do

! END SUBROUTINE
end subroutine SP_LU

! ---------------------------------------------------------------------------- !

subroutine SP_LU_PP(n, p, A)
    ! "In place" LU factorization with PARTIAL pivoting of an n-by-n single
    ! precision matrix A. Inputs and outputs should be as follows
    !
    ! Inputs:
    ! -------
    ! integer :: n - the size of A
    ! integer :: p - a dummy vector with n entries to store the row permutation of A
    ! real    :: A - the matrix A, datatype real (single precision)
    !
    ! Outputs:
    ! --------
    ! integer :: p - a permutation vector such that A(p,:) = L*U, produced by
    !                the partial pivoting procedure.
    ! real    :: A - a real n-by-n array, overwritten onto the input matrix A,
    !                containing the LU factors of the permuted matrix A(p,p).
    !                Precisely:
    !                (1) The strictly lower triangular part of L is stored in the
    !                    strictly lower triangular part of A
    !                (2) The upper triangular factor U is stored in the upper
    !                    triangular part of A.

    ! ********************************************************************** !
    ! Write the code for this subroutine. Include clear and informative comments.
    ! Add arguments to the subroutine declaration as needed.
    ! Replace these comment lines with a description of the subroutine, and
    ! include insructions for the user.
    ! ********************************************************************** !
	implicit none
    integer, intent(in) :: n
    real(4) :: A(n,n), p(1,n), A_temp(n,n), p_temp(1,n) , A_max
    integer :: k, j , r , i  
	
    do k = 1 , n-1
		A_max = abs(A(k,k));
		 r = k;
			do i = k , n
				if (abs(A(i,k)) > A_max) then
						r = i;
						A_max = abs(A(i,k));
				end if 
			enddo
		if (A(r,k)  == 0) then
				print *, "Error: A is not invertiable"
		else
				A_temp = A;
				A(k,:) = A(r,:);
				A(r,:) = A_temp(k,:);
				
				p_temp = p;
				p(1,k) = p(1,r);
				p(1,r) = p_temp(1,k);
				
				do i = k+1 , n
						A(i,k) = A(i,k) / A(k,k);
						do j = k+1 , n
								A(i,j) = A(i,j) - A(i,k)*A(k,j);
						enddo
				enddo	
		end if
	enddo
	
	return

! END SUBROUTINE
end subroutine SP_LU_PP

! ---------------------------------------------------------------------------- !

subroutine SP_LU_CP(n, p, q, A)
    ! "In place" LU factorization with COMPLETE pivoting of an n-by-n single
    ! precision matrix A. Inputs and outputs should be as follows
    !
    ! Inputs:
    ! -------
    ! integer :: n - the size of A
    ! integer :: p - a dummy vector with n entries to store the row permutation of A
    ! integer :: q - a dummy vector with n entries to store the column permutation of A
    ! real    :: A - the matrix A, datatype real (single precision)
    !
    ! Outputs:
    ! --------
    ! integer :: p, q - permutation vectors such that A(p,q) = L*U. The vectors
    !                   result from the complete pivoting procedure
    ! real    :: A    - a real n-by-n array, overwritten onto the input matrix A,
    !                   containing the LU factors of the permuted matrix A(p,p).
    !                   Precisely:
    !                   (1) The strictly lower triangular part of L is stored in
    !                       the strictly lower triangular part of A
    !                   (2) The upper triangular factor U is stored in the upper
    !                       triangular part of A.

    ! ********************************************************************** !
    ! Write the code for this subroutine. Include clear and informative comments.
    ! Add arguments to the subroutine declaration as needed.
    ! ********************************************************************** !
	implicit none
    integer, intent(in) :: n
    real(4) :: A(n,n), p(1,n), q(1,n) , A_temp(n,n), p_temp(1,n) , q_temp(1,n) , A_max
    integer :: k, j , r , i , m 
	
    do k = 1 , n-1
		A_max = abs(A(k,k));
					r = k;
					m = k;
			do i = k , n
				do j = k , n
						if (abs(A(i,j)) > A_max) then
							r = i;
							m = j;
							A_max = abs(A(i,j));
						end if 
				 enddo
			enddo
			
			
			
		if (A(r,m)  == 0) then
				print *, "Error: A is not invertiable"
		else
				A_temp = A;
				A(k,:) = A(r,:);
				A(r,:) = A_temp(k,:);
				
				A_temp = A;
				A(:,k) = A(:,m);
				A(:,m) = A_temp(:,k);
				
				p_temp = p;
				p(1,k) = p(1,r);
				p(1,r) = p_temp(1,k);
				
				q_temp = q;
				q(1,k) = q(1,m);
				q(1,m) = q_temp(1,k);
			    !print *, q;
				
				do i = k+1 , n
						A(i,k) = A(i,k) / A(k,k);
						do j = k+1 , n
								A(i,j) = A(i,j) - A(i,k)*A(k,j);
						enddo
				enddo	
		end if
	enddo
	
	return

! END SUBROUTINE
end subroutine SP_LU_CP

! ---------------------------------------------------------------------------- !

subroutine SP_UNIT_LSOLVE(n, L, b)
    ! Subroutine to solve a unit lower triangular system Lx=b with L unit lower
    ! triangular. Only the strictly lower triangular part of L is used since
    ! the diagonal entries are known to be 1. The rest of L is ignored (so it
    ! can be anything!). Solution is overwritten onto b
    !
    ! Inputs:
    ! -------
    ! integer :: n - size of L and b
    ! real    :: L - the coefficient matrix (an n-by-n real array).
    ! real    :: b - the right-hand-side vector (an 1D real array with n elements)
    !
    ! Outputs:
    ! -------
    ! real :: b - the solution vector (overwrites b)

    ! Declare all variables, no implicits
    implicit none
    integer, intent(in) :: n
    real(4), intent(in) :: L(n,n)
    real(4) :: b(n)
    integer :: i, j

    ! Loop, starting from second entry: x(1)=b(1) since L is unit lower triang.
    ! General formula: x(i) = b(i) - L(i,1)*x(1) - ... - L(i,i-1)x(i-1)
    do i = 2, n
        do j = 1, i-1
            b(i) = b(i) - L(i,j)*b(j)
        end do
    end do

! END SUBROUTINE
end subroutine SP_UNIT_LSOLVE

! ---------------------------------------------------------------------------- !

subroutine SP_USOLVE(n, U, b)
    ! Subroutine to solve an upper triangular system Ux=b with U upper
    ! triangular. The input U comes from low-memory versions of LU factorization
    ! algorithms, implemented in the following functions:
    ! SP_LU: LU factorization with no pivoting (single precision)
    ! SP_LU_PP: LU factorization with partial pivoting (single precision)
    ! SP_LU_CP: LU factorization with complete pivoting (single precision)
    !
    ! Inputs:
    ! -------
    ! integer :: n - size of L and b
    ! real    :: U - the coefficient matrix (an n-by-n real array).
    ! real    :: b - the right-hand-side vector (an 1D real array with n elements)
    !
    ! Outputs:
    ! -------
    ! real :: b - the solution vector (overwrites b)

    ! Declare all variables, no implicits
    implicit none
    integer, intent(in) :: n
    real(4), intent(in) :: U(n,n)
    real(4) :: b(n)
    integer :: i, j

    ! Loop moving backwards
    ! General formula: x(i) = ( b(i) - U(i,i+1)*x(i+1) - ... - U(i,n)x(n) ) / U(i,i)
    b(n) = b(n)/U(n,n)
    do i = n-1, 1, -1
        do j = i+1, n
            b(i) = b(i) - U(i,j)*b(j)
        end do
        b(i) = b(i)/U(i,i)
    end do

! END SUBROUTINE
end subroutine SP_USOLVE
