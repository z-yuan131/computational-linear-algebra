function R = cholFact(A)

% R = cholFact(A) returns the upper triangular Cholesky factor of a
%     positive semidefinite symmetric matrix A. This function only operates
%     on the upper-triangular part of A and implements the CHolesky
%     algorithm.
%
% G. Fantuzzi, 28 Nov 2019

% Check inputs
[m, n] = size(A);
assert(m==n, 'Incorrect input size: A must be n-by-n')

% Zero all entries below the main diagonal
R = triu(A);

% Operate on the upper triangular part using Cholesky algorithm (see
% Algorithm 4.11 of the lecture notes)
for i = 1:n
   if R(i,i) > 0
      for j = (i+1):min(i+2,n)             %%%tridiagonal matrix R
         R(j,j:n) = R(j,j:n) - R(i,j:n).*( R(i,j)/R(i,i) );
      end

      R(i,i:n) = R(i,i:n)./sqrt(R(i,i));
   else
       error('A is not positive definite!')
   end
end