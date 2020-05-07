function x = Rsolve(R,b)

% x = Rsolve(R,b) solves the upper triangular system R*x=b using backward
% substitution. Only the upper triangular part of R is used.
%
% G. Fantuzzi, 28 Nov 2019

% Preliminary stuff: size of R and preallocate solution (column vector)
% for speed
n = size(R,1);
x = zeros(n,1);

% Compute: backward substitution with no error checks
% Implement the backward substitution formula using vector operations.
x(n) =  b(n)/R(n,n);
for i = n-1:-1:1
%     x(i) =  ( b(i) - dot( R(i,i+1:n), x(i+1:n) ) ) / R(i,i);
     x(i) =  ( b(i) - dot( R(i,i+1), x(i+1) ) ) / R(i,i);   %%%tridiagonal matrix
end