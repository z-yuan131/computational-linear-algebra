function x = Rtsolve(R,b)

% x = Rtsolve(R,b) solves the lower triangular system R'*x=b using forward
% substitution. This code avoids transposing R.
%
% G. Fantuzzi, 28 Nov 2019

% Preliminary stuff: size of R and preallocate solution (column vector)
% for speed
n = size(R,1);
x = zeros(n,1);

% Compute: forward substitution for R'*x = b with no error checks
% Must take into account that R in this function has not been transposed,
% we work with its columns instead of its rows. Use vector operations as
% much as possible
x(1) =  b(1)/R(1,1);    % Use R'(1,1) = R(1,1)
for i = 2:n
    % Use R'(i,1:i-1) = R(1:i-1,i) and R'(i,i) = R(i,i)
%     x(i) =  ( b(i) - dot( R(1:i-1,i), x(1:i-1) ) ) / R(i,i);
    x(i) =  ( b(i) - dot( R(i-1,i), x(i-1) ) ) / R(i,i);   %%tridiagonal matrix
    
end