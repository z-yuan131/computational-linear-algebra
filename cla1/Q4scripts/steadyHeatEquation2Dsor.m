%%% Function for SOR method
function [T,j] = steadyHeatEquation2Dsor(A,b,TOL)

%%L lower diagonal matrix of A;U is upper diagonal matrix of A;
L = tril(A,-1);
D = tril(A) - L;
U = A - D - L;
N = U + L;
%%Find the optimal w for SOR iteration; function eigs gives largest six
%%eigenvaalues;rho_j means biggest eigenvalues of -M\N for Jacobi method.
%%from lecture note page 65.
rho_j = max(eigs(-D\N));      
w = 2/(1+sqrt(1-rho_j^2));


%%initialization
T = sparse(length(A),1);
k=10;
j=0;

while (k > TOL)
    
    T = (L + 1/w*D)\(b - (U + (1 - 1/w)*D)*T);
    
    k = norm(A*T - b)/norm(b);
    
    %%%Still, functions on the lecture note on page 63 can be used here to
    %%%calculate SOR. But it have to loop every index in the calculation
    %%%domain which is not very efficient way to do. So I choose to use
    %%%back slash to do SOR.
%     c = b - U * T - (1 - 1/w )* D * T;
%     a = 0;
%     T(1,1) = c(1,1)/D(1,1);
%         for i = 2 : length(A)
%             T(i,1) = w/D(i,i)*(c(i) -L(i,:)*T(:,1));
%         end;
% 
%         k = norm(A*T - b)/norm(b);
%         TT = T;
        

%         
j = j + 1;
end;
