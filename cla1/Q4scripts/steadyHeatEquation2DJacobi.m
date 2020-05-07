%%%% Function for Jacobi method
function [T,j] = steadyHeatEquation2DJacobi(A,b,TOL)

%%produce diagonal entries
D = diag(diag(A));
N = A - D;

%%initialization
T = sparse(size(A,1),1);
k=10;
j=0;


while (k > TOL)
    
    
    T = D\(b -N*T);
    
    k = norm(A*T - b)/norm(b);

    
    
j = j + 1;
end
