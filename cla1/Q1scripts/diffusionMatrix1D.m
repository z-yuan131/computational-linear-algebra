function B = diffusionMatrix1D(n)
x = 0 : 1/(n+1) : 1;
B = sparse(n);

for i = 1:n
    B(i,i) = 0.01 + exp(-10*(x(i+1)+x(i))/2) + 0.01 + exp(-10*(x(i+2)+x(i+1))/2);    %k = 0.01 + exp(-10*x);x loop from 0-1,N+2 items.
    if(i ~= n)
        B(i,i+1) = -(0.01 + exp(-10*(x(i+2)+x(i+1))/2));       %%%spalloc to speed up;
    end
end
% B = B' + triu(B,1);   %%B is symmetric, only store half to save memmory.
