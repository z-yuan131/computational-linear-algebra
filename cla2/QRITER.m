%% QR
function [x,sigma] = QRITER(n,M,K)
    %%cholechskey
    R = cholFact(M);                %% M = R'*R;
    A = R'\K;
    A = R\A;
    AA = A;
    er = 100;

    U = zeros(n);
    e = zeros(n);
    while(er > 10^-4)
        for i = 1:n
             U(:,i) = A(:,i);
                for j = 1:i-1
                    U(:,i) = U(:,i) - dot(A(:,i),e(:,j))/norm(e(:,j))*e(:,j);
                end;
             e(:,i) = U(:,i)/norm(U(:,i));
        end;
        Q = e;
        A_t = A;
        A = Q'*A*Q;
        er = norm(A-A_t);
    end;

    sigma = diag(U);
    %%inverse itr
        x = ones(n,n);
        for i = 1:n
            er = 100;
            b = (AA-sigma(i)*eye(n));
            while(er > 10^-4)
                y = b\x(:,i);
                x(:,i) = y/norm(y);
                er = norm((AA-sigma(i)*eye(n))*x(:,i));
            end;
        end;
end