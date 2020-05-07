function [x,sigma] = SSI(n,M,K)

    %% SSI and RITZ
    U = rand(n,n);
    %%cholechskey
    R = cholFact(M);                %% M = R'*R;
    a = R'\K;
    a = R\a;
    A = a;

    er = 100;
    k=0;
    while(er > 10^-4)
        U = a*U;
        %%RITZ
        M_bar = U'*M*U;
        K_bar = U'*K*U;
        [p,Lamda]=eig(K_bar,M_bar);

        q = U*p;
        U = q;
        k=k+1;
        er = norm(K_bar*p - M_bar*p*Lamda);
    end;

    sigma = diag(Lamda);
    %Inverse iteration to calculate eigenvectors
        x = ones(n,n);
        for i = 1:n
            er = 100;			%%initial error
            b = (A-sigma(i)*eye(n));
            while(er > 10^-4)		%%error control
                y = b\x(:,i);			%%inverse iteration
                x(:,i) = y/norm(y);	%%normalization
                er = norm((A-sigma(i)*eye(n))*x(:,i));
            end;
        end;
end  
        