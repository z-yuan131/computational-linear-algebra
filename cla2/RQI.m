%% subroutine RQI and HHD
function [x,sigma]=RQI(n,M,K)

    %%cholechskey
    R = cholFact(M);                %% M = R'*R;
    A = R'\K;
    A = R\A;
    
    
    AA = A;             %%for HHD will overwrite matrix A, AA is for calculating eigenvectors
    sigma = zeros(n,1);
%     xx = zeros(n,n);
    k = n;
    for i=1:4
        %%RQI
        x = ones(k,1);
        er = 100;
        while(er > 10^-4)
            sigma(i) = (x'*A*x)/(x'*x);   %%calculate eigenvalue
            b= A-sigma(i)*eye(k);
            y = b\x;                      %%inverse interation with shift for eigenvector
            x = y/norm(y,2);
            er = norm(A*x - sigma(i)*eye(k)*x);     %%error control to stop loop
        end

        %%HHD
    %     v = zeros(n,1);
    %     alpha = sqrt(sum(x.*x));
    %     if x(1) > 0
    %         alpha = -alpha;
    %     end
    %     v(1) = sqrt(0.5*(1-x(1)/alpha));
    %     for j = 2:n
    %         v(j) = sqrt(x(j)^2/(2*alpha^2*(1-x(1)/alpha)));
    %     end
    %     H = eye(n) - 2*(v)*(v');
    %     A = H*A/H;
        %%HHD
        if (k > 1)
            alpha = abs(norm(x));
            if x(1) > 0                 %%opposite sign to prevent leak from FLOPS
                alpha = -alpha;
            end
            e = zeros(k,1); e(1) = 1;   %%normal algrithem from lecture
            u = x - alpha*e;
            H = eye(k) - 2* (u)*(u')/(u'*u);
            A = H*A/H;
            
%             b = A(1,:);                 %%for calculating eigenvectors

            A(1,:)=[];                  %%for each iteration, size of A reduced by 1
            A(:,1)=[];
            
%             if k == n
%                xx(1:k,i) = x;              
%             else
%                xx(n-k+1:n,i) = x;         %%x here is just part of eigenvector;
%                yita = b*x/(sigma(i)-sigma(i-1));
%                xx(n-k,i) = yita;
%                xx(:,i) = HH\xx(:,i);
%             end
            k = k - 1;
%             HH = H;                       %%store H matrix for calculating eigenvectors
        end
    end


%Inverse iteration to calculate eigenvector
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