function [delta_rho,delta_e] = SENS(eigval,eigvec,n)
    m = zeros(2,2,n);
    m(:,:,1) = [2.2484,1.0087;1.0087,1.8052];
    m(:,:,2) = [1.4473,0.6486;0.6486,1.1594];
    m(:,:,3) = [0.9273,0.4151;0.4151,0.7410];
    m(:,:,4) = [0.5910,0.2642;0.2642,0.4709];

    k = zeros(2,2,n);
    k(:,:,1) = [5.4931,-5.4931;-5.4931,5.4931];
    k(:,:,2) = [4.0024,-4.0024;-4.0024,4.0024];
    k(:,:,3) = [2.9024,-2.9024;-2.9024,2.9024];
    k(:,:,4) = [2.0932,-2.0932;-2.0932,2.0932];

    sigma = diag(eigval);                           %%get eigenvalues
    x = eigvec;

    delta_rho = zeros(n,n);                         %%output data 
    delta_e = zeros(n,n);
    %%lamda to rho(M)                               %%operate dlambda/drho
    for j = 1:n
        for i = 1:n
            temp=zeros(n+1,n+1);                    %%global matrix dimension should be (n+1,n+1)
            temp(i:i+1,i:i+1) = m(:,:,i);
            temp(1,:)=[];temp(:,1)=[];              %%for fixed end boundary condition
            delta_rho(j,i) = -x(:,j)'*(sigma(j)*temp/3)*x(:,j);
        end
    end


    %%lamda to E(K)                                 %%operate dlambda/dE
    for j = 1:n
        for i = 1:n
            temp=zeros(n+1,n+1);
            temp(i:i+1,i:i+1) = k(:,:,i);
            temp(1,:)=[];temp(:,1)=[];
            delta_e(j,i) = x(:,j)'*(temp)*x(:,j)/4;
        end
    end
end
