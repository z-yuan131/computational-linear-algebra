%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all


%% Question 1
%%initialization
M = [3.2525,0.6486,0,0;0.6486,2.0868,0.4151,0;0,0.4151,1.3320,0.2642;0,0,0.2642,0.4709];
K = [9.4955,-4.0024,0,0;-4.0024,6.9048,-2.9024,0;0,-2.9024,4.9956,-2.0932;0,0,-2.0932,2.0932];

M = sparse(M);
K = sparse(K);
%%solver
% profile clear
% profile -memory on
[eigvec,eigval]=eigs(K,M);   %%inbuilt function eig is used to calculate eigenvalue problem
wn = sqrt(diag(eigval));    %%natural frequency  is squre root of eigvalue

n = 4;

% % n = input('system size n = ');




%% RQI && HHD
[x_RQI,sigma_RQI] = RQI(n,M,K);


%% QRITER
[x_QR,sigma_QR] = QRITER(n,M,K);


%% SSI && RITZ
[x_SSI,sigma_SSI] = SSI(n,M,K);
% profile report
%% SENS
[delta_rho,delta_e] = SENS(eigval,eigvec,n);


