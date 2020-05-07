%%%%%%CLA course work Qestion 1%%%%%%%%
%%%%%%By Zhenyang Yuan 01793990
%%%%%%21/01/2020

clc;
clear;
close all;

%%
n = 100;
deltax = 1/(n+1);
deltat = 0.01;
x = 0+deltax:deltax:1-deltax;
B = diffusionMatrix1D(n);

%%
sigema = deltat/deltax^2;
E = sparse(eye(n));
A = E + sigema/2*B;
b0 = E - sigema/2*(B' + triu(B,1));

%%
[m, n] = size(A);
assert(m==n, 'Incorrect input size: A must be n-by-n')

R = cholFact(A);
u = 1-cos(4*pi*x);
% plot(x,u);    %t = 0.0
u = u';
for i = 0:deltat:1.5  %t = 0.1,0.4,1.5
    b = b0*u;
    y = Rtsolve(R,b);         %% R'y=b
    u = Rsolve(R,y);          %% Ru=y

    figure(1)
    plot(x,u);
    ylim([0,2]);
    drawnow
end
