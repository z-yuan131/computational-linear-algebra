%%%%%% CLA course work question 4 %%%%%%
%%%%%% Zhenyang Yuan 01793990 %%%%%%
clc
clear
close all

a = input('Input a: ');
nx = input('Input n_x: ');
ny = input('Input n_y: ');
choice = input('Input method choice: 1 = Jacobi; 2 = SOR ');
TOL = input('input tolerance (i.e. 10^(-8) or 10e-8) = ');

%%only for debug
% a = 1;
% nx = 5;     %%%nx
% ny = 10;       %%%ny
% choice = 2;

%%initial matrix b
b = sparse(3*nx*ny + nx + ny ,1);      %%size of elements
%%call function to calculate A
A = diffusionMatrix2D(a,nx,ny);

dx = 2*a/(2*nx+2);
dy = 2*a/(2*ny+2);
x1 = 0 : dx : 2 * a;
y1 = 0 : dy : 2 * a;
x = x1(2:length(x1) - 1);       %%without boundaries
y = y1(2:length(y1) - 1);

q = sparse(nan(2*ny+1 ,2*nx+1 ));

%%% q matrix is matrix version of vector b with nan elements in the domian
%%% that we don't calculate
for j = 1 : 2*nx+1
    if j <= nx
        for i = 1 : 2*ny+1      %%domain omega 1
            q(i,j) =  100 * cos(pi*x(j)/a) * sin(3*pi*y(i)/2/a) * exp(-0.5*((x(j)-a)^2+(y(i)-a)^2));
        end
    else
        for i = ny+2 : 2*ny+1   %%domain omega 2
            q(i,j) =  100 * cos(pi*x(j)/a) * sin(3*pi*y(i)/2/a) * exp(-0.5*((x(j)-a)^2+(y(i)-a)^2));
        end
    end
end

%%get rid of nan elements to produce vector b
q_t = reshape(q,[(2*nx+1)*(2*ny+1),1]);

k = 1;
for i = 1 : (2*nx+1)*(2*ny+1)
    if isnan(q_t(i))
        i = i+1;
    else
        b(k,1) = q_t(i);
        k = k + 1; 
    end;
end;


%% display number of iteration and running time
if (choice == 1)
    tic
    [T,n] = steadyHeatEquation2DJacobi(A,b,TOL);            %%call jacobi method
    toc
    disp('number of iterations for Jacobi')
    disp(n)
else
    tic
    [T,m] = steadyHeatEquation2Dsor(A,b,TOL);               %%call SOR method
    toc
    disp('number of iterations for SOR')
    disp(m)
end


%% contour
%%T_plot is for plot use, with boundaries; T_temp is for converting vector
%%T to oridinary matrix with normal dimensions, other elements stay nan.
T_plot = sparse(2*ny+3,2*nx+3);
T_temp = sparse(nan(2*ny+1,2*nx+1));

T_temp(:,1:nx) = reshape(T(1:(2*ny+1)*nx),[2*ny+1,nx]);       %%the index is a little tricky but okay with diagram.
T_temp(ny+2:2*ny+1,nx+1:2*nx+1) = reshape(T((2*ny+1)*nx+1:length(T)),[ny,nx+1]);
T_temp(1:ny,nx+1) = 0;                  %%boundary is equal to 0 rahter than nan;
T_temp(ny+1,nx+1:2*nx+1) = 0;

T_plot(2:2*ny+2,2:2*nx+2) = T_temp;     %%convert T_temp to whole domain for plot

[L, R] = meshgrid(x1,y1);
figure(2)
contourf(L,R,T_plot);
% mesh(L,R,T_plot);
colorbar('eastoutside')



