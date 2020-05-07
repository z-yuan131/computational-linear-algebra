%% Function for building matrix A

function A = diffusionMatrix2D(a,nx,ny)

dx = 2*a/(2*nx+2); %%elements in x direction
dy = 2*a/(2*ny+2); %%elements in x direction




%%construct matrix A
%%% initialization matrix A
        B1 = sparse(2 * ny + 1, 2 * ny + 1);        %%diagonal elements in domain omega 1 which has elements 2 * ny + 1
        B2 = sparse(ny , ny);                       %%diagonal elements in domain omega 2 which has elements ny
        A = sparse(3*nx*ny + nx + ny , 3*nx*ny + nx + ny);   %%materix A which has elements 3*nx*ny + nx + ny
        A_ini1 = cell(nx , nx);                     %%cell materix in domain omega 1
        A_ini2 = cell(nx + 2 , nx + 2);             %%cell materix in domain omega 2
        
        %% construct tridiagonal matrix B in two domain omega 1(B1) and omega 2(B2) 
        for i = 1 : 2 * ny + 1                          %% 2 * ny + 1 elements each column in domain omega 1
             for j = 1 : 2 * ny + 1
                 if i == j
                   B1(i , j) = 2/dx^2+2/dy^2;          %%diagonal elements
                     if i < 2 * ny + 1                  %%locate non-diagonal locations
                            B1(i + 1 , j) = -1/dy^2;              %non-diagonal elements
                            B1(i , j + 1) = -1/dy^2;
                     end
                 end
             end
        end
        
        for i = 1 : ny                                  %% ny elements each column in domain omega 2
             for j = 1 : ny
                 if i == j
                   B2(i , j) = 2/dx^2+2/dy^2;
                     if i < ny     
                            B2(i + 1 , j) = -1/dy^2;       
                            B2(i , j + 1) = -1/dy^2;
                     end
                 end
             end
        end
%          figure(2)   %%debug
%          spy(B1);
%          figure(3)
%          spy(B2);
        
        %% construct tridiagonal cell matrix A_ini in two domain omega 1(A_ini1) and omega 2(A_ini2) 
        A_ini1(: , :) = {sparse(2 * ny + 1,2 * ny + 1)};%% 2 * ny + 1 elements each column in domain omega 1
        A_ini2(: , :) = {sparse(ny,ny)};                %% ny elements each column in domain omega 2
        for i = 1 : nx
                for j = 1 : nx
                        if i == j
                            A_ini1{i , j} = B1;                                   %diagonal cells 
                            if i < nx
                                A_ini1{i + 1 , j} = -sparse(1/dx^2*eye(2 * ny + 1));     %non-diagonal cells, use sparse matrix     
                                A_ini1{i , j + 1} = -sparse(1/dx^2*eye(2 * ny + 1));
                            end
                        end
                end
        end
                
        for i = 1 : nx + 2                                                         %extend one more cell for convience in building materix A
                for j = 1 : nx + 2
                        if i == j && i ~= 1
                            A_ini2{i , j} = B2;
                                A_ini2{i - 1 , j} = -sparse(1/dx^2*eye(ny));        
                                A_ini2{i , j - 1} = -sparse(1/dx^2*eye(ny));
                        end
                end

        end
        
        A1 = cell2mat(A_ini1);                      %% comvert cell matrix A_ini to oridinary matrix to calculate, same as follow
        A2 = cell2mat(A_ini2);
                
        %% build matrix A
        A(2*nx*ny+nx-ny+1 : 3*nx*ny + nx + ny , 2*nx*ny+nx-ny+1 : 3*nx*ny + nx + ny) = A2;         %%tricky index
        A(1:2*nx*ny+nx , 1:2*nx*ny+nx) = A1; 
         
        figure(1)     %%produce structure of A      
        spy(A);
     
        
        
        
        
        
 