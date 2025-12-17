% Builds M and K matrices 
function [M_mat,K_mat] = construct_2nd_order_matrices(string_params)
    Tf = string_params.Tf; 
    dx = string_params.dx; 
    n = string_params.n; 
    M = string_params.M; 
    
    % Construct Discrete Laplacian 
    I_n = eye(n);
    M_left = circshift(I_n,[0,-1]);
    M_right = circshift(I_n,[0,1]);
    my_Laplacian = M_left - 2 * I_n + M_right;

    % Fix circular shift artifacts
    my_Laplacian(1,n) = 0;
    my_Laplacian(n,1) = 0;

    K_mat = (Tf/dx) * my_Laplacian;
    M_mat = (M/n) * I_n;
end