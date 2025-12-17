% This function addresses "From Discrete to Continuous" and 
% "Modal Analysis of the Wave Equation". It compares the 
% mode shapes and resonant frequencies of the discrete system 
% to the continuous system as n increases.
function limiting_behavior()

    string_params = struct();
    mode_index = 3;

    string_params.M = 10;
    string_params.Tf = 2;
    string_params.L = 7;
    string_params.c = 0.0001;

    rho = string_params.M/string_params.L;
    c = sqrt(string_params.Tf/rho);

    % Continuous theoretical frequency and shape
    omega_n_spatial = pi * mode_index / string_params.L;
    omega_n = c * omega_n_spatial;

    x_list_continuous = linspace(0, string_params.L,1000);
    mode_shape_WA = sin(omega_n_spatial * x_list_continuous);

    % List of masses to test for convergence
    n_list = [5,10,30,50];
    omega_list = [];

    figure(1);
    for k = 1:length(n_list)
        n = n_list(k);
        string_params.n = n;
        string_params.dx = string_params.L/(n+1);
    
        [M_mat,K_mat] = construct_2nd_order_matrices(string_params);
        [Ur_mat,lambda_mat] = eig(K_mat,M_mat);
    
        % Extract discrete mode shape and rescale for comparison
        mode_shape_LA = [0;Ur_mat(:,n+1-mode_index);0];
        mode_shape_LA = mode_shape_LA/max(abs(mode_shape_LA));
    
        % Capture discrete resonant frequency
        omega_n_WE = sqrt(-lambda_mat(n+1-mode_index, n+1-mode_index));
        omega_list(end+1) = omega_n_WE;
   
        x_list = linspace(0,string_params.L,n+2)';

        % Subplots comparing discrete markers to continuous line
        subplot(length(n_list),1,k)
        hold on
        plot(x_list_continuous,mode_shape_WA,'b', 'LineWidth',4);
        plot(x_list, mode_shape_LA,'o-', 'color', 'k', 'LineWidth', 2, 'MarkerFaceColor','r', 'MarkerEdgeColor','r', 'markersize', 4);
    end

    % Plot frequency approximation error vs number of masses
    figure(2);
    semilogx(n_list, omega_list, 'o-', 'color', 'k', 'LineWidth', 2, 'MarkerFaceColor','r', 'MarkerEdgeColor','r', 'markersize', 4);
    hold on
    semilogx(n_list, ones(length(n_list),1)*omega_n, 'b--', 'linewidth', 2);
    title('Convergence of Discrete Frequency to Continuous Frequency');
    xlabel('Number of Masses (n)');
    ylabel('Frequency (\omega)');
end
