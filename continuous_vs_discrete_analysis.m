function continuous_vs_discrete_analysis(initial_params, r_mode)

    M = initial_params.M;
    Tf = initial_params.Tf;
    L = initial_params.L;
    
    rho = M / L;
    cw = sqrt(Tf / rho);
    
    r = r_mode;          
    n_low = 10;        
    n_high = 50;        
    n_list_error = 5:5:100; 
    
    omega_r_cont = r * pi * cw / L; 
    
    figure();
    set(gcf, 'Position', [100, 100, 1200, 900]); 
    
    
    % ------------------ SUBPLOT: Mode Shape Comparison ------------------

    subplot(2, 2, 1);
    
    % Analysis for n_low 
    params_low = create_string_params(n_low, M, Tf, L, initial_params.c, 0, 0);
    [M_mat_low, K_mat_low] = construct_2nd_order_matrices(params_low);

    [U_modes_low, ~] = eig(K_mat_low, M_mat_low);
    Ur_discrete_low = U_modes_low(:, r);
    x_discrete_low = (1:n_low)' * params_low.dx;
    
    % Analysis for n_high 
    params_high = create_string_params(n_high, M, Tf, L, initial_params.c, 0, 0);
    [M_mat_high, K_mat_high] = construct_2nd_order_matrices(params_high);

    [U_modes_high, ~] = eig(K_mat_high, M_mat_high);
    Ur_discrete_high = U_modes_high(:, r);
    x_discrete_high = (1:n_high)' * params_high.dx;
    
    % --- Continuous Solution ---
    x_cont = linspace(0, L, 500);
    Ur_continuous = sin(r * pi * x_cont / L);
    
    % Rescale discrete mode shapes for comparison
    Ur_discrete_low = Ur_discrete_low / max(abs(Ur_discrete_low));
    Ur_discrete_high = Ur_discrete_high / max(abs(Ur_discrete_high));
    
    % Plotting
    plot([0; x_cont'; L], [0; Ur_continuous'; 0], '-', 'LineWidth', 2, 'DisplayName', 'Continuous');
    hold on;
    plot([0; x_discrete_low; L], [0; Ur_discrete_low; 0], 'r', 'DisplayName', ['Discrete (n=', num2str(n_low), ')']);
    plot([0; x_discrete_high; L], [0; Ur_discrete_high; 0], 'k.', 'DisplayName', ['Discrete (n=', num2str(n_high), ')']);
    
    xlabel('x');
    ylabel('Mode Shape u_r(x)');
    title(['Mode Shape Comparison (r=', num2str(r), ')']);
    legend('show');
    grid on;
    
    
    %------------ SUBPLOT 2: Frequency Approximation Error ---------------

    subplot(2, 2, 2);
    
    error_list = zeros(size(n_list_error));
    
    % Calculate convergence error
    for i = 1:length(n_list_error)
        n = n_list_error(i);
        
        params = create_string_params(n, M, Tf, L, initial_params.c, 0, 0);
        [M_mat, K_mat] = construct_2nd_order_matrices(params);
        
        % Find eigenvalues
        lambda = eig(K_mat, M_mat);
        omega = sqrt(lambda);
        omega_sorted = sort(omega);
        
        % Select frequency
        omega_r_discrete = omega_sorted(r);
        
        % Calculate Error
        error_list(i) = abs(omega_r_discrete - omega_r_cont);
    end
    
    semilogy(n_list_error, error_list, 'k.-', 'MarkerSize', 15);
    xlabel('Number of Masses (n)');
    ylabel('Absolute Frequency Error');
    title(['Approximation Error of omega_r vs. n (r=', num2str(r), ')']);
    grid on;
    
    
    %------------------ SUBPLOT 3: Harmonics Comparison -------------------

    subplot(2, 2, 3);
    
    % Calculate discrete frequencies
    n = n_low;
    params = create_string_params(n, M, Tf, L, initial_params.c, 0, 0);
    [M_mat, K_mat] = construct_2nd_order_matrices(params);
    lambda = eig(K_mat, M_mat);
    omega_discrete = sort(sqrt(lambda)); % The n discrete frequencies
    
    % Calculate continuous frequencies 
    k_harmonics = (1:n)';
    omega_cont = k_harmonics * pi * cw / L; 
    
    % Plot
    plot(k_harmonics, omega_discrete, 'r.', 'MarkerSize', 15, 'DisplayName', 'Discrete omega_k');
    hold on;
    plot(k_harmonics, omega_cont, 'b.', 'MarkerSize', 15, 'DisplayName', 'Continuous omega_k,cont');
    xlabel('Harmonic Index (k)');
    ylabel('Frequency (omega)');
    title(['Harmonics Comparison (n=', num2str(n), ')']);
    legend('show');
    grid on;
    
    
    %------------------ SUBPLOT 4: Harmonics Comparison ------------------

    subplot(2, 2, 4);
    
    % Calculate discrete frequencies
    n = n_high;
    params = create_string_params(n, M, Tf, L, initial_params.c, 0, 0);
    [M_mat, K_mat] = construct_2nd_order_matrices(params);
    lambda = eig(K_mat, M_mat);
    omega_discrete = sort(sqrt(lambda)); 
    
    % Calculate continuous frequencies 
    k_harmonics = (1:n)';
    omega_cont = k_harmonics * pi * cw / L; 
    
    % Plot
    plot(k_harmonics, omega_discrete, 'r.', 'MarkerSize', 15, 'DisplayName', 'Discrete omega_k');
    hold on;
    plot(k_harmonics, omega_cont, 'b.', 'MarkerSize', 15, 'DisplayName', 'Continuous omega_k,cont');
    xlabel('Harmonic Index (k)');
    ylabel('Frequency (omega)');
    title(['Harmonics Comparison (n=', num2str(n), ')']);
    legend('show');
    grid on;
    
end