function traveling_wave_setup(initial_params, explicit_RK_step_func, DormandPrince)

    % New string params
    string_params = initial_params;
    
    % Fix right point
    string_params.amplitude_Uf = 0;
    string_params.omega_Uf = 0; 
    
    % Redefine the function handles 
    string_params.Uf_func = @(t_in) 0;
    string_params.dUfdt_func = @(t_in) 0;

    % Recalculate wave speed for initial velocity calculation
    rho = string_params.M / string_params.L;
    c = sqrt(string_params.Tf / rho);
    
    n = string_params.n;
    L = string_params.L;
    dx = string_params.dx;
    
    % Define x-coords
    x_masses = dx : dx : L - dx;
    
    % % Initial conditions
    % 
    % % Pulse params
    % A = 1;      % Amplitude 
    % xc = L/4;   % Center
    % sigma = L/15; % Width 
    % 
    % % Pulse time width for tracking line
    % w = sigma / c;
    % 
    % % Initial displacement 
    % U0 = A * exp(-(x_masses - xc).^2 / (2 * sigma^2))';
    % 
    % % Initial velocity
    % dUdt0_prime = zeros(n, 1);
    % 
    % % Central difference for interior points
    % for i = 2:n-1
    %     dUdt0_prime(i) = (U0(i+1) - U0(i-1)) / (2 * dx);
    % end
    % 
    % % First mass diff
    % dUdt0_prime(1) = (U0(2) - U0(1)) / dx;
    % 
    % % Last mass diff
    % dUdt0_prime(n) = (U0(n) - U0(n-1)) / dx;
    % 
    % % Final initial velocity 
    % dUdt0 = -c * dUdt0_prime;
    % 
    % V0 = [U0; dUdt0];
    
    %  Run simulation
    t0 = 0;
    t_end = 2 * L / c * 1.5; 
    h = 0.005;               
    tspan = t0:h:t_end;
    Nt = length(tspan);
    
    my_rate = @(t, V) string_rate_func01(t, V, string_params);
    
    XA = V0;
    tlist = tspan;
    Vlist = zeros(length(V0), Nt);
    Vlist(:, 1) = XA;
    Energy_list = zeros(1, Nt);
    
    % Calculate initial energy
    [M_mat, K_mat] = construct_2nd_order_matrices(string_params);
    Energy_list(1) = 0.5 * (dUdt0' * M_mat * dUdt0 + U0' * K_mat * U0);

    % Time step
    for k = 1:Nt-1
        t = tspan(k);
        [XB, ~] = explicit_RK_step_func(my_rate, t, XA, h, DormandPrince);
        XA = XB;
        Vlist(:, k+1) = XA;
        
        U = XA(1:n);
        dUdt = XA(n+1:2*n);
        Energy_list(k+1) = 0.5 * (dUdt' * M_mat * dUdt + U' * K_mat * U);
    end
    
    % Animate
    animate_title = 'Traveling Wave Animation';
    animate_string(tlist, Vlist, string_params, animate_title, c, w);
    
    % Plot 
    figure();
    plot(tlist, Energy_list);
    xlabel('Time (s)');
    ylabel('Total Mechanical Energy (J)');
    title('Total Energy vs. Time for Traveling Wave');
    grid on;

end