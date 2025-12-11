function traveling_wave_setup(initial_params, DormandPrince, V0, c, w)


    % New string params
    string_params = initial_params;
    n = string_params.n;
    L = string_params.L;
    dx = string_params.dx;
  
    
    % Fix right point
    string_params.amplitude_Uf = 0;
    string_params.omega_Uf = 0; 
    
    % Redefine the function handles 
    string_params.Uf_func = @(t_in) 0;
    string_params.dUfdt_func = @(t_in) 0;
    
    %  Run simulation
    L = string_params.L;

    t0 = 0;
    t_end = 2 * L / c * 1.5; 
    h = 0.005;               
    tspan = t0:h:t_end;
    Nt = length(tspan);
    
    my_rate = @(t, V) string_rate_func01(t, V, string_params);
    
    % XA = V0;
    % tlist = tspan;
    % Vlist = zeros(length(V0), Nt);
    % Vlist(:, 1) = XA;
    % Energy_list = zeros(1, Nt);
    % 

    % Initial conditions
    U0 = zeros(n, 1);
    dUdt0 = zeros(n, 1);
    V0 = [U0; dUdt0];

    % Calculate initial energy
    [M_mat, K_mat] = construct_2nd_order_matrices(string_params);
    Energy_list(1) = 0.5 * (dUdt0' * M_mat * dUdt0 + U0' * K_mat * U0);
    % 
    % % Time step
    % for k = 1:Nt-1
    %     t = tspan(k);
    %     [XB, ~] = explicit_RK_step_func(my_rate, t, XA, h, DormandPrince);
    %     XA = XB;
    %     Vlist(:, k+1) = XA;
    % 
    %     U = XA(1:n);
    %     dUdt = XA(n+1:2*n);
    %     Energy_list(k+1) = 0.5 * (dUdt' * M_mat * dUdt + U' * K_mat * U);
    % end

    XA = V0;                
    t0 = 0;
    h = 0.005;              
    tspan = [0,20];
    Nt = length(tspan);
    
    tlist = tspan;                  
    Vlist = zeros(length(V0), Nt);  
    Vlist(:,1) = XA;
    xlist = 0:string_params.dx:string_params.L;

    % Run integration
    [t_list,Vlist,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
    (my_rate,tspan,V0,h,DormandPrince);


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