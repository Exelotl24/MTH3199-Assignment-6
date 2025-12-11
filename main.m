

function main()

    % Create string_params

    % Width for tracking line
    w = .1;
    h = 0.005;

    num_masses = 100;
    total_mass = 1;                    % kg
    tension_force = 10;                % N
    string_length = 5;                  % m
    damping_coeff = 0.001;
    amplitude_Uf = 1;
    omega_Uf = 10;
    string_params = create_string_params(num_masses, total_mass, tension_force, string_length, damping_coeff, amplitude_Uf, omega_Uf, w, h);
        
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



    % Define Dormand Prince BT Struct
    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
        1/5, 0, 0, 0,0,0,0;...
        3/40, 9/40, 0, 0, 0, 0,0;...
        44/45, -56/15, 32/9, 0, 0, 0,0;...
        19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
        9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
        35/384, 0, 500/1113, 125/192, -2187/6784,11/84,0];


    U0 = zeros(num_masses, 1);
    dUdt0 = zeros(num_masses, 1);
    V0 = [U0; dUdt0];
    tspan = [0 10];
    string_simulation_template01(string_params, @string_rate_func01, U0, dUdt0, tspan)

%-------------------------- MODAL ANALYSIS --------------------------------
    
    [M_mat, K_mat] = construct_2nd_order_matrices(string_params);
   
    [U_modes, lambda_mat] = eig(K_mat, M_mat);
    lambda = diag(lambda_mat);
    omega = sqrt(lambda);

    % Pick mode shape
    r = 2;
    Ur = U_modes(:, r);     % mode shape
    wr = omega(r);          % resonant frequency

    % Normalize
    Ur = Ur / max(abs(Ur));

    omega_Uf = wr * 0.99;   % slight detuning to show resonance


    % Initial conditions
    U0 = zeros(num_masses, 1);
    dUdt0 = zeros(num_masses, 1);
    V0 = [U0; dUdt0];


    my_rate = @(t, V) string_rate_func01(t, V, string_params);

    XA = V0;                
    t0 = 0;
    h = 0.005;              
    tspan = [0,20];
    Nt = length(tspan);
    
    tlist = [0: 0.01: 20];                  
    Vlist = zeros(length(V0), Nt);  
    Vlist(:,1) = XA;
    xlist = 0:string_params.dx:string_params.L;

    % Run integration
    [t_list,Vlist,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
    (my_rate,tspan,V0,h,DormandPrince);

    % % Time-stepping
    % for k = 1:Nt-1
    %     t = tspan(k);                        
    %     [XB, ~] = explicit_RK_step(my_rate, t, XA, h, DormandPrince);        
    %     XA = XB;                              
    %     V_list(:, k+1) = XA;                 
    % end
    % [XB, num_evals] = explicit_RK_step(my_rate, t, XA, h, DormandPrince); 

    animate_title = "Vibrating String Animation";
    animate_string(tlist, Vlist, string_params, animate_title, c, w);
    subplot(1,2,1)
    plot(xlist(2:end-1), Ur)

    
% ------------------ CONTINUOUS VS. DISCRETE ANALYSIS ------------------

    % continuous_vs_discrete_analysis(string_params, 3);     

% ------------------------- TRAVELLING WAVES ---------------------------
    
    % traveling_wave_setup(string_params, DormandPrince, V0, c, w);
    

end