% This function addresses the "Discrete Modal Analysis" task. 
% It identifies a specific mode of a mass-string system using the 
% generalized eigenvalue problem and simulates resonance by forcing 
% the system at that resonant frequency.
function discrete_modal_analysis()

    string_params = struct();

    n = 10; % Number of masses 
    mode_index = 6; % Specific harmonic mode to excite

    string_params.n = n; 
    string_params.M = 10; % Total mass M_total 
    string_params.Tf = 2; % Tension force Tf 
    string_params.L = 7; % Total length L 
    string_params.c = 0.0001; % Damping "fudge factor"
    string_params.dx = string_params.L/(n+1); % Spacing delta_x 

    [M_mat,K_mat] = construct_2nd_order_matrices(string_params);

    [Ur_mat,lambda_mat] = eig(K_mat,M_mat);

    % Extract the mode shape and calculate resonant frequency omega_n 
    mode_shape_LA = Ur_mat(:,mode_index);
    omega_n = sqrt(-lambda_mat(mode_index, mode_index)); 

    omega = omega_n; 
    A = 3; % Forcing amplitude

    % Define forcing function uf(t) and its derivative
    Uf_func = @(t_in) A*sin(omega*t_in);
    dUfdt_func = @(t_in) -omega*A*cos(omega*t_in);
    string_params.Uf_func = Uf_func; 
    string_params.dUfdt_func = dUfdt_func;
    string_params.Uf_amplitude = A;

    % Initial conditions: zero displacement and velocity
    U0 = zeros(n,1);
    dUdt0 = zeros(n,1);
    V0 = [U0;dUdt0];

    tlist = linspace(0,20*(2*pi)/omega,5000+1);

    % Integration using the rate function derived from Newton's Second Law 
    my_rate_func = @(t_in, V_in) string_rate_func02(t_in, V_in, string_params);
    [tlist,Vresult] = ode45(my_rate_func, tlist, V0);

    % Visualization comparing predicted mode shape to simulation 
    animate_spring_with_mode_shape(tlist, Vresult, string_params, mode_shape_LA);

end