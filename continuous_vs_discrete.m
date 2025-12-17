% This function demonstrates the analytical solution to the wave equation.
% It uses the theoretical resonant frequency 
% to excite a high-resolution discrete model.
function continuous_vs_discrete()

    string_params = struct();
    n = 200;
    mode_index = 10; 

    string_params.n = n;
    string_params.M = 10;
    string_params.Tf = 2;
    string_params.L = 7;
    string_params.c = 0.0001;
    string_params.dx = string_params.L/(n+1);

    rho = string_params.M/string_params.L;
    c = sqrt(string_params.Tf/rho);

    % Spatial frequency 
    omega_n_spatial = pi * mode_index / string_params.L;
    % Temporal frequency
    omega_n = c * omega_n_spatial;

    x_list = linspace(0,string_params.L,n+2)';
    x_list = x_list(2:end-1);

    % Theoretical mode shape
    mode_shape_WA = sin(omega_n_spatial * x_list);

    omega = omega_n; 
    A = 3;

    string_params.Uf_func = @(t_in) A*sin(omega*t_in);
    string_params.dUfdt_func = @(t_in) -omega*A*cos(omega*t_in);
    string_params.Uf_amplitude = A;

    U0 = zeros(n,1);
    dUdt0 = zeros(n,1);
    V0 = [U0;dUdt0];

    tlist = linspace(0,20*(2*pi)/omega,2000+1);
    my_rate_func = @(t_in, V_in) string_rate_func02(t_in, V_in, string_params);
    [tlist,Vresult] = ode45(my_rate_func, tlist, V0);

    animate_spring_with_mode_shape(tlist, Vresult, string_params, mode_shape_WA);

end
