% This function addresses the "Traveling Wave" task.
% It uses a large number of masses to approximate a continuous 
% string and demonstrates how a pulse translates and reflects at 
% the wavespeed c = sqrt(Tf/rho).
function traveling_waves()

    string_params = struct();
    n = 200; % High n to approximate continuous string 

    string_params.n = n;
    string_params.M = 10;
    string_params.Tf = 2;
    string_params.L = 7;
    string_params.c = 0.0001;
    string_params.dx = string_params.L/(n+1);

    % Boundary conditions for traveling wave experiment 
    string_params.Uf_func = @(t_in) 0; 
    string_params.dUfdt_func = @(t_in) 0;
    string_params.Uf_amplitude = 0;

    % Calculate wavespeed
    rho = string_params.M/string_params.L;
    c = sqrt(string_params.Tf/rho);

    x_list = linspace(0,string_params.L,n+2)';
    x_list = x_list(2:end-1);

    w_pulse = string_params.L/5;
    h_pulse = 7;
    x_shift = 0;

    % Set initial conditions to a pulse shape
    U0 = b_spline_pulse(x_list-x_shift,w_pulse,h_pulse);
    % Initial velocity set to move pulse in one direction
    dUdt0 = -c * b_spline_pulse_derivative(x_list-x_shift,w_pulse,h_pulse);

    figure(1);
    hold on
    plot(x_list,U0, 'o-', 'color', 'r');
    plot(x_list,dUdt0, 'o-','color','b');

    V0 = [U0;dUdt0];
    tlist = linspace(0,6*string_params.L/c,10000+1);

    my_rate_func = @(t_in, V_in) string_rate_func02(t_in, V_in, string_params);
    [tlist,Vresult] = ode45(my_rate_func, tlist, V0);

    % Animate with a tracking line moving at speed c [Page 10, 11]
    animate_spring_travelling_wave(tlist, Vresult, string_params, x_shift+w_pulse/2, c);

end