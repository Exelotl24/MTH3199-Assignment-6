function string_simulation_template01()
    num_masses = 4;
    total_mass = 10;                    % kg
    tension_force = 100;                % N
    string_length = 5;                  % m
    damping_coeff = 0.05;
    dx = string_length/(num_masses+1);
    amplitude_Uf = %your code here
    omega_Uf = %your code here
    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);
    Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);
    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);
    %initial conditions
    U0 = %your code here
    dUdt0 = %your code here
    V0 = [U0;dUdt0];
    tspan = %your code here
    %run the integration
    % [tlist,Vlist] = your_integrator(my_rate_func,tspan,V0,...);
    %your code to generate an animation of the system
end