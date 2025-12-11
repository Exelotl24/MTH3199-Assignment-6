function string_params = create_string_params(num_masses, total_mass, tension_force, string_length, damping_coeff, amplitude_Uf, omega_Uf, w, h)
    dx = string_length/(num_masses+1);
    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);

    %Vibrating Wave
   % Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
   % dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);

    %Triangle wave
    Uf_func = @(t_in) triangle_pulse(t_in,w,h);
    dUfdt_func = @(t_in) triangle_pulse_derivative(t_in,w,h);


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
    string_params.m = total_mass/num_masses;

end