function main()
    % Create string_params
    num_masses = 4;
    total_mass = 10;                    % kg
    tension_force = 100;                % N
    string_length = 5;                  % m
    damping_coeff = 0.05;
    amplitude_Uf = 1;
    omega_Uf = 1;
    string_params = create_string_params(num_masses, total_mass, tension_force, string_length, damping_coeff, amplitude_Uf, omega_Uf);


    U0 = 0;
    dUdt0 = 0.2;
    tspan = [0 10];
    string_simulation_template01(string_params, @string_rate_func01, U0, dUdt0, tspan)

end