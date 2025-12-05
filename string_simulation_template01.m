function string_simulation_template01(string_params, string_rate_func, U0, dUdt0, tspan)
    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func(t_in,V_in,string_params);
    %initial conditions
    V0 = [U0;dUdt0];
    %run the integration
    % [tlist,Vlist] = your_integrator(my_rate_func,tspan,V0,...);
    %your code to generate an animation of the system
end