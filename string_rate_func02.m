% Implementation of the rate function
function dVdt = string_rate_func02(t,V,string_params)
    n = string_params.n; 
    M = string_params.M; 
    Uf_func = string_params.Uf_func; 
    dUfdt_func = string_params.dUfdt_func;
    Tf = string_params.Tf; 
    c = string_params.c; 
    dx = string_params.dx; 

    % Unpack state 
    U = V(1:n);
    dUdt = V((n+1):(2*n));
    Uf = Uf_func(t);
    dUfdt = dUfdt_func(t);

    % Handle boundary conditions u1 and un 
    U_left = [0;U(1:end-1)];
    U_right = [U(2:end);Uf];

    dUdt_left = [0;dUdt(1:end-1)];
    dUdt_right = [dUdt(2:end);dUfdt];

    % Acceleration term
    term1 = (Tf/dx)*(U_left-2*U+U_right);
    term2 = (c/dx)*(dUdt_left-2*dUdt+dUdt_right);

    % d2U/dt2
    d2Udt2 = (term1+term2)/(M/n); 
   
    dVdt = [dUdt;d2Udt2]; % Assemble state derivative 
end