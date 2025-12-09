%INPUTS
%t: current time
%V: system state. V = [U;dUdt] where
% U and dUdt are n x 1 column vectors
%string_params: a struct containing the system parameters describing the string
% string_params.n: number of masses
% string_params.M: total mass attached to the string
% string_params.Uf_func: function describing motion of end point
% string_params.dUfdt_func: time derivative of Uf
% string_params.Tf: %tension in string
% string_params.L: %length of string
% string_params.c: %damping coefficient
% string_params.dx: %horizontal spacing between masses
function dVdt = string_rate_func01(t,V,string_params)
    
    % Unpack essential parameters
    n = string_params.n;
    
    % Unpack function handles (CRITICAL FIX)
    Uf_func = string_params.Uf_func;
    dUfdt_func = string_params.dUfdt_func;

    % unpack state variable
    U = V(1:n);
    dUdt = V((n+1):(2*n));
    
    Uf = Uf_func(t);
    dUfdt = dUfdt_func(t);
    
    % compute acceleration
    d2Udt2 = zeros(n,1);  
    
    % left boundary (i=1)
    d2Udt2(1) = (string_params.Tf/string_params.dx * (-2*U(1) + U(2)) + string_params.c/string_params.dx  * (-2*dUdt(1) + dUdt(2))) / string_params.m;
    
    % interior points (i=2 to n-1)
    for i = 2:n-1
        d2Udt2(i) = (string_params.Tf/string_params.dx * (U(i-1) - 2*U(i) + U(i+1)) + string_params.c/string_params.dx * (dUdt(i-1) - 2*dUdt(i) + dUdt(i+1))) / string_params.m;
    end
    
    % right boundary (i=n) - Using 'n' which is now defined
    d2Udt2(n) = (string_params.Tf/string_params.dx * (U(n-1) - 2*U(n) + Uf) + string_params.c/string_params.dx * (dUdt(n-1) - 2*dUdt(n) + dUfdt)) / string_params.m;
    
    % assemble state derivative
    dVdt = [dUdt;d2Udt2];
end