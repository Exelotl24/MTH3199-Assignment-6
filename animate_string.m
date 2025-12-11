function animate_string(tlist, Vlist, string_params, animate_title)
    
    % Recalculate wave speed for initial velocity calculation
    rho = string_params.M / string_params.L;
    c = sqrt(string_params.Tf / rho);
    
    n = string_params.n;
    L = string_params.L;
    dx = string_params.dx;
    
    % Define x-coords
    x_masses = dx : dx : L - dx;
   
    % Initial conditions

    % Pulse params
    A = 1;      % Amplitude 
    xc = L/4;   % Center
    sigma = L/15; % Width 

    % Pulse time width for tracking line
    w = sigma / c;

    % Initial displacement 
    U0 = A * exp(-(x_masses - xc).^2 / (2 * sigma^2))';

    % Initial velocity
    dUdt0_prime = zeros(n, 1);

    % Central difference for interior points
    for i = 2:n-1
        dUdt0_prime(i) = (U0(i+1) - U0(i-1)) / (2 * dx);
    end

    % First mass diff
    dUdt0_prime(1) = (U0(2) - U0(1)) / dx;

    % Last mass diff
    dUdt0_prime(n) = (U0(n) - U0(n-1)) / dx;

    % Final initial velocity 
    dUdt0 = -c * dUdt0_prime;

    V0 = [U0; dUdt0];


    % Extract data
    Ulist = Vlist(1:string_params.n, :);     
    
    xlist = linspace(0, string_params.L, string_params.n+2);  
    
    % Y-lims for animation
    maxU = max(abs(Ulist(:)));
    if maxU == 0
        maxU = 1;   
    end
    
    % Create plot
    % figure();
    % h = plot(xlist, zeros(size(xlist)), 'LineWidth', 2);
    % grid on;
    % xlabel('x');
    % ylabel('displacement');
    % title(animate_title);
    % ylim(1.5 * [-maxU, maxU]);
   
    figure();
    
    h_string = plot(xlist, zeros(size(xlist)), 'LineWidth', 2, 'DisplayName', 'String');
    hold on;

    % Initialize tracking line 
    h_track = plot(0, 0, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Tracking Line'); 
    grid on;

    xlabel('$x$');
    ylabel('displacement');
    title(animate_title); 
    ylim(1.5 * [-maxU, maxU]);
    legend('show');

    % Animation loop
    Nt = length(tlist);
    skip = max(1, floor(Nt/800)); 
    
    for k = 1:skip:Nt

        t = tlist(k);

        % Reconstruct string shape
        U_full = [0 ; Ulist(:,k) ; string_params.Uf_func(tlist(k))];

        % Calculate tracking line (description of what I'm doing is shown
        % in commented code below

        % x = x-coord of tracking line, t = current time
        % c = wave speed, w = pulse width (in time), L = string length

        % Calculate the position pf the wave center assuming an infinitely unrolled string
        x_track_unmod = string_params.L - c * t + 0.5 * w * c;

        % Use mod 2L to account for reflections (mod returns the remainder after division of a by m, where a is the dividend and m is the divisor)
        x_track = mod(x_track_unmod, 2 * string_params.L);

        % "Fold" the position back into the 0-L domain due to fixed end
        if x_track > string_params.L
            x_track = 2 * string_params.L - x_track;
        end

        h = 0.005;

        % Update plot
        set(h_string, 'YData', U_full);

        drawnow;
        pause(0.005);
    end

end
