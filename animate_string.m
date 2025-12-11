function animate_string(tlist, Vlist, string_params, animate_title, c, w)
    
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
        
        % Update plot
        set(h, 'YData', U_full);
        drawnow;
        pause(0.005);
    end

end
