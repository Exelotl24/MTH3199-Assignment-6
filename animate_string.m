function animate_string(tlist, Vlist, string_params, animate_title)
    
    % Extract data
    Ulist = Vlist(1:string_params.n, :);     
    
    xlist = linspace(0, string_params.L, string_params.n+2);  
    
    % Y-lims for animation
    maxU = max(abs(Ulist(:)));
    if maxU == 0
        maxU = 1;   
    end
    
    % Create plot
    figure;
    h = plot(xlist, zeros(size(xlist)), 'LineWidth', 2);
    grid on;
    xlabel('x');
    ylabel('displacement');
    title(animate_title);
    ylim(1.5 * [-maxU, maxU]);
    
    % Animation loop
    Nt = length(tlist);
    skip = max(1, floor(Nt/800)); 
    
    for k = 1:skip:Nt
        % Reconstruct string shape
        U_full = [0 ; Ulist(:,k) ; string_params.Uf_func(tlist(k))];
        
        % Update plot
        set(h, 'YData', U_full);
        drawnow;
        pause(0.01);
    end

end
