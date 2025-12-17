% Visualizes a pulse and a tracking line 
function animate_spring_travelling_wave(tlist, Vresult, string_params, x_centroid0, c)
    n = string_params.n;
    L = string_params.L;
    x_list = linspace(0,L,n+2);
    maxU = max(max(Vresult(:,1:n)));
    minU = min(min(Vresult(:,1:n)));
    h = max([abs(maxU), abs(minU),string_params.Uf_amplitude]);

    axis([0,L, -1.1*h,1.1*h]);
    hold on
    string_plot = plot(0,0, 'o-', 'color', 'k', 'LineWidth', 2, 'MarkerFaceColor','r', 'MarkerEdgeColor','r', 'markersize', 4);
    centroid_plot = plot(0,0, 'b','linewidth',2);

    xlabel('x'); ylabel('U(t,x)'); title('Traveling Wave Simulation');

    for k = 1:length(tlist)
        t = tlist(k);
        U = Vresult(k,1:n);
        Uf = string_params.Uf_func(t);
        U_padded = [0, U, Uf]; % boundary conditions

        % Tracking line
        x_centroid = c*t+x_centroid0;
        x_centroid = mod(x_centroid,2*L);
        if x_centroid > L
            x_centroid = 2*L - x_centroid;
        end

        set(centroid_plot, 'xdata', [x_centroid,x_centroid],'ydata',[-1.1*h,1.1*h]);
        set(string_plot, 'xdata', x_list, 'ydata', U_padded);
        drawnow;
    end
end