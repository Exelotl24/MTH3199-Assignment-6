% Visualizes the string and overlays the theoretical mode shape 
function animate_spring_with_mode_shape(tlist, Vresult, string_params, mode_shape_LA)
    n = string_params.n;
    L = string_params.L;
    x_list = linspace(0,L,n+2);
    maxU = max(max(Vresult(:,1:n)));
    minU = min(min(Vresult(:,1:n)));
    h = max([abs(maxU), abs(minU),string_params.Uf_amplitude]);

    axis([0,L, -1.1*h,1.1*h]);
    hold on
    string_plot = plot (0,0, 'o-', 'color', 'k', 'LineWidth', 2, 'MarkerFaceColor','r', 'MarkerEdgeColor','r', 'markersize', 4);
    mode_shape_plot = plot (0,0, 'o-', 'color', 'b', 'LineWidth', 2, 'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'markersize', 4);

    % Scale mode shape for visual comparison
    mode_shape_padded = [0;mode_shape_LA;0];
    scale_factor = max(abs(maxU),abs(minU))/max(abs(mode_shape_padded));
    mode_shape_padded = mode_shape_padded * scale_factor;
    set(mode_shape_plot, 'xdata', x_list, 'ydata', mode_shape_padded);

    xlabel('x'); ylabel('U(t,x)'); title('Resonant Mode Simulation');

    for k = 1:length(tlist)
        t = tlist(k);
        U = Vresult(k,1:n);
        Uf = string_params.Uf_func(t);
        U_padded = [0, U, Uf];
        set(string_plot, 'xdata', x_list, 'ydata', U_padded);
        drawnow;
    end
end