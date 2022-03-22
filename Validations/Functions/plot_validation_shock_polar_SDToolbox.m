function [fig1, fig2] = plot_validation_shock_polar_SDToolbox(results_SDToolbox, config)
    % Plot numerical results obtained with SDToolbox, which use CANTERA as a thermochemical kernel.
    %   * Pressure ratio with the deflection angle [deg]
    %   * Wave angle [deg] with the deflection angle [deg]
    nfrec = [2, 4, 7, 10];
    
    fig1 = figure(1); ax1 = gca;
    fig2 = figure(2); ax2 = gca;
    set(fig1, 'position', [1921 -471 1080 1795]);
    set(fig2, 'position', [1921 -471 1080 1795]);
    for i=length(results_SDToolbox):-1:1
        plot(ax1, results_SDToolbox{i}.theta(1:nfrec(i):end) * 180/pi, results_SDToolbox{i}.P2P1(1:nfrec(i):end), 'ko', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
        plot(ax2, results_SDToolbox{i}.theta(1:nfrec(i):end) * 180/pi, results_SDToolbox{i}.beta(1:nfrec(i):end) * 180/pi, 'ko', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    end
end