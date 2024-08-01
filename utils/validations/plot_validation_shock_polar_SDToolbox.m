function [ax1, ax2] = plot_validation_shock_polar_SDToolbox(mixArray1, mixArray2, results_SDToolbox, config)
    % Plot numerical results obtained with SDToolbox, which use CANTERA as a thermochemical kernel.
    %
    %   * Pressure ratio with the deflection angle [deg]
    %   * Wave angle [deg] with the deflection angle [deg]
    
    % Import packages
    import combustiontoolbox.utils.display.plotPolar
    
    % Shock polars from Combustion Toolbox
    [ax1, ax2] = plotPolar(mixArray1, mixArray2);

    % Shock polars from SDToolbox
    nfrec = [2, 4, 7, 10]; % frec points SD Toolbox per case

    for i = length(results_SDToolbox):-1:1
        plot(ax1, results_SDToolbox{i}.theta(1:nfrec(i):end) * 180 / pi, results_SDToolbox{i}.P2P1(1:nfrec(i):end), 'ko', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
        plot(ax2, results_SDToolbox{i}.theta(1:nfrec(i):end) * 180 / pi, results_SDToolbox{i}.beta(1:nfrec(i):end) * 180 / pi, 'ko', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    end

    % Set x-axis limits for wave angle-deflection polar
    xlim(ax2, [0, ax2.XLim(2)]);
end
