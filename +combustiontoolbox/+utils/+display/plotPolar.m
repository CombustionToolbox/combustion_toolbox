function [ax1, ax2, ax3] = plotPolar(mixArray1, mixArray2, varargin)
    % Routine to obtain shock polar plots, namely:
    %
    %   * Plot (pressure, deflection)
    %   * Plot (wave angle, deflection)
    %   * Plot (velocity_x velocity_y)
    %
    % Args:
    %     mixArray1 (Mixture): Mixture object with the pre-shock state
    %     mixArray2 (Mixture): Mixture object with the post-shock state
    %
    % Optional Args:
    %     * mix2_case (Mixture): Mixture object with the post-shock state (case 2)
    %     * mix0 (Mixture): Mixture object with the post-shock state (case 0)
    %
    % Returns:
    %     Tuple containing
    %     
    %     * ax1 (object): Handle of the axes with (pressure, deflection) 
    %     * ax2 (object): Handle of the axes with (wave angle, deflection)
    %     * ax3 (object): Handle of the axes with (velocity_x velocity_y)
    %
    % Examples:
    %     * [ax1, ax2, ax3] = plotPolar(mixArray1, mixArray2)
    %     * [ax1, ax2, ax3] = plotPolar(mixArray1, mixArray2, mix2_case, mix0)
    
    % Import packages
    import combustiontoolbox.utils.display.PlotConfig
    
    % Unpack additional inputs
    if nargin > 2
        mix2_case = varargin{1};
        mix0 = varargin{2};
    else
        mix2_case = zeros(size(mixArray1));
        mix0 = mix2_case;
    end

    % Definitions
    config = PlotConfig();
    N = length(mixArray1);
    
    % Plots
    for i = N:-1:1
        % Plot (pressure, deflection)
        ax1 = plot_shock_polar_pressure(mixArray1(i), mixArray2(i), config, mix2_case(i), mix0(i));
        % Plot (wave angle, deflection)
        ax2 = plot_shock_polar_wave(mixArray1(i), mixArray2(i), config);
        % Plot (velocity_x velocity_y)
        ax3 = plot_shock_polar_velocities(mixArray2(i), config);
    end

end

% SUB-PASS FUNCTIONS
function ax = set_fixed_figure(fixed_number, config)
    % Generate a figure with a defined object identifier (fixed_number)
    fig = figure(fixed_number);
    ax = gca;
    set(fig, 'units', 'normalized', 'innerposition', config.innerposition, 'outerposition', config.outerposition);
    set(ax, 'LineWidth', config.linewidth, 'FontSize', config.fontsize - 2, 'BoxStyle', 'full');
    grid(ax, config.grid); box(ax, config.box); hold(ax, config.hold); ax.Layer = 'Top';
    set(ax, 'TickLabelInterpreter', 'latex');
    xlim(ax, config.axis_x);
    ylim(ax, config.axis_y);
end

function ax = plot_shock_polar_pressure(mix1, mix2, config, mix2_case, mix0)
    % Plot (pressure, deflection) - shock polar
    ax = set_fixed_figure(config.id_polar1, config);

    M1 = mix1.mach;
    p1 = pressure(mix1);
    Mach_label = '$\mathcal{M}_1 = ';
    xaxis = 0;

    if isobject(mix2_case)
        M1 = mix2_case.mach;
        p1 = mix0.p;
        Mach_label = '$\mathcal{M}_2 = ';
        xaxis = mix2_case.theta;
        y_lim = get(ax, 'YLim'); y_lim(2) = max(mix2.polar.p) / p1;
        plot(ax, [0, 0], y_lim, '-.k', 'LineWidth', config.linewidth)
    end

    plot(ax, xaxis + [- flip(mix2.polar.theta), mix2.polar.theta], [flip(mix2.polar.p), mix2.polar.p] / p1, 'k', 'LineWidth', config.linewidth, 'LineStyle', config.linestyle);
    % plot(ax, mix2.theta_max, mix2.polar.p(mix2.ind_max), 'kd', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    % plot(ax, -mix2.theta_max, mix2.polar.p(mix2.ind_max), 'kd', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    % plot(ax, mix2.polar.theta(mix2.ind_sonic), mix2.polar.p(mix2.ind_sonic), 'ks', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    text(ax, xaxis, max(mix2.polar.p), sprintf('%s%.2g$', Mach_label, M1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', config.fontsize - 4, 'Interpreter', 'latex');
    % text(ax, mix2.theta_max, str2.polar.p2(str2.ind_max), 'detachment point \rightarrow m', 'HorizontalAlignment', 'right', 'FontSize', config.fontsize-6)
    % text(ax, mix2.theta_sonic, str2.polar.p2(str2.ind_sonic), 'sonic point \rightarrow s', 'HorizontalAlignment', 'right', 'FontSize', config.fontsize-6)
    xlabel(ax, 'Deflection angle [deg]', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    ylabel(ax, '$p/p_1$', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    % ylim(ax, [1, max(str2.polar.p2)]);
end

function ax = plot_shock_polar_wave(mix1, mix2, config)
    % Plot (wave angle, deflection) - shock polar
    ax = set_fixed_figure(config.id_polar2, config);

    M1 = velocity_relative(mix1) / soundspeed(mix1);
    plot(ax, mix2.polar.theta, mix2.polar.beta, 'k', 'LineWidth', config.linewidth, 'LineStyle', config.linestyle);
    % plot(ax, mix2.theta_max, mix2.polar.beta(str2.ind_max), 'kd', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    % plot(ax, str2.polar.theta(str2.ind_sonic), mix2.polar.beta(str2.ind_sonic), 'ks', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    % text(ax, mix2.theta_max, mix2.polar.beta(str2.ind_max),
    % strcat('$\mathcal{M}_1 = ', sprintf('%.2g', M1), '$'), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', config.fontsize-6, 'Interpreter', 'latex');
    xlabel(ax, 'Deflection angle [deg]', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    ylabel(ax, 'Wave angle [deg]', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    %     xlim(ax, [0, mix2.theta_max]);
    ylim(ax, [0, 90.001]);
end

function ax = plot_shock_polar_velocities(mix2, config)
    % Plot velocity components
    ax = set_fixed_figure(config.id_polar3, config);

    plot(ax, mix2.polar.ux, mix2.polar.uy, 'k', 'LineWidth', config.linewidth, 'LineStyle', config.linestyle);
    xlabel(ax, '$u_{2x}$ [m/s]', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    ylabel(ax, '$u_{2y}$ [m/s]', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    %     ylim(ax, [0, max(mix2.polar.uy)]);
end
