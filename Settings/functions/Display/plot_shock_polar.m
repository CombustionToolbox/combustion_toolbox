function [ax1, ax2, ax3] = plot_shock_polar(varargin)
    % Routine to obtain shock polar plots
    %   * Plot (pressure, deflection) 
    %   * Plot (wave angle, deflection)
    %   * Plot velocity components
    
    self = varargin{1};
    mix1 = varargin{2};
    mix2 = varargin{3};
    if nargin > 3
        mix2_case = varargin{4};
        mix0 = varargin{5};
    else
        mix2_case = cell(1, length(mix1));
        mix0 = cell(1, length(mix1));
    end
    % Abbreviations
    config = self.Misc.config;
    % Definitions
    N = length(mix1);
    % Plots
    for i = N:-1:1
        % Plot (pressure, deflection) - shock polar
	    ax1 = plot_shock_polar_pressure(mix1{i}, mix2{i}, config, mix2_case{i}, mix0{i});
        % Plot (wave angle, deflection) - shock polar
        ax2 = plot_shock_polar_wave(mix1{i}, mix2{i}, config);
        % Plot velocity components
        ax3 = plot_shock_polar_velocities(mix1{i}, mix2{i}, config);
    end
end

% SUB-PASS FUNCTIONS
function ax = set_fixed_figure(fixed_number, config)
    % Generate a figure with a defined object identifier (fixed_number)
    fig = figure(fixed_number);
    ax = gca;
    set(fig, 'position', [1921 -471 1080 1795]);
    set(ax, 'LineWidth', config.linewidth, 'FontSize', config.fontsize-2, 'BoxStyle', 'full');
    grid(ax, 'off'); box(ax, 'off'); hold(ax, 'on'); ax.Layer = 'Top'; axis(ax, 'tight');
end

function ax = plot_shock_polar_pressure(mix1, mix2, config, mix2_case, mix0)
    % Plot (pressure, deflection) - shock polar
    ax = set_fixed_figure(1, config);
    
    M1 = velocity_relative(mix1) / soundspeed(mix1);
    p1 = pressure(mix1);
    Mach_label = '$\mathcal{M}_1 = ';
    xaxis = 0;
    if ~isempty(mix2_case)
        M1 = velocity_relative(mix2_case) / soundspeed(mix2_case);
        p1 = mix0.p;
        Mach_label = '$\mathcal{M}_2 = ';
        xaxis = mix2_case.theta;
        y_lim = get(ax, 'YLim'); y_lim(2) = max(mix2.polar.p) / p1;
        plot(ax, [0, 0], y_lim, '-.k', 'LineWidth', config.linewidth)
    end
    plot(ax, xaxis + [- flip(mix2.polar.theta), mix2.polar.theta], [flip(mix2.polar.p), mix2.polar.p] / p1, 'k', 'LineWidth', config.linewidth);
%     plot(ax, theta_max, mix2.polar.p(mix2.ind_max), 'ko', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
%     plot(ax, mix2.polar.theta(mix2.ind_sonic), mix2.polar.p(mix2.ind_sonic), 'ks', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    text(ax, xaxis, max(mix2.polar.p), strcat(Mach_label, sprintf('%.2g', M1), '$'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', config.fontsize-4, 'Interpreter', 'latex');
%     text(ax, theta_max, str2.polar.p2(str2.ind_max), 'detachment point \rightarrow m', 'HorizontalAlignment', 'right', 'FontSize', config.fontsize-6)
%     text(ax, theta_sonic, str2.polar.p2(str2.ind_sonic), 'sonic point \rightarrow s', 'HorizontalAlignment', 'right', 'FontSize', config.fontsize-6)
    xlabel(ax, 'Deflection angle $[^\circ]$', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    ylabel(ax, '$p/p_1$','FontSize', config.fontsize, 'Interpreter', 'latex');
%     ylim(ax, [1, max(str2.polar.p2)]);
end

function ax = plot_shock_polar_wave(mix1, mix2, config)
    % Plot (wave angle, deflection) - shock polar
    ax = set_fixed_figure(2, config);
    
    M1 = velocity_relative(mix1) / soundspeed(mix1);
    plot(ax, mix2.polar.theta, mix2.polar.beta, 'k', 'LineWidth', config.linewidth);  
%     plot(ax, mix2.theta_max, mix2.polar.beta(str2.ind_max), 'ko', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
%     plot(ax, str2.polar.theta(str2.ind_sonic), mix2.polar.beta(str2.ind_sonic), 'ks', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
%     text(ax, mix2.theta_max, mix2.polar.beta(str2.ind_max),
%     strcat('$\mathcal{M}_1 = ', sprintf('%.2g', M1), '$'), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', config.fontsize-6, 'Interpreter', 'latex');
    xlabel(ax, 'Deflection angle $[^\circ]$','FontSize', config.fontsize, 'Interpreter', 'latex');
    ylabel(ax, 'Wave angle $[^\circ]$','FontSize', config.fontsize, 'Interpreter', 'latex');
%     xlim(ax, [0, mix2.theta_max]);
    ylim(ax, [0, 90.001]);
end

function ax = plot_shock_polar_velocities(mix1, mix2, config)
    % Plot velocity components
    ax = set_fixed_figure(3, config);

    plot(ax, mix2.polar.ux, mix2.polar.uy, 'k', 'LineWidth', config.linewidth);  
    xlabel(ax, '$u_{2x}$ [m/s]', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    ylabel(ax, '$u_{2y}$ [m/s]', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    %     ylim(ax, [0, max(mix2.polar.uy)]);
end