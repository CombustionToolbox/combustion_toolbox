function [str1, str2, str5, str5_2] = shock_polar_reflected(varargin)
    % Solve oblique shock

    % Unpack input data
    [self, str1, str2, str5] = unpack(varargin);
    % Abbreviations
    config = self.Misc.config;
    % Definitions
    a2 = soundspeed(str2);
    u2 = str2.u;
    M2 = u2/a2;
    step = 5; % [m/s]
    u2n = linspace(a2 + 1, u2, (u2 - a2) / step);
    N = length(u2n);
    beta = asin(u2n ./ u2);
    beta_min = asin(a2 / u2);
    theta = self.PD.theta.value * 180/pi;
    % Solve first case for initialization
    [~, ~, str5, str5_2] = shock_oblique_reflected_theta(self, str1, str2.u, str2.theta_range(end) * pi/180, str2, str5);
    % Loop
    for i = N:-1:1
        [~, ~, str5, str5_2] = shock_oblique_reflected_theta(self, str1, str2.u, str2.theta_range(end) * pi/180, str2, str5);
        a5(i) = soundspeed(str5);
        u5n(i) = str5.v_shock;
        p2(i) = pressure(str5);
        u5(i) = u5n(i) * csc(beta(i) - theta(i));
    end

    % Velocity components downstream - laboratory fixed
    u5x = u5 .* cos(theta);
    u5y = u5 .* sin(theta);
    % Sonic point
    M5 = u5 ./ a5;
    ind_sonic = find(M5 > 1, 1, 'last');
    theta_sonic = theta(ind_sonic) * 180/pi;
    % Max deflection angle
    [theta_max, ind_max] = max(theta * 180/pi); % [deg]
    % Save results
    str5.beta = beta(end) * 180/pi;             % [deg]
    str5.theta = theta(end) * 180/pi;           % [deg]
    str5.theta_max = theta_max;                 % [deg]
    str5.beta_min = beta_min * 180/pi;          % [deg]
    str5.beta_max = beta(ind_max) * 180/pi;     % [deg]
    str5.theta_sonic = theta_sonic;             % [deg]
    str5.beta_sonic = beta(ind_sonic) * 180/pi; % [deg]
    
    str5.un = u5n;
    % Plot (pressure, deflection) - shock polar
	figure(1); ax = gca;
    set(ax, 'LineWidth', config.linewidth, 'FontSize', config.fontsize-2, 'BoxStyle', 'full')
    grid(ax, 'off'); box(ax, 'off'); hold(ax, 'on'); ax.Layer = 'Top'; axis(ax, 'tight')

	plot(ax, [-flip(theta), theta] * 180/pi, [flip(p2), p2] / pressure(str2), 'k', 'LineWidth', config.linewidth);
%     plot(ax, theta_max, p2(ind_max), 'ko', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
%     plot(ax, theta(ind_sonic) * 180/pi, p2(ind_sonic), 'ks', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    text(ax, 0, max(p2), strcat('$\mathcal{M}_1 = ', sprintf('%.2g', M2), '$'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', config.fontsize-8, 'Interpreter', 'latex')
%     text(ax, theta_max, p2(ind_max), 'detachment point \rightarrow m', 'HorizontalAlignment', 'right', 'FontSize', config.fontsize-8)
%     text(ax, theta_sonic * 180/pi, p2(ind_sonic), 'sonic point \rightarrow s', 'HorizontalAlignment', 'right', 'FontSize', config.fontsize-8)
    xlabel(ax, 'Deflection angle $[^\circ]$', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    ylabel(ax, '$p_2/p_1$','FontSize', config.fontsize, 'Interpreter', 'latex');
%     ylim(ax, [1, max(p2)]);
    % Plot (wave angle, deflection) - shock polar
    figure(2); ax = gca;
    set(ax, 'LineWidth', config.linewidth, 'FontSize', config.fontsize-2, 'BoxStyle', 'full')
    grid(ax, 'off'); box(ax, 'off'); hold(ax, 'on'); ax.Layer = 'Top'; axis(ax, 'tight');

    plot(ax, theta * 180/pi, beta * 180/pi, 'k', 'LineWidth', config.linewidth);  
%     plot(ax, theta_max, beta(ind_max) * 180/pi, 'ko', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
%     plot(ax, theta(ind_sonic) * 180/pi, beta(ind_sonic) * 180/pi, 'ks', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
%     text(ax, theta_max, beta(ind_max) * 180/pi, strcat('$\mathcal{M}_1 = ', sprintf('%.2g', M1), '$'), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', config.fontsize-8, 'Interpreter', 'latex')
    xlabel(ax, 'Deflection angle $[^\circ]$','FontSize', config.fontsize, 'Interpreter', 'latex');
    ylabel(ax, 'Wave angle $[^\circ]$','FontSize', config.fontsize, 'Interpreter', 'latex');
%     xlim(ax, [0, theta_max]);
    ylim(ax, [0, 90]);
    % Plot velocity components
    figure(3); ax = gca;
    set(ax, 'LineWidth', config.linewidth, 'FontSize', config.fontsize-2, 'BoxStyle', 'full')
    grid(ax, 'off'); box(ax, 'off'); hold(ax, 'on'); ax.Layer = 'Top'; axis(ax, 'tight');

    plot(ax, u5x, u5y, 'LineWidth', config.linewidth);  
    xlabel(ax, '$u_{2x}$ [m/s]', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    ylabel(ax, '$u_{2y}$ [m/s]', 'FontSize', config.fontsize, 'Interpreter', 'latex');
%     ylim(ax, [0, max(u2y)]);
end

% SUB-PASS FUNCTIONS
function [self, str1, str2, str5] = unpack(x)
    % Unpack input data
    self  = x{1};
    str1  = x{2};
    u2    = x{3};
    str2  = x{4};
    str2.u = u2;        % velocity preshock [m/s] - laboratory fixed
    str2.v_shock = u2;  % velocity preshock [m/s] - shock fixed
    if length(x) > 4
        str5 = x{5};
    else
        str5 = [];
    end
end