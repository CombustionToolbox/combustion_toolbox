function [str1, str2] = shock_oblique(varargin)
    % Solve oblique shock

    % Unpack input data
    [self, str1, str2] = unpack(varargin);
    % Abbreviations
    config = self.Misc.config;
    % Definitions
    a1 = soundspeed(str1);
    u1 = str1.u;
    M1 = u1/a1;
    step = 1;
    w1 = linspace(a1, u1, round((u1 - a1) / step));
    N = length(w1);
    beta = asin(w1 ./ u1);
    beta_min = asin(a1 / u1);
    v = u1 .* cos(beta);

    % Solve first case for initialization
    [str1, str2] = shock_incident(self, str1, w1(end), str2);
    % Loop
    for i = N:-1:1
        [str1, str2] = shock_incident(self, str1, w1(i), str2);
        a2(i) = soundspeed(str2);
        w2(i) = str2.v_shock;
        p2(i) = pressure(str2);
        theta(i) = beta(i) - atan(w2(i) / v(i));
        u2(i) = w2(i) * csc(beta(i) - theta(i));
    end

    % Velocity components downstream - laboratory fixed
    u2x = u2 .* cos(theta);
    u2y = u2 .* sin(theta);
    % Sonic point
    M2 = u2 ./ a2;
    ind_sonic = find(M2 > 1, 1, 'last');
    theta_sonic = theta(ind_sonic) * 180/pi;
    % Max deflection angle
    [theta_max, ind_max] = max(theta * 180/pi); % [deg]
    % Save results
    str1.beta = beta(end) * 180/pi;   % [deg]
    str1.theta = theta(end) * 180/pi; % [deg]
    str1.theta_max = theta_max;       % [deg]
    str1.theta_sonic = theta_sonic;   % [deg]
    % Plot (pressure, deflection) - shock polar
	figure(1); ax = gca;
    set(ax, 'LineWidth', config.linewidth, 'FontSize', config.fontsize-2, 'BoxStyle', 'full')
    grid(ax, 'off'); box(ax, 'off'); hold(ax, 'on'); ax.Layer = 'Top';

	plot(ax, [-flip(theta), theta] * 180/pi, [flip(p2), p2], 'k', 'LineWidth', config.linewidth);
    plot(ax, theta_max, p2(ind_max), 'ko', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    plot(ax, theta(ind_sonic) * 180/pi, p2(ind_sonic), 'ks', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    text(ax, 0, max(p2), strcat('$\mathcal{M}_1 = ', sprintf('%.2f', M1), '$'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', config.fontsize-8, 'Interpreter', 'latex')
%     text(ax, theta_max, p2(ind_max), 'detachment point \rightarrow m', 'HorizontalAlignment', 'right', 'FontSize', config.fontsize-8)
%     text(ax, theta_sonic * 180/pi, p2(ind_sonic), 'sonic point \rightarrow s', 'HorizontalAlignment', 'right', 'FontSize', config.fontsize-8)
    xlabel(ax, 'Deflection angle $[^\circ]$', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    ylabel(ax, 'Pressure [bar]','FontSize', config.fontsize, 'Interpreter', 'latex');
    ylim(ax, [1, max(p2)]);
    % Plot (wave angle, deflection) - shock polar
    figure(2); ax = gca;
    set(ax, 'LineWidth', config.linewidth, 'FontSize', config.fontsize-2, 'BoxStyle', 'full')
    grid(ax, 'off'); box(ax, 'off'); hold(ax, 'on'); ax.Layer = 'Top';

    plot(ax, [-flip(theta), theta] * 180/pi, [flip(beta), beta] * 180/pi, 'LineWidth', config.linewidth);  
    plot(ax, theta_max, beta(ind_max) * 180/pi, 'ko', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    plot(ax, theta(ind_sonic) * 180/pi, beta(ind_sonic) * 180/pi, 'ks', 'LineWidth', config.linewidth, 'MarkerFaceColor', 'auto');
    xlabel(ax, 'Deflection angle $[^\circ]$','FontSize', config.fontsize, 'Interpreter', 'latex');
    ylabel(ax, 'Wave angle $[^\circ]$','FontSize', config.fontsize, 'Interpreter', 'latex');
    % Plot velocity components
    figure(3); ax = gca;
    set(ax, 'LineWidth', config.linewidth, 'FontSize', config.fontsize-2, 'BoxStyle', 'full')
    grid(ax, 'off'); box(ax, 'off'); hold(ax, 'on'); ax.Layer = 'Top';

    plot(ax, u2x, u2y, 'LineWidth', config.linewidth);  
    xlabel(ax, '$u_{2x}$ [m/s]', 'FontSize', config.fontsize, 'Interpreter', 'latex');
    ylabel(ax, '$u_{2y}$ [m/s]', 'FontSize', config.fontsize, 'Interpreter', 'latex');
end

% SUB-PASS FUNCTIONS
function [self, str1, str2] = unpack(x)
    % Unpack input data
    self = x{1};
    str1 = x{2};
    u1   = x{3};
    str1.u = u1;       % velocity preshock [m/s] - laboratory fixed
    str1.v_shock = u1; % velocity preshock [m/s] - shock fixed
    if length(x) > 3
        str2 = x{4};
    else
        str2 = [];
    end
end