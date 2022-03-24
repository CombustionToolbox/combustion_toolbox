function [mix1, mix2] = shock_polar(varargin)
    % Solve oblique shock

    % Unpack input data
    [self, mix1, mix2] = unpack(varargin);
    % Abbreviations
    TN = self.TN;
    % Definitions
    a1 = soundspeed(mix1);
    u1 = mix1.u;

    beta_min = asin(a1 / u1);
    beta = linspace(beta_min, pi/2, TN.N_points_polar);
    u1n = u1 * sin(beta);

%     % Compute as SDToolbox for comparison
%     step = 5; % [m/s]
%     u1n = linspace(a1 + 1, u1, (u1 - a1) / step);
%     beta = asin(u1n ./ u1);
%     beta_min = asin(a1 / u1);

    N = length(u1n);
    ut = u1 .* cos(beta);
    % Loop
    for i = N:-1:1
        [~, mix2] = shock_incident(self, mix1, u1n(i), mix2);
        a2(i) = soundspeed(mix2);
        u2n(i) = mix2.v_shock;
        p2(i) = pressure(mix2);
        theta(i) = beta(i) - atan(u2n(i) / ut(i));
        u2(i) = u2n(i) * csc(beta(i) - theta(i));
        Xi(:, i) = mix2.Xi;
    end
    % Velocity components downstream - laboratory fixed
    u2x = u2 .* cos(theta); % [m/s]
    u2y = u2 .* sin(theta); % [m/s]
    % Sonic point
    M2 = u2 ./ a2;
    ind_sonic = find(M2 > 1, 1, 'last');
    theta_sonic = theta(ind_sonic) * 180/pi;
    % Max deflection angle
    [theta_max, ind_max] = max(theta * 180/pi); % [deg]
    % Save results
    mix2.beta = beta(end) * 180/pi;             % [deg]
    mix2.theta = theta(end) * 180/pi;           % [deg]
    mix2.theta_max = theta_max;                 % [deg]
    mix2.beta_min = beta_min * 180/pi;          % [deg]
    mix2.beta_max = beta(ind_max) * 180/pi;     % [deg]
    mix2.theta_sonic = theta_sonic;             % [deg]
    mix2.beta_sonic = beta(ind_sonic) * 180/pi; % [deg]
    mix2.ind_max = ind_max;
    mix2.ind_sonic = ind_sonic;
    % Range values for the shock polar
    mix2.polar.p = p2;    % [bar]
    mix2.polar.Mach = M2; % [-]
    mix2.polar.un = u2n;  % [m/s]
    mix2.polar.ux = u2x;  % [m/s]
    mix2.polar.uy = u2y;  % [m/s]
    mix2.polar.beta = beta * 180/pi;   % [deg]
    mix2.polar.theta = theta * 180/pi; % [deg]
    mix2.polar.Xi = Xi;   % [-]
end

% SUB-PASS FUNCTIONS
function [self, str1, str2] = unpack(x)
    % Unpack inpv data
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