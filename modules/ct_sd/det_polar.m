function [mix1, mix2] = det_polar(self, mix1, drive_factor, varargin)
    % Compute detonation polars
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the mixture in the pre-shock state
    %     u1 (float):    Pre-shock velocity [m/s]
    %
    % Optional Args:
    %     mix2 (struct): Properties of the mixture in the post-shock state (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     - mix1 (struct): Properties of the mixture in the pre-shock state
    %     - mix2 (struct): Properties of the mixture at the post-shock state with the shock polar results

    % Unpack input data
    [self, mix1, mix2] = unpack(self, mix1, drive_factor, varargin);
    % Abbreviations
    TN = self.TN;
    % Compute CJ state
    [mix1, ~] = det_cj(self, mix1);
    mix1.cj_speed = mix1.u;
    % Set initial velocity for the given drive factor
    mix1.u = mix1.u * drive_factor; % [m/s]
    % Definitions
    u1 = mix1.u; % [m/s]
    beta_min = asin(mix1.cj_speed / u1); % [rad]
    beta_max = pi / 2; % [rad]
    beta_over = linspace(beta_max, beta_min, TN.N_points_polar / 2); % [rad]
    beta_under = fliplr(beta_over); % [rad]
    % Compute overdriven branch
    [results_over, mix2] = compute_detonation(self, @det_overdriven, mix1, beta_over, mix2, false);
    % Compute underdriven branch
    results_under = compute_detonation(self, @det_underdriven, mix1, beta_under, mix2, true);
    % Append results polar curves
    results = append_structs(results_over, results_under);
    % Get results polar curves
    u2n = results.u2n; % [m/s]
    p2 = results.p2; % [bar]
    theta = results.theta; % [rad]
    a2 = results.a2; % [m/s]
    u2 = results.u2; % [m/s]
    beta = results.beta; % [rad]
    M2 = u2 ./ a2; % [-]
    % Velocity components downstream - laboratory fixed
    u2x = u2 .* cos(theta); % [m/s]
    u2y = u2 .* sin(theta); % [m/s]
    % Sonic point - CJ point (minimum beta)
    ind_sonic = TN.N_points_polar / 2;
    theta_sonic = theta(ind_sonic) * 180 / pi; % [deg]
    % Max deflection angle
    [theta_max, ind_max] = max(theta * 180 / pi); % [deg]
    % Save results
    mix2.beta = beta(end) * 180 / pi; % [deg]
    mix2.theta = theta(end) * 180 / pi; % [deg]
    mix2.theta_min = theta_sonic; % [deg]
    mix2.theta_max = theta_max; % [deg]
    mix2.theta_sonic = theta_sonic; % [deg]
    mix2.beta_min = beta_min * 180 / pi; % [deg]
    mix2.beta_max = beta(ind_max) * 180 / pi; % [deg]
    mix2.beta_sonic = beta(ind_sonic) * 180 / pi; % [deg]
    mix2.ind_min = ind_sonic;
    mix2.ind_max = ind_max;
    mix2.ind_sonic = ind_sonic;
    % Save results polar curves
    mix2.polar.p = p2; % [bar]
    mix2.polar.Mach = M2; % [-]
    mix2.polar.sound = a2; % [m/s]
    mix2.polar.u = u2; % [m/s]
    mix2.polar.un = u2n; % [m/s]
    mix2.polar.ux = u2x; % [m/s]
    mix2.polar.uy = u2y; % [m/s]
    mix2.polar.beta = beta * 180 / pi; % [deg]
    mix2.polar.theta = theta * 180 / pi; % [deg]
    mix2.polar.Xi = results.Xi; % [-]
    mix2.polar.T = results.T; % [K]
end

% SUB-PASS FUNCTIONS
function [self, mix1, mix2] = unpack(self, mix1, drive_factor, x)
    % Unpack input data
    mix1.drive_factor = drive_factor;

    if ~isempty(x)
        mix2 = x{1};
    else
        mix2 = [];
    end

    if isempty(mix2)
        mix1.cj_speed = [];
    else
        mix1.cj_speed = mix2.cj_speed;
    end

end

function [results, mix2] = compute_detonation(self, fun_det, mix1, beta, mix2, FLAG_UNDER)
    % Compute oblique detonation for a range of deflection angles (beta)
    % and a given drive_factor = M1 / M_cj

    % Definitions
    N = length(beta);
    drive_factor_n = mix1.drive_factor * sin(beta); % [-]
    ut = mix1.u .* cos(beta); % [m/s]
    % Initialization
    index = [];
    u2n = zeros(1, N);
    p2 = zeros(1, N);
    theta = zeros(1, N);
    a2 = zeros(1, N);
    u2 = zeros(1, N);
    T = zeros(1, N);
    Xi = zeros(self.S.NS, N);
    % Loop
    for i = 1:N
        [~, mix2] = fun_det(self, mix1, drive_factor_n(i), mix2);
        u2n(i) = mix2.v_shock; % [m/s]
        p2(i) = pressure(mix2); % [bar]
        theta(i) = beta(i) - atan(u2n(i) / ut(i)); % [rad]
        a2(i) = soundspeed(mix2); % [m/s]
        u2(i) = u2n(i) * csc(beta(i) - theta(i)); % [m/s]
        T(i) = mix2.T; % [K]
        Xi(:, i) = mix2.Xi; % [-]
        % Get index spurious results
        if (FLAG_UNDER && u2n(i) / a2(i) < 1) || (~FLAG_UNDER && u2n(i) / a2(i) > 1)
            index = [index, i];
            mix2 = [];
        end

    end

    % Remove spurious results
    beta(index) = [];
    u2n(index) = [];
    p2(index) = [];
    theta(index) = [];
    a2(index) = [];
    u2(index) = [];
    Xi(:, index) = [];
    T(index) = [];
    % Save results polar curves
    results.u2n = u2n;
    results.p2 = p2;
    results.theta = theta;
    results.a2 = a2;
    results.u2 = u2;
    results.Xi = Xi;
    results.T = T;
    results.beta = beta;
end
