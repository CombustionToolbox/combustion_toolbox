function [mix1, mix2_1, mix2_2] = det_oblique_theta(self, mix1, drive_factor, theta, varargin)
    % Compute pre-shock and post-shock states of an oblique detonation wave
    % given the deflection angle.
    %
    % Two solutions:
    %   * Weak detonation
    %   * Strong detonation
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the mixture in the pre-shock state
    %     drive_factor (float): Drive factor [-]
    %     theta (float): Deflection angle [deg]
    %
    % Optional Args:
    %     - mix2_1 (struct): Properties of the mixture in the post-shock state - weak detonation (previous calculation)
    %     - mix2_2 (struct): Properties of the mixture in the post-shock state - strong detonation (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     - mix1 (struct): Properties of the mixture in the pre-shock state
    %     - mix2_1 (struct): Properties of the mixture in the post-shock state - weak detonation
    %     - mix2_2 (struct): Properties of the mixture in the post-shock state - strong detonation

    % Unpack input data
    [self, mix1, mix2_1, mix2_2] = unpack(self, mix1, drive_factor, theta, varargin);
    % Abbreviations
    TN = self.TN;
    % Compute CJ state
    [mix1, ~] = det_cj(self, mix1);
    mix1.cj_speed = mix1.u;
    % Set initial velocity for the given drive_factor factor
    mix1.u = mix1.u * drive_factor; % [m/s]
    % Definitions
    a1 = soundspeed(mix1); % [m/s]
    u1 = mix1.u; % [m/s]
    M1 = u1 / a1; % [-]
    theta = mix1.theta * pi / 180; % [rad]
    beta_min = asin(1 / M1); % [rad]
    beta_max = pi / 2; % [rad]
    lambda = 0.7;
    m = 4;
    % Solve first branch  (weak shock)
    beta_guess = 0.5 * (beta_min + beta_max); % [rad]
    mix2_1 = solve_det_oblique(mix2_1);
    % Solve second branch (strong shock)
    beta_guess = beta_max * 0.97; % [rad]
    mix2_2 = solve_det_oblique(mix2_2);
    % NESTED FUNCTIONS
    function mix2 = solve_det_oblique(mix2)
        % Obtain wave angle, pre-shock state and post-shock states for an
        % overdriven oblique detonation
        STOP = 1; it = 0; itMax = TN.it_oblique;

        while STOP > TN.tol_oblique && it < itMax
            % Update iteration
            it = it + 1;
            % Compute incident normal velocity
            drive_factor_n = drive_factor * sin(beta_guess); % [-]
            % Obtain post-shock state
            [~, mix2] = det_overdriven(self, mix1, drive_factor_n, mix2);
            u2n = mix2.v_shock;
            % Compute f0 , df0, and d2f0
            f0 = f0_beta(beta_guess, theta, u2n, u1);
            df0 = df0_beta(beta_guess, u2n, u1);

            % Apply correction
            % beta = beta_guess - f0 / df0;

            d2f0 = d2f0_beta(beta_guess, u2n, u1);
            beta = beta_guess - (f0 * df0) / (df0^2 - m^-1 * f0 * d2f0);
            % Check value
            if beta < beta_min || beta > beta_max
                % beta = beta_guess - lambda * f0 / df0;
                beta = beta_guess - lambda * (f0 * df0) / (df0^2 - 0.5 * f0 * d2f0);
            end

            % Compute error
            STOP = max(abs((beta - beta_guess) / beta), abs(f0));
            % Update guess
            beta_guess = beta;
        end

        a2 = mix2.sound;
        u2 = u2n * csc(beta - theta);
        M2 = u2 / a2;

        if beta < beta_min || beta > pi / 2 || theta < 0 || M2 >= M1
            error('There is not solution for the given deflection angle %.2g', mix1.theta);
        end

        if STOP > TN.tol_oblique
            print_error_root(it, itMax, mix2.T, STOP);
        end

        % Save results
        mix2.beta = beta * 180 / pi; % [deg]
        mix2.theta = theta * 180 / pi; % [deg]
        mix2.beta_min = beta_min * 180 / pi; % [deg]
        mix2.beta_max = beta_max * 180 / pi; % [deg]
        mix2.sound = a2; % [m/s]
        mix2.u = u2; % [m/s]
        mix2.un = u2n; % [m/s]
        mix2.v_shock = u2; % [m/s]
    end

end

% SUB-PASS FUNCTIONS
function [self, mix1, mix2_1, mix2_2] = unpack(self, mix1, drive_factor, theta, x)
    % Unpack input data
    mix1.drive_factor = drive_factor;
    mix1.theta = theta; % deflection angle  [deg]

    if length(x) > 1
        mix2_1 = x{1};
        mix2_2 = x{2};
    else
        mix2_1 = [];
        mix2_2 = [];
    end

    if isempty(mix2_2)
        mix1.cj_speed = [];
    else
        mix1.cj_speed = mix2_1.cj_speed;
    end

end

function value = f0_beta(beta, theta, u2n, u1)
    % Function to find the roots of beta
    value = theta + atan(u2n / (u1 .* cos(beta))) - beta;
end

function value = df0_beta(beta, u2n, u1)
    % Derivative of the function to find the roots of beta
    value = (u2n * u1 * sin(beta)) / (u2n^2 + u1^2 * cos(beta)^2) - 1;
end

function value = d2f0_beta(beta, u2n, u1)
    % Derivative of the function to find the roots of beta
    value = 0.5 * (u1 * u2n * cos(beta) * (3 * u1^2 + 2 * u2n^2 - u1^2 * cos(2 * beta))) / (u2n^2 + u1^2 * cos(beta)^2)^2;
end
