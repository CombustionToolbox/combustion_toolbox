function [mix1, mix2_1, mix2_2] = detonationObliqueTheta(obj, mix1, driveFactor, theta, varargin)
    % Compute pre-shock and post-shock states of an oblique detonation wave
    % given the deflection angle (only overdriven solution).
    %
    % Two solutions:
    %   * Weak detonation
    %   * Strong detonation
    %
    % Args:
    %     obj (DetonationSolver): Detonation solver object
    %     mix1 (Mixture): Properties of the mixture in the pre-shock state
    %     driveFactor (float): Drive factor [-]
    %     theta (float): Deflection angle [deg]
    %
    % Optional Args:
    %     * mix2_1 (Mixture): Properties of the mixture in the post-shock state - weak detonation (previous calculation)
    %     * mix2_2 (Mixture): Properties of the mixture in the post-shock state - strong detonation (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (Mixture): Properties of the mixture in the pre-detonation state
    %     * mix2_1 (Mixture): Properties of the mixture in the post-detonation state - weak detonation
    %     * mix2_2 (Mixture): Properties of the mixture in the post-detonation state - strong detonation
    %
    % Examples:
    %     * [mix1, mix2_1, mix2_2] = detonationObliqueTheta(DetonationSolver(), mix1, 2, 30)
    %     * [mix1, mix2_1, mix2_2] = detonationObliqueTheta(DetonationSolver(), mix1, 2, 30, mix2)

    % Unpack input data
    [mix1, mix2_1, mix2_2] = unpack(mix1, driveFactor, theta, varargin);

    % Compute CJ state
    [mix1, ~] = detonationCJ(obj, mix1);
    mix1.cjSpeed = mix1.u;

    % Set initial velocity for the given driveFactor factor
    mix1.u = mix1.u * driveFactor;   % pre-shock velocity [m/s] - laboratory fixed
    mix1.uShock = mix1.u;            % pre-shock velocity [m/s] - shock fixed
    mix1.mach = mix1.u / mix1.sound; % pre-shock Mach number [-]

    % Definitions
    a1 = soundspeed(mix1); % [m/s]
    u1 = mix1.u; % [m/s]
    M1 = u1 / a1; % [-]
    theta = mix1.theta * pi / 180; % [rad]
    betaMin = asin(mix1.cjSpeed / u1); % [rad]
    betaMax = pi / 2; % [rad]
    lambda = 0.7;
    m = 4;

    % Solve first branch  (weak shock)
    beta_guess = 0.5 * (betaMin + betaMax); % [rad]
    mix2_1 = solve_det_oblique(mix2_1);

    % Solve second branch (strong shock)
    beta_guess = betaMax * 0.97; % [rad]
    mix2_2 = solve_det_oblique(mix2_2);
    
    % NESTED FUNCTIONS
    function mix2 = solve_det_oblique(mix2)
        % Obtain wave angle, pre-shock state and post-shock states for an
        % overdriven oblique detonation
        STOP = 1; it = 0; itMax = obj.itMax;

        while STOP > obj.tol0 && it < itMax
            % Update iteration
            it = it + 1;
            % Compute incident normal velocity
            driveFactor_n = driveFactor * sin(beta_guess); % [-]
            % Obtain post-shock state
            [~, mix2] = detonationOverdriven(obj, mix1, driveFactor_n, mix2);
            u2n = mix2.uShock;
            % Compute f0 , df0, and d2f0
            f0 = f0_beta(beta_guess, theta, u2n, u1);
            df0 = df0_beta(beta_guess, u2n, u1);

            % Apply correction
            % beta = beta_guess - f0 / df0;

            d2f0 = d2f0_beta(beta_guess, u2n, u1);
            beta = beta_guess - (f0 * df0) / (df0^2 - m^-1 * f0 * d2f0);
            % Check value
            if beta < betaMin || beta > betaMax
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

        if beta < betaMin || beta > pi / 2 || theta < 0 || M2 >= M1
            error('There is not solution for the given deflection angle %.2g', mix1.theta);
        end

        if STOP > obj.tol0
            print_error_root(it, itMax, mix2.T, STOP);
        end

        % Save results
        mix2.beta = beta * 180 / pi; % [deg]
        mix2.theta = theta * 180 / pi; % [deg]
        mix2.betaMin = betaMin * 180 / pi; % [deg]
        mix2.betaMax = betaMax * 180 / pi; % [deg]
        mix2.sound = a2; % [m/s]
        mix2.u = u2; % [m/s]
        mix2.uNormal = u2n; % [m/s]
        mix2.uShock = u2; % [m/s]
    end

end

% SUB-PASS FUNCTIONS
function [mix1, mix2_1, mix2_2] = unpack(mix1, driveFactor, theta, x)
    % Unpack input data
    mix1.driveFactor = driveFactor;
    mix1.theta = theta; % deflection angle  [deg]

    if length(x) > 1
        mix2_1 = x{1};
        mix2_2 = x{2};
    else
        mix2_1 = [];
        mix2_2 = [];
    end

    if isempty(mix2_2)
        mix1.cjSpeed = [];
    else
        mix1.cjSpeed = mix2_1.cjSpeed;
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
