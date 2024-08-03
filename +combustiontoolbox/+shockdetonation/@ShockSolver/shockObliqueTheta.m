function [mix1, mix2_1, mix2_2] = shockObliqueTheta(obj, mix1, u1, theta, varargin)
    % Compute pre-shock and post-shock states of an oblique shock wave
    % given the deflection angle.
    %
    % Two solutions:
    %   * Weak shock
    %   * Strong shock
    %
    % Args:
    %     obj (ShockSolver): ShockSolver object
    %     mix1 (Mixture):    Properties of the mixture in the pre-shock state
    %     u1 (float):        Pre-shock velocity [m/s]
    %     theta (float):     Deflection angle [deg]
    %
    % Optional Args:
    %     * mix2_1 (Mixture): Properties of the mixture in the post-shock state - weak shock (previous calculation)
    %     * mix2_2 (Mixture): Properties of the mixture in the post-shock state - strong shock (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (Mixture): Properties of the mixture in the pre-shock state
    %     * mix2_1 (Mixture): Properties of the mixture in the post-shock state - weak shock
    %     * mix2_2 (Mixture): Properties of the mixture in the post-shock state - strong shock
    %
    % Examples:
    %     * [mix1, mix2_1, mix2_2] = shockObliqueTheta(ShockSolver(), mix1, u1, theta)
    %     * [mix1, mix2_1, mix2_2] = shockObliqueTheta(ShockSolver(), mix1, u1, theta, mix2_1, mix2_2)

    % Unpack input data
    mix1 = unpack(mix1, u1, theta);

    % Definitions
    a1 = mix1.sound; % [m/s]
    u1 = mix1.u; % [m/s]
    M1 = u1 / a1; % [-]
    theta = mix1.theta * pi / 180; % [rad]
    betaMin = asin(1 / M1); % [rad]
    betaMax = pi / 2; % [rad]
    lambda = 0.7;
    m = 4;
    
    % Create a temporal shallow copy of mix1 to avoid overwrite results
    temp_mix1 = mix1.copy;

    % Solve first branch  (weak shock)
    betaGuess = 0.5 * (betaMin + betaMax); % [rad]
    % beta_guess = compute_guess_beta(M1, mix1.gamma, theta, true);
    if nargin > 4
        mix2_1 = solve_shock_oblique(temp_mix1, betaGuess, varargin{1});
    else
        mix2_1 = solve_shock_oblique(temp_mix1, betaGuess);
    end

    % Solve second branch (strong shock)
    if nargout > 2
        betaGuess = betaMax * 0.97; % [rad]
        % beta_guess = compute_guess_beta(M1, mix1.gamma, theta, false);
        if nargin > 4
            mix2_2 = solve_shock_oblique(temp_mix1, betaGuess, varargin{2});
        else
            mix2_2 = solve_shock_oblique(temp_mix1, betaGuess);
        end

    else
        mix2_2 = [];
    end
    
    % NESTED FUNCTIONS
    function mix2 = solve_shock_oblique(mix1, betaGuess, varargin)
        % Obtain wave angle, pre-shock state and post-shock states for an oblique shock
        STOP = 1; it = 0; itMax = obj.itOblique;

        while STOP > obj.tolOblique && it < itMax
            % Update iteration
            it = it + 1;
            % Compute incident normal velocity
            u1n = u1 * sin(betaGuess); % [m/s]
            % Obtain post-shock state
            [~, mix2] = shockIncident(obj, mix1, u1n, varargin{:});
            u2n = mix2.uShock;
            % Compute f0 , df0, and d2f0
            f0 = f0_beta(betaGuess, theta, u2n, u1);
            df0 = df0_beta(betaGuess, u2n, u1);

            % Apply correction
            beta = betaGuess - f0 / df0;
            df0_2 = df0_beta(beta, u2n, u1);
            beta = betaGuess - 2 * f0 / (df0 + df0_2);
            % d2f0 = d2f0_beta(betaGuess, u2n, u1);
            % beta = betaGuess - (f0 * df0) / (df0^2 - m^-1 * f0 * d2f0);
            % Check value
            if beta < betaMin || beta > betaMax
                beta = betaGuess - lambda * f0 / df0;
                beta = max(1.0001 * betaMin, 0.9999 * min(beta, betaMax));
                % beta = betaGuess - lambda * (f0 * df0) / (df0^2 - 0.5 * f0 * d2f0);
            end

            % Compute error
            STOP = max(abs(beta - betaGuess) / abs(beta), abs(f0));
            % Update guess
            betaGuess = beta;
            % Debug
            % aux_lambda(it) = betaGuess * 180/pi;
            % aux_STOP(it) = STOP;
        end

        % debug_plot_error(it, aux_STOP, aux_lambda);
        a2 = mix2.sound;
        u2 = u2n * csc(beta - theta);
        %M2 = u2 / a2;

        % if beta < betaMin || beta > pi/2 || theta < 0 || M2 >= M1
        %     error('There is not solution for the given deflection angle %.2g', mix1.theta);
        % end
        % Check convergence
        combustiontoolbox.utils.printConvergence(it, obj.itOblique, mix2.T, STOP, obj.itOblique);

        % Save results
        mix2.beta = beta * 180 / pi; % [deg]
        mix2.theta = theta * 180 / pi; % [deg]
        mix2.betaMin = betaMin * 180 / pi; % [deg]
        mix2.betaMax = betaMax * 180 / pi; % [deg]
        mix2.sound = a2; % [m/s]
        mix2.u = u2; % [m/s]
        mix2.mach = u2 / a2; % [m/s]
        mix2.uNormal = u2n; % [m/s]
        mix2.uShock = u2; % [m/s]
    end

end

% SUB-PASS FUNCTIONS
function [mix1, mix2_1, mix2_2] = unpack(mix1, u1, theta, varargin)
    % Unpack input data
    mix1.u = u1; % pre-shock velocity [m/s] - laboratory fixed
    mix1.uShock = u1; % pre-shock velocity [m/s] - shock fixed
    mix1.mach = mix1.u / mix1.sound; % pre-shock Mach number [-]
    mix1.theta = theta; % deflection angle  [deg]

    if nargin > 3
        mix2_1 = varargin{1};
        mix2_2 = varargin{2};
    else
        mix2_1 = [];
        mix2_2 = [];
    end

end

function value = f0_beta(beta, theta, u2n, u1)
    % Function to find the roots of beta
    value = theta + atan(u2n / (u1 .* cos(beta))) - beta;
end

function value = df0_beta(beta, u2n, u1)
    % First derivative of the function to find the roots of beta
    value = (u2n * u1 * sin(beta)) / (u2n^2 + u1^2 * cos(beta)^2) - 1;
end

function value = d2f0_beta(beta, u2n, u1)
    % Second derivative of the function to find the roots of beta
    value = 0.5 * (u1 * u2n * cos(beta) * (3 * u1^2 + 2 * u2n^2 - u1^2 * cos(2 * beta))) / (u2n^2 + u1^2 * cos(beta)^2)^2;
end

function beta = compute_guess_beta(M1, gamma, theta, FLAG_WEAK)
    % Compute guess of wave angle [rad]
    l = sqrt((M1^2 -1)^2 - 3 * (1 + (gamma - 1) / 2 * M1^2) * (1 + (gamma + 1) / 2 * M1^2) * tan(theta)^2);
    m = ((M1^2 - 1)^3 - 9 * (1 + (gamma - 1) / 2 * M1^2) * (1 + (gamma - 1) / 2 * M1^2 + (gamma + 1) / 4 * M1^4) * tan(theta)^2) / l^3;
    beta = atan(((M1^2 - 1) + 2 * l * cos((4 * pi * FLAG_WEAK + acos(m)) / 3)) / (3 * (1 + (gamma - 1) / 2 * M1^2) * tan(theta)));
end
