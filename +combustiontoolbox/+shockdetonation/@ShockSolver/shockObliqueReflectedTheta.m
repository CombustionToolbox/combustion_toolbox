function [mix1, mix2, mix5_1, mix5_2] = shockObliqueReflectedTheta(obj, mix1, u2, theta, mix2, varargin)
    % Compute pre-shock and post-shock states of an oblique reflected shock wave
    % given the deflection angle.
    %
    % Two solutions:
    %   * Weak shock
    %   * Strong shock
    %
    % Args:
    %     obj (ShockSolver): ShockSolver object
    %     mix1 (Mixture): Properties of the mixture in the pre-shock state of the incident shock
    %     u2 (float): Post-shock velocity [m/s] of the incident shock
    %     theta (float): Deflection angle [deg]
    %     mix2 (Mixture): Properties of the mixture in the post-shock state of the incident shock
    %
    % Optional Args:
    %     * mix5_1 (Mixture): Properties of the mixture in the post-shock state of the reflected shock - weak shock (previous calculation)
    %     * mix5_2 (Mixture): Properties of the mixture in the post-shock state of the reflected shock - strong shock (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (Mixture): Properties of the mixture in the pre-shock state of the incident shock
    %     * mix2 (Mixture): Properties of the mixture in the post-shock state of the incident shock
    %     * mix5_1 (Mixture): Properties of the mixture in the post-shock state of the reflected shock - weak shock
    %     * mix5_2 (Mixture): Properties of the mixture in the post-shock state of the reflected shock - strong shock
    %
    % Examples:
    %     * [mix1, mix2, mix5_1, mix5_2] = shockObliqueReflectedTheta(ShockSolver(), mix1, u2, theta, mix2)
    %     * [mix1, mix2, mix5_1, mix5_2] = shockObliqueReflectedTheta(ShockSolver(), mix1, u2, theta, mix2, mix5_1, mix5_2)

    % Unpack input data
    [obj, mix1, mix2, mix5_1, mix5_2] = unpack(obj, mix1, u2, theta, mix2, varargin);

    % Definitions
    a2 = mix2.sound; % [m/s]
    u2 = mix2.u; % [m/s]
    M2 = u2 / a2; % [-]
    theta = mix2.theta * pi / 180; % [rad]
    betaMin = asin(1 / M2); % [rad]
    betaMax = pi / 2; % [rad]
    
    % Miscellaneous
    temp = obj.FLAG_RESULTS;
    obj.FLAG_RESULTS = false;

    % Obtain maximum deflection angle
    [~, mix5_polar] = shockPolar(obj, mix2, u2);
    thetaMaxReflected = mix5_polar.thetaMax * pi / 180;
    
    % Miscellaneous
    obj.FLAG_RESULTS = temp;
    
    % Check if theta is in the range for a Regular Reflection (RR)
    if theta > thetaMaxReflected
        error(['There is not Regular Reflection (RR) solution for the given conditions:\n' ...
                '   * incident velocity [m/s] %.2f\n' ...
                '   * wave angle [deg]        %.2f\n'], mix1.u, mix2.beta);
    end

    % Solve first branch  (weak shock)
    betaGuess = 0.5 * (betaMin + betaMax); % [rad]
    mix5_1 = solve_shock_oblique(mix5_1);

    % Solve second branch (strong shock)
    betaGuess = betaMax * 0.95; % [rad]
    mix5_2 = solve_shock_oblique(mix5_2);

    % NESTED FUNCTIONS
    function mix5 = solve_shock_oblique(mix5)
        % Obtain wave angle, pre-shock state and post-shock states for an oblique shock
        STOP = 1; it = 0; itMax = obj.itOblique;

        while STOP > obj.tolOblique && it < itMax
            it = it + 1;

            u2n = u2 * sin(betaGuess); % [m/s]
            [~, mix5] = shockIncident(obj, mix2, u2n);
            u5n = mix5.uShock;
            % Compute f0 and df0
            f0 = f0_beta(betaGuess, theta, u5n, u2);
            df0 = df0_beta(betaGuess, u5n, u2);
            % Compute new value
            beta = betaGuess - f0 / df0;
            % Compute error
            STOP = max(abs(beta - betaGuess) / abs(beta), abs(f0));
            % Update guess
            betaGuess = beta;
        end

        a5 = mix5.sound * csc(beta - theta);
        u5 = u5n * csc(beta - theta);
        M5 = u5 / a5;

        if STOP > obj.tolOblique || beta < betaMin || beta > pi / 2 || theta < 0 || M5 >= M2
            error('There is not solution for the given deflection angle %.2g', mix5.theta);
        end

        % Save results
        mix5.beta = beta * 180 / pi; % [deg]
        mix5.theta = theta * 180 / pi; % [deg]
        mix5.betaMin = betaMin * 180 / pi; % [deg]
        mix5.betaMax = betaMax * 180 / pi; % [deg]
        mix5.sound = a5; % [m/s]
        mix5.u = u5; % [m/s]
        mix5.mach = u5 / a5; % [m/s]
        mix5.uNormal = u5n; % [m/s]
        mix5.uShock = u5; % [m/s]
    end

end

% SUB-PASS FUNCTIONS
function [self, mix1, mix2, mix5_1, mix5_2] = unpack(self, mix1, u2, theta, mix2, x)
    % Unpack input data
    mix2.u = u2; % pre-shock velocity [m/s] - laboratory fixed
    mix2.uShock = u2; % pre-shock velocity [m/s] - shock fixed
    mix2.mach = mix2.uShock / mix2.sound; % pre-shock Mach number [-]
    mix2.theta = theta; % deflection angle  [deg]

    if length(x) > 5
        mix5_1 = x{6};
        mix5_2 = x{7};
    else
        mix5_1 = [];
        mix5_2 = [];
    end

end

function value = f0_beta(beta, theta, w5, u2)
    % Function to find the roots of beta
    value = theta + atan(w5 / (u2 .* cos(beta))) - beta;
end

function value = df0_beta(beta, w5, u2)
    % Derivative of the function to find the roots of beta
    value = (w5 * u2 * sin(beta)) / (w5^2 + u2^2 * cos(beta)^2) - 1;
end
