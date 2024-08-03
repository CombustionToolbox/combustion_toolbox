function [mix1, mix2_1, mix2_2] = detonationObliqueBeta(obj, mix1, driveFactor, beta, varargin)
    % Compute pre-shock and post-shock states of an oblique detonation wave
    % given the wave angle.
    %
    % Two solutions:
    %   * Underdriven detonation
    %   * Overdriven detonation
    %
    % Args:
    %     obj (DetonationSolver): Detonation solver object
    %     mix1 (struct): Properties of the mixture in the pre-shock state
    %     driveFactor (float): Driven factor [-]
    %     beta (float): Wave angle [deg] of the incident oblique detonation
    %
    % Optional Args:
    %     mix2 (Mixture): Properties of the mixture in the post-shock state (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (Mixture): Properties of the mixture in the pre-detonation state
    %     * mix2_1 (Mixture): Properties of the mixture in the post-detonation state - underdriven detonation
    %     * mix2_2 (Mixture): Properties of the mixture in the post-detonation state - overdriven detonation
    %
    % Examples:
    %     * [mix1, mix2_1, mix2_2] = detonationObliqueBeta(DetonationSolver(), mix1, 2, 60)
    %     * [mix1, mix2_1, mix2_2] = detonationObliqueBeta(DetonationSolver(), mix1, 2, 60, mix2)

    % Unpack input data
    [mix1, mix2_1, mix2_2] = unpack(mix1, driveFactor, beta, varargin);

    % Compute CJ state
    [mix1, ~] = detonationCJ(obj, mix1);
    mix1.cjSpeed = mix1.u;

    % Set initial velocity for the given drive factor
    mix1.u = mix1.u * driveFactor;   % pre-shock velocity [m/s] - laboratory fixed
    mix1.uShock = mix1.u;            % pre-shock velocity [m/s] - shock fixed
    mix1.mach = mix1.u / mix1.sound; % pre-shock Mach number [-]

    % Definitions
    u1 = mix1.u; % [m/s]
    beta = mix1.beta * pi / 180; % [rad]
    betaMin = asin(mix1.cjSpeed / u1); % [rad]
    betaMax = pi / 2; % [rad]
    driveFactor_n = mix1.driveFactor * sin(beta); % [-]
    
    % Check range beta
    if beta < betaMin || beta > betaMax
        error([' \ nERROR !The given wave angle beta = %.2g is not in the '...
            'range of possible solutions [%.2g, %2.g]'], mix1.beta, ...
            betaMin * 180 / pi, betaMax * 180 / pi);
    end

    % Create a temporal shallow copy of mix1 to avoid overwrite results
    temp_mix1 = mix1.copy;
    
    % Solve underdriven detonation
    mix2_1 = solve_det_oblique(@detonationUnderdriven, mix2_1);

    % Solve overdriven detonation
    mix2_2 = solve_det_oblique(@detonationOverdriven, mix2_2);

    % NESTED FUNCTIONS
    function mix2 = solve_det_oblique(det_fun, mix2)
        % Obtain deflection angle, pre-shock state and post-shock states for
        % an overdriven oblique detonation
        [~, mix2] = det_fun(obj, temp_mix1, driveFactor_n, mix2);
    
        u2n = mix2.uShock;
        theta = beta - atan(u2n / (u1 .* cos(beta)));
        a2 = mix2.sound;
        u2 = u2n * csc(beta - theta);
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
function [mix1, mix2_1, mix2_2] = unpack(mix1, driveFactor, beta, x)
    % Unpack input data
    mix1.driveFactor = driveFactor;
    mix1.beta = beta; % wave angle  [deg]

    if ~isempty(x)
        mix2_1 = x{1};
        mix2_2 = x{2};
    else
        mix2_1 = mix1.copy;
        mix2_2 = [];
    end

    if isempty(mix2_2)
        mix1.cjSpeed = [];
    else
        mix1.cjSpeed = mix2_1.cjSpeed;
    end

end