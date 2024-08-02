function [mix1, mix2] = detonationPolar(obj, mix1, driveFactor, varargin)
    % Compute detonation polar diagrams
    %
    % Args:
    %     obj (DetonationSolver): DetonationSolver object
    %     mix1 (Mixture): Properties of the mixture in the pre-shock state
    %     driveFactor (float): Driven factor
    %
    % Optional Args:
    %     mix2 (Mixture): Properties of the mixture in the post-shock state (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (Mixture): Properties of the mixture in the pre-shock state
    %     * mix2 (Mixture): Properties of the mixture at the post-shock state with the shock polar results
    %
    % Examples:
    %     * [mix1, mix2] = detonationPolar(DetonationSolver(), mix1, 3000)
    %     * [mix1, mix2] = detonationPolar(DetonationSolver(), mix1, 3000, mix2)

    % Unpack input data
    [mix1, mix2] = unpack(mix1, driveFactor, varargin);
    
    % Compute CJ state
    [mix1, ~] = detonationCJ(obj, mix1);
    mix1.cjSpeed = mix1.u;
    
    % Set initial velocity for the given drive factor
    mix1.u = mix1.u * driveFactor;   % pre-shock velocity [m/s] - laboratory fixed
    mix1.uShock = mix1.u;            % pre-shock velocity [m/s] - shock fixed
    mix1.mach = mix1.u / mix1.sound; % pre-shock Mach number [-]

    % Definitions
    u1 = mix1.u; % [m/s]
    betaMin = asin(mix1.cjSpeed / u1); % [rad]
    betaMax = pi / 2; % [rad]
    betaOver = linspace(betaMax, betaMin, round(obj.numPointsPolar / 2) ); % [rad]
    betaUnder = fliplr(betaOver); % [rad]

    % Compute overdriven branch
    [results_over, mix2] = compute_detonation(obj, @detonationOverdriven, mix1, betaOver, mix2, false);
    
    % Compute underdriven branch
    results_under = compute_detonation(obj, @detonationUnderdriven, mix1, betaUnder, mix1.copy, true);
    
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
    indexSonic = round(obj.numPointsPolar / 2);
    thetaSonic = theta(indexSonic) * 180 / pi; % [deg]

    % Max deflection angle
    [thetaMax, indexMax] = max(theta * 180 / pi); % [deg]

    % Save results
    mix2.beta = beta(end) * 180 / pi; % [deg]
    mix2.theta = theta(end) * 180 / pi; % [deg]
    mix2.thetaMin = thetaSonic; % [deg]
    mix2.thetaMax = thetaMax; % [deg]
    mix2.thetaSonic = thetaSonic; % [deg]
    mix2.betaMin = betaMin * 180 / pi; % [deg]
    mix2.betaMax = beta(indexMax) * 180 / pi; % [deg]
    mix2.betaSonic = beta(indexSonic) * 180 / pi; % [deg]
    mix2.indexMin = indexSonic;
    mix2.indexMax = indexMax;
    mix2.indexSonic = indexSonic;

    % Save results polar curves
    mix2.polar.p = p2; % [bar]
    mix2.polar.Mach = M2; % [-]
    mix2.polar.sound = a2; % [m/s]
    mix2.polar.u = u2; % [m/s]
    mix2.polar.uShock = u2n; % [m/s]
    mix2.polar.ux = u2x; % [m/s]
    mix2.polar.uy = u2y; % [m/s]
    mix2.polar.beta = beta * 180 / pi; % [deg]
    mix2.polar.theta = theta * 180 / pi; % [deg]
    mix2.polar.Xi = results.Xi; % [-]
    mix2.polar.T = results.T; % [K]
end

% SUB-PASS FUNCTIONS
function [mix1, mix2] = unpack(mix1, driveFactor, x)
    % Unpack input data
    mix1.driveFactor = driveFactor;
    
    if ~isempty(x)
        mix2 = x{1};
    else
        mix2 = [];
    end

    if isempty(mix2)
        mix1.cjSpeed = [];
    else
        mix1.cjSpeed = mix2.cjSpeed;
    end

end

function [results, mix2] = compute_detonation(obj, fun_det, mix1, beta, mix2, FLAG_UNDER)
    % Compute oblique detonation for a range of deflection angles (beta)
    % and a given driveFactor = M1 / M_cj

    % Definitions
    N = length(beta);
    driveFactor_n = mix1.driveFactor * sin(beta); % [-]
    ut = mix1.u .* cos(beta); % [m/s]

    % Initialization
    index = [];
    u2n = zeros(1, N);
    p2 = zeros(1, N);
    theta = zeros(1, N);
    a2 = zeros(1, N);
    u2 = zeros(1, N);
    T = zeros(1, N);
    Xi = zeros(mix1.chemicalSystem.numSpecies, N);
    
    % Create a temporal shallow copy of mix1 to avoid overwrite results
    temp_mix1 = mix1.copy;

    % Loop
    for i = 1:N
        [~, mix2] = fun_det(obj, temp_mix1, driveFactor_n(i), mix2);
        u2n(i) = mix2.uShock; % [m/s]
        p2(i) = pressure(mix2); % [bar]
        theta(i) = beta(i) - atan(u2n(i) / ut(i)); % [rad]
        a2(i) = soundspeed(mix2); % [m/s]
        u2(i) = u2n(i) * csc(beta(i) - theta(i)); % [m/s]
        T(i) = mix2.T; % [K]
        Xi(:, i) = mix2.Xi; % [-]

        % Get index spurious results
        if (FLAG_UNDER && round(u2n(i) / a2(i), 1) < 1) || (~FLAG_UNDER && round(u2n(i) / a2(i), 1) > 1)
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
