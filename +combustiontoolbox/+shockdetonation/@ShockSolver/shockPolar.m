function [mix1, mix2] = shockPolar(obj, mix1, u1)
    % Compute shock polar diagrams
    %
    % Args:
    %     obj (ShockSolver): ShockSolver object
    %     mix1 (Mixture):    Properties of the mixture in the pre-shock state
    %     u1 (float):        Pre-shock velocity [m/s]
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (Mixture): Properties of the mixture in the pre-shock state
    %     * mix2 (Mixture): Properties of the mixture at the post-shock state with the shock polar results
    %
    % Example:
    %     * [mix1, mix2] = shockPolar(ShockSolver(), mix1, u1)

    % Unpack input data
    [mix1, mix2] = unpack(mix1, u1);
    
    % Definitions
    a1 = mix1.sound;
    u1 = mix1.u;

    betaMin = asin(a1 / u1);
    beta = linspace(betaMin, pi / 2, obj.numPointsPolar);
    u1n = u1 * sin(beta);

    % Compute as SDToolbox for comparison
    % step = 5; % [m/s]
    % u1n = linspace(a1 + 1, u1, (u1 - a1) / step);
    % beta = asin(u1n ./ u1);
    % betaMin = asin(a1 / u1);

    N = length(u1n);
    ut = u1 .* cos(beta);
    
    % Create a temporal shallow copy of mix1 to avoid overwrite results
    temp_mix1 = mix1.copy;

    % Loop
    for i = N:-1:2
        [~, mix2] = shockIncident(obj, temp_mix1, u1n(i), mix2);
        a2(i) = mix2.sound;
        u2n(i) = mix2.uShock;
        p2(i) = mix2.p;
        theta(i) = beta(i) - atan(u2n(i) / ut(i));
        u2(i) = u2n(i) * csc(beta(i) - theta(i));
        Xi(:, i) = mix2.Xi;

        % Print results
        if obj.FLAG_RESULTS
            temp_mix1.problemType = 'SHOCK_OBLIQUE';
            mix2.problemType = 'SHOCK_OBLIQUE';
            mix2.beta = beta(i) * 180 / pi;    % [deg]
            mix2.theta = theta(i) * 180 / pi;  % [deg]
            mix2.betaMin = betaMin * 180 / pi; % [deg]
            print(temp_mix1, mix2);
        end

    end

    % Calculate last case (beta = betaMin), which corresponds with the sonic condition
    i = 1;
    mix2 = mix1.copy;
    a2(i) = mix2.sound;
    u2n(i) = a2(i);
    p2(i) = mix2.p;
    theta(i) = beta(i) - atan(u2n(i) / ut(i));
    u2(i) = u2n(i) * csc(beta(i) - theta(i));
    Xi(:, i) = mix2.Xi;
    mix2.u = 0;
    mix2.mach = 0;
    mix2.uShock = mix2.sound;
    mix2.dVdT_p = 1;
    mix2.dVdp_T = -1;
    mix2.errorMoles = 0; mix2.errorMolesIons = 0; mix2.errorProblem = 0;
    
    % Velocity components downstream - laboratory fixed
    u2x = u2 .* cos(theta); % [m/s]
    u2y = u2 .* sin(theta); % [m/s]

    % Sonic point
    M2 = u2 ./ a2;
    indexSonic = find(M2 >= 1, 1, 'last');
    thetaSonic = theta(indexSonic) * 180 / pi;

    % Max deflection angle
    [thetaMax, indexMax] = max(theta * 180 / pi); % [deg]

    % Save results
    mix2.beta = beta(end) * 180 / pi; % [deg]
    mix2.theta = theta(end) * 180 / pi; % [deg]
    mix2.thetaMax = thetaMax; % [deg]
    mix2.betaMin = betaMin * 180 / pi; % [deg]
    mix2.betaMax = beta(indexMax) * 180 / pi; % [deg]
    mix2.thetaSonic = thetaSonic; % [deg]
    mix2.betaSonic = beta(indexSonic) * 180 / pi; % [deg]
    mix2.indexMax = indexMax;
    mix2.indexSonic = indexSonic;
    
    % Range values for the shock polar
    mix2.polar.p = p2; % [bar]
    mix2.polar.mach = M2; % [-]
    mix2.polar.u = u2; % [m/s]
    mix2.polar.uNormal = u2n; % [m/s]
    mix2.polar.ux = u2x; % [m/s]
    mix2.polar.uy = u2y; % [m/s]
    mix2.polar.beta = beta * 180 / pi; % [deg]
    mix2.polar.theta = theta * 180 / pi; % [deg]
    mix2.polar.Xi = Xi; % [-]
end

% SUB-PASS FUNCTIONS
function [mix1, mix2] = unpack(mix1, u1)
    % Unpack input data
    mix1.u = u1; % pre-shock velocity [m/s] - laboratory fixed
    mix1.uShock = u1; % pre-shock velocity [m/s] - shock fixed
    mix1.mach = mix1.u / mix1.sound; % pre-shock Mach number [-]
    mix2 = [];
end
