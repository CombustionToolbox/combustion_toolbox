function [mix1, mix2] = shockObliqueBeta(obj, mix1, u1, beta, varargin)
    % Compute pre-shock and post-shock states of an oblique shock wave
    % given the wave angle (one solution)
    %
    % Args:
    %     obj (ShockSolver): ShockSolver object
    %     mix1 (Mixture):    Properties of the mixture in the pre-shock state
    %     u1 (float):        Pre-shock velocity [m/s]
    %     beta (float):      Wave angle [deg] of the incident oblique shock
    %
    % Optional Args:
    %     mix2 (Mixture): Properties of the mixture in the post-shock state (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (Mixture): Properties of the mixture in the pre-shock state
    %     * mix2 (Mixture): Properties of the mixture at the post-shock state
    %
    % Examples:
    %     * [mix1, mix2] = shockObliqueBeta(ShockSolver(), mix1, u1, beta)
    %     * [mix1, mix2] = shockObliqueBeta(ShockSolver(), mix1, u1, beta, mix2)


    % Unpack input data
    [mix1, mix2] = unpack(mix1, u1, beta, varargin{:});

    % Definitions
    a1 = soundspeed(mix1); % [m/s]
    u1 = mix1.u; % [m/s]
    M1 = u1 / a1; % [-]
    beta = mix1.beta * pi / 180; % [rad]
    betaMin = asin(1 / M1); % [rad]
    betaMax = pi / 2; % [rad]
    u1n = u1 * sin(beta); % [m/s]

    % Check range beta
    if beta < betaMin || beta > betaMax
        error([' \ nERROR !The given wave angle beta = %.2g is not in the '...
               'range of possible solutions [%.2g, %2.g]'], mix1.beta, ...
               betaMin * 180 / pi, betaMax * 180 / pi);
    end

    % Obtain deflection angle, pre-shock state and post-shock states for
    % an oblique shock
    [~, mix2] = shockIncident(obj, mix1, u1n, mix2);

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

% SUB-PASS FUNCTIONS
function [mix1, mix2] = unpack(mix1, u1, beta, varargin)
    % Unpack input data
    mix1.u = u1; % velocity preshock [m/s] - laboratory fixed
    mix1.uShock = u1; % velocity preshock [m/s] - shock fixed
    mix1.beta = beta; % wave angle  [deg]

    if nargin > 3
        mix2 = varargin{1}.copy();
        return
    end

    mix2 = mix1.copy();
end
