function [mix1, mix2] = det_oblique_beta(self, mix1, overdriven, beta, varargin)
    % Compute pre-shock and post-shock states of an oblique detonation wave 
    % given the wave angle (one solution) 
    %
    % Args:
    %     self (struct):      Data of the mixture, conditions, and databases
    %     mix1 (struct):      Properties of the mixture in the pre-shock state
    %     overdriven (float): Overdriven factor [-]
    %     beta (float):       Wave angle [deg] of the incident oblique detonation
    %
    % Optional Args:
    %     mix2 (struct):      Properties of the mixture in the post-shock state (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     - mix1 (struct):      Properties of the mixture in the pre-shock state
    %     - mix2 (struct):      Properties of the mixture at the post-shock state
    
    % Unpack input data
    [self, mix1, mix2] = unpack(self, mix1, overdriven, beta, varargin);
    % Compute CJ state
    [mix1, ~] = cj_detonation(self, mix1);
    mix1.cj_speed = mix1.u;
    % Definitions
    a1 = soundspeed(mix1);     % [m/s]
    u1 = mix1.u;               % [m/s]
    M1 = u1 / a1;              % [-]
    beta = mix1.beta * pi/180; % [rad]
    beta_min = asin(1/M1);     % [rad]
    beta_max = pi/2;           % [rad]
    overdriven_n = mix1.overdriven * sin(beta); % [-]
    % Check range beta
    if beta < beta_min || beta > beta_max
        error(['\nERROR! The given wave angle beta = %.2g is not in the '...
               'range of possible solutions [%.2g, %2.g]'], mix1.beta,...
               beta_min * 180/pi, beta_max * 180/pi);
    end
    
    % Obtain deflection angle, pre-shock state and post-shock states for
    % an oblique detonation
    [~, mix2] = overdriven_detonation(self, mix1, overdriven_n, mix2);
    
    u2n = mix2.v_shock;
    theta = beta - atan(u2n / (u1 .* cos(beta)));
    u2 = u2n * csc(beta - theta);
    % Save results
    mix2.beta = beta * 180/pi;   % [deg]
    mix2.theta = theta * 180/pi; % [deg]
    mix2.beta_min = beta_min * 180/pi; % [deg]
    mix2.beta_max = beta_max * 180/pi; % [deg]
    mix2.u = u2; % [m/s]
    mix2.u2n = u2n; % [m/s]
    mix2.v_shock = u2; % [m/s]
end

% SUB-PASS FUNCTIONS
function [self, mix1, mix2] = unpack(self, mix1, overdriven, beta, x)
    % Unpack input data
    mix1.overdriven = overdriven;
    mix1.beta = beta; % wave angle  [deg]
    if length(x) > 0
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