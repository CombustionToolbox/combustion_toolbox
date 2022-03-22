function [str1, str2] = shock_oblique_beta(varargin)
    % Solve oblique shock for a given shock wave angle
    
    % Unpack input data
    [self, str1, str2] = unpack(varargin);
    % Definitions
    a1 = soundspeed(str1);     % [m/s]
    u1 = str1.u;               % [m/s]
    M1 = u1 / a1;              % [-]
    beta = str1.beta * pi/180; % [rad]
    beta_min = asin(1/M1);     % [rad]
    beta_max = pi/2;           % [rad]
    w1 = u1 * sin(beta);       % [m/s]
    % Check range beta
    if beta < beta_min || beta > beta_max
        error('\nERROR! The given wave angle beta = %.2g is not in the range of possible solutions [%.2g, %2.g]', str1.beta, beta_min * 180/pi, beta_max * 180/pi);
    end
    
    % Obtain deflection angle, pre-shock state and post-shock states for an oblique shock
    if ~contains(self.PD.ProblemType, '_R')
        [~, str2] = shock_incident(self, str1, w1, str2);
    else
        [~, ~, str2] = shock_reflected(self, str1, str2.w2, str2);
    end
    
    w2 = str2.v_shock;
    theta = beta - atan(w2 / (u1 .* cos(beta)));
    u2 = w2 * csc(beta - theta);
    
    % Save results
    str2.beta = beta * 180/pi;   % [deg]
    str2.theta = theta * 180/pi; % [deg]
    str2.beta_min = beta_min * 180/pi; % [deg]
    str2.beta_max = beta_max * 180/pi; % [deg]
    str2.u = u2; % [m/s]
    str2.w2 = w2; % [m/s]
    str2.v_shock = u2; % [m/s]
end

% SUB-PASS FUNCTIONS
function [self, str1, str2] = unpack(x)
    % Unpack input data
    self  = x{1};
    str1  = x{2};
    u1    = x{3};
    beta = x{4};       
    str1.u = u1;        % velocity preshock [m/s] - laboratory fixed
    str1.v_shock = u1;  % velocity preshock [m/s] - shock fixed
    str1.beta = beta;   % wave angle  [deg]
    if length(x) > 4
        str2 = x{5};
    else
        str2 = [];
    end
end