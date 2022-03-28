function [mix1, mix2, mix5_1, mix5_2] = shock_oblique_reflected_theta(self, mix1, u2, theta, mix2, varargin)
    % Compute pre-shock and post-shock states of an oblique reflected shock wave 
    % given the deflection angle.
    %
    % Two solutions:
    %   * Weak shock
    %   * Strong shock
    %
    % Args:
    %     self (struct):   Data of the mixture, conditions, and databases
    %     mix1 (struct):   Properties of the mixture in the pre-shock state of the incident shock
    %     u2 (float):      Post-shock velocity [m/s] of the incident shock
    %     theta (float):   Deflection angle [deg]
    %     mix2 (struct):   Properties of the mixture in the post-shock state of the incident shock
    %
    % Optional Args:
    %     mix5_1 (struct): Properties of the mixture in the post-shock state of the reflected shock - weak shock (previous calculation)
    %     mix5_2 (struct): Properties of the mixture in the post-shock state of the reflected shock - strong shock (previous calculation)
    %
    % Returns:
    %     mix1 (struct):   Properties of the mixture in the pre-shock state of the incident shock 
    %     mix2 (struct):   Properties of the mixture in the post-shock state of the incident shock 
    %     mix5_1 (struct): Properties of the mixture in the post-shock state of the reflected shock - weak shock
    %     mix5_2 (struct): Properties of the mixture in the post-shock state of the reflected shock - strong shock
    
    % Unpack input data
    [self, mix1, mix2, mix5] = unpack(self, mix1, u2, theta, mix2, varargin);
    % Abbreviations
    TN = self.TN;
    % Definitions
    a2 = soundspeed(mix2);     % [m/s]
    u2 = mix2.u;               % [m/s]
    M2 = u2 / a2;              % [-]
    theta = mix2.theta * pi/180; % [rad]
    beta_min = asin(1/M2);     % [rad]
    beta_max = pi/2;           % [rad]
    % Obtain maximum deflection angle
    [~, mix5_polar] = shock_polar(self, mix2, u2);
    theta_max_reflected = mix5_polar.theta_max * pi/180;
    % Check if theta is in the range for a Regular Reflection (RR)
    if theta > theta_max_reflected
        error(['There is not Regular Reflection (RR) solution for the given conditions:\n' ...
              '   * incident velocity [m/s] %.2f\n' ...
              '   * wave angle [deg]        %.2f\n'], mix1.u, mix2.beta);
    end
    % Solve first branch  (weak shock)
    beta_guess = 0.5 * (beta_min + beta_max); % [rad]
    mix5_1 = solve_shock_oblique(mix5);
    % Solve second branch (strong shock)
    beta_guess = beta_max * 0.95; % [rad]
    mix5_2 = solve_shock_oblique(mix5);
    % NESTED FUNCTIONS
    function mix5 = solve_shock_oblique(mix5) 
        % Obtain wave angle, pre-shock state and post-shock states for an oblique shock
        STOP = 1; it = 0; itMax = TN.it_oblique;
        while STOP > TN.tol_shocks && it < itMax
            it = it + 1;
            
            u2n = u2 * sin(beta_guess); % [m/s]
            [~, mix5] = shock_incident(self, mix2, u2n);
            u5n = mix5.v_shock;
            % Compute f0 and df0
            f0 = f0_beta(beta_guess, theta, u5n, u2);
            df0 = df0_beta(beta_guess, u5n, u2);
            % Compute new value
            beta = beta_guess - f0 / df0;
            % Compute error
            STOP = max(abs(beta - beta_guess) / abs(beta), abs(f0));
            % Update guess
            beta_guess = beta;
        end
        
        u5 = u5n * csc(beta - theta);
        M5 = u5 / soundspeed(mix5);
    
        if STOP > TN.tol_oblique || beta < beta_min || beta > pi/2 || theta < 0 || M5 >= M2
            error('There is not solution for the given deflection angle %.2g', mix5.theta);
        end
    
        % Save results
        mix5.beta = beta * 180/pi;   % [deg]
        mix5.theta = theta * 180/pi; % [deg]
        mix5.beta_min = beta_min * 180/pi; % [deg]
        mix5.beta_max = beta_max * 180/pi; % [deg]
        mix5.u = u5; % [m/s]
        mix5.un = u5n; % [m/s]
        mix5.v_shock = u5; % [m/s]
    end
end

% SUB-PASS FUNCTIONS
function [self, mix1, mix2, mix5] = unpack(self, mix1, u2, theta, mix2, x)
    % Unpack input data
    mix2.u = u2;        % velocity preshock [m/s] - laboratory fixed
    mix2.v_shock = u2;  % velocity preshock [m/s] - shock fixed
    mix2.theta = theta; % deflection angle  [deg]
    if length(x) > 5
        mix5 = x{6};
    else
        mix5 = [];
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