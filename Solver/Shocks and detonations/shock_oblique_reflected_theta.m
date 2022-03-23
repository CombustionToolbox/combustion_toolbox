function [str1, str2, str5_1, str5_2] = shock_oblique_reflected_theta(varargin)
    % Solve oblique shock

    % Unpack input data
    [self, str1, str2, str5] = unpack(varargin);
    % Abbreviations
    TN = self.TN;
    % Definitions
    a2 = soundspeed(str2);     % [m/s]
    u2 = str2.u;               % [m/s]
    M2 = u2 / a2;              % [-]
    theta = str2.theta * pi/180; % [rad]
    beta_min = asin(1/M2);     % [rad]
    beta_max = pi/2;           % [rad]
    % Obtain maximum deflection angle
    [~, str5_polar] = shock_polar(self, str2, u2);
    theta_max_reflected = str5_polar.theta_max * pi/180;
    % Check if theta is in the range for a Regular Reflection (RR)
    if theta > theta_max_reflected
        error(['There is not Regular Reflection (RR) solution for the given conditions:\n' ...
              '   * incident velocity [m/s] %.2f\n' ...
              '   * wave angle [deg]        %.2f\n'], str1.u, str2.beta * 180/pi);
    end
    % Solve first branch  (weak shock)
    beta_guess = 0.5 * (beta_min + beta_max); % [rad]
    str5_1 = solve_shock_oblique(str5);
    % Solve second branch (strong shock)
    beta_guess = beta_max * 0.95; % [rad]
    str5_2 = solve_shock_oblique(str5);
    % NESTED FUNCTIONS
    function str5 = solve_shock_oblique(str5) 
        % Obtain wave angle, pre-shock state and post-shock states for an oblique shock
        STOP = 1; it = 0; itMax = TN.it_oblique;
        while STOP > TN.tol_shocks && it < itMax
            it = it + 1;
            
            u2n = u2 * sin(beta_guess); % [m/s]
            [~, str5] = shock_incident(self, str2, u2n);
            u5n = str5.v_shock;
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
        M5 = u5 / soundspeed(str5);
    
        if STOP > TN.tol_oblique || beta < beta_min || beta > pi/2 || theta < 0 || M5 >= M2
            error('There is not solution for the given deflection angle %.2g', str5.theta);
        end
    
        % Save results
        str5.beta = beta * 180/pi;   % [deg]
        str5.theta = theta * 180/pi; % [deg]
        str5.beta_min = beta_min * 180/pi; % [deg]
        str5.beta_max = beta_max * 180/pi; % [deg]
        str5.u = u5; % [m/s]
        str5.un = u5n; % [m/s]
        str5.v_shock = u5; % [m/s]
    end
end

% SUB-PASS FUNCTIONS
function [self, str1, str2, str5] = unpack(x)
    % Unpack input data
    self  = x{1};
    str1  = x{2};
    u2    = x{3};
    theta = x{4};  
    str2  = x{5};
    str2.u = u2;        % velocity preshock [m/s] - laboratory fixed
    str2.v_shock = u2;  % velocity preshock [m/s] - shock fixed
    str2.theta = theta; % deflection angle  [deg]
    if length(x) > 5
        str5 = x{6};
    else
        str5 = [];
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