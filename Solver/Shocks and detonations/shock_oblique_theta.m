function [str1, str2_1, str2_2] = shock_oblique_theta(varargin)
    % Solve oblique shock

    % Unpack input data
    [self, str1, str2] = unpack(varargin);
    % Abbreviations
    TN = self.TN;
    % Definitions
    a1 = soundspeed(str1);     % [m/s]
    u1 = str1.u;               % [m/s]
    M1 = u1 / a1;              % [-]
    theta = str1.theta * pi/180; % [rad]
    beta_min = asin(1/M1);     % [rad]
    beta_max = pi/2;           % [rad]
    % Solve first branch  (weak shock)
    beta_guess = 0.5 * (beta_min + beta_max); % [rad]
    str2_1 = solve_shock_oblique(str2);
    % Solve second branch (strong shock)
    beta_guess = beta_max * 0.95; % [rad]
    str2_2 = solve_shock_oblique(str2);
    % NESTED FUNCTIONS
    function str2 = solve_shock_oblique(str2) 
        % Obtain wave angle, pre-shock state and post-shock states for an oblique shock
        STOP = 1; it = 0; itMax = TN.it_oblique;
        while STOP > TN.tol_shocks && it < itMax
            it = it + 1;
            % Compute incident normal velocity
            u1n = u1 * sin(beta_guess); % [m/s]
            % Obtain post-shock state
            [~, str2] = shock_incident(self, str1, u1n, str2);
            u2n = str2.v_shock;
            % Compute f0 and df0
            f0 = f0_beta(beta_guess, theta, u2n, u1);
            df0 = df0_beta(beta_guess, u2n, u1);
            % Apply correction
            beta = beta_guess - f0 / df0;
            % Compute error
            STOP = max(abs(beta - beta_guess) / abs(beta), abs(f0));
            % Update guess
            beta_guess = beta;
        end
        
        u2 = u2n * csc(beta - theta);
        M2 = u2 / soundspeed(str2);
    
        if STOP > TN.tol_oblique || beta < beta_min || beta > pi/2 || theta < 0 || M2 >= M1
            error('There is not solution for the given deflection angle %.2g', str1.theta);
        end
    
        % Save results
        str2.beta = beta * 180/pi;   % [deg]
        str2.theta = theta * 180/pi; % [deg]
        str2.beta_min = beta_min * 180/pi; % [deg]
        str2.beta_max = beta_max * 180/pi; % [deg]
        str2.u = u2; % [m/s]
        str2.un = u2n; % [m/s]
        str2.v_shock = u2; % [m/s]
    end
end

% SUB-PASS FUNCTIONS
function [self, str1, str2] = unpack(x)
    % Unpack input data
    self  = x{1};
    str1  = x{2};
    u1    = x{3};
    theta = x{4};       
    str1.u = u1;        % velocity preshock [m/s] - laboratory fixed
    str1.v_shock = u1;  % velocity preshock [m/s] - shock fixed
    str1.theta = theta; % deflection angle  [deg]
    if length(x) > 4
        str2 = x{5};
    else
        str2 = [];
    end
end

function value = f0_beta(beta, theta, u2n, u1)
    % Function to find the roots of beta
    value = theta + atan(u2n / (u1 .* cos(beta))) - beta;
end

function value = df0_beta(beta, u2n, u1)
    % Derivative of the function to find the roots of beta
    value = (u2n * u1 * sin(beta)) / (u2n^2 + u1^2 * cos(beta)^2) - 1;
end