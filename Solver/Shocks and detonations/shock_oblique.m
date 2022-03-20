function [str1, str2] = shock_oblique(varargin)
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

    w2_guess = 0.2 * u1;       % [m/s]
    beta_guess = 0.5 * (beta_min + beta_max); % [rad]
    beta_guess = theta + atan(w2_guess / (u1 * cos(beta_guess)));
    beta = beta_guess;

    % Solve first case for initialization
    STOP = 1; it = 0; itMax = TN.it_oblique;
    while STOP > TN.tol_shocks && it < itMax
        it = it + 1;

        w1 = u1 * sin(beta);       % [m/s]
        v = u1 .* cos(beta);       % [m/s]
        if ~contains(self.PD.ProblemType, '_R')
            [~, str2] = shock_incident(self, str1, w1, str2);
        else
            [~, ~, str2] = shock_reflected(self, str1, str2.w2, str2);
        end
    
        w2 = str2.v_shock;
        beta_guess = theta + atan(w2 / v);
        beta = theta + atan(w2 / (u1 * cos(beta_guess)));
        STOP = abs(beta - beta_guess);
    end
    
    u2 = w2 * csc(beta - theta);
    M2 = u2 / soundspeed(str2);

    if STOP > TN.tol_oblique || beta < beta_min || beta > pi/2 || theta < 0 || M2 >= M1
        error('There is not solution for the given deflection angle');
    end

%     if theta < 0 || theta > theta_max
%         error('There is not solution for the given deflection angle');
%     end

% NEWWW
    w2_range = linspace(0, u1, 200);
    beta_range = linspace(beta_min, beta_max, 200);
    v_range = u1 * cos(beta_range);
    
    ind_max = find(beta_min - atan(w2_range ./ v_range) > 0, 1, 'last');
    theta_range = beta_range - atan(w2_range(ind_max) ./ v_range);
    [theta_max, ind_max]= max(theta_range);
    beta_max = beta_range(ind_max);
% NEWWW



%     w2_range = linspace(0, u1, 200);
%     beta_range = linspace(beta_min, beta_max, 200);
%     v_range = u1 * cos(beta_range);
    
%     figure; hold on;
%     for i=1:1:200
%         plot(theta_range(i, :) * 180/pi, beta_range* 180/pi, 'k');
%     end
%     i = length(w2_range); FLAG = true;
%     while FLAG
%         i = i - 1;
%         theta_range = beta_range - atan(w2_range(i) ./ v_range);
%         FLAG = any(theta_range < 0);
%         plot(theta_range * 180/pi, beta_range* 180/pi, 'k');
%         pause(0.1)
%     end
%     theta_max = max(max(theta_range));
%     
%     plot(theta_range * 180/pi, beta_range* 180/pi, 'r');
    
    % Sonic point and Max deflection angle
%     w2_range = linspace(0, w1, 200);
%     theta_range = beta - atan(w2_range / v);
%     u2_range = w2_range .* csc(beta - theta_range);
%     M2_range = u2_range ./ a2;
%     ind_sonic = find(M2_range > 1, 1, 'last');
%     theta_sonic = theta(ind_sonic) * 180/pi;
%     [theta_max, ind_max] = max(theta_range * 180/pi);
    % Save results
    str1.beta = beta * 180/pi;   % [deg]
    str1.theta = theta * 180/pi; % [deg]
    str1.theta_max = theta_max * 180/pi; % [deg]
    str1.beta_min = beta_min * 180/pi; % [deg]
    str1.beta_max = beta_max * 180/pi; % [deg]
    str1.theta_sonic = 0;        % [deg]
    str1.beta_sonic = 0;         % [deg]
    str2.u2 = u2; % [m/s]
    str2.w2 = w2; % [m/s]
end

% SUB-PASS FUNCTIONS
function [self, str1, str2] = unpack(x)
    % Unpack input data
    self  = x{1};
    str1  = x{2};
    u1    = x{3};
    theta = x{4};
    str1.u = u1;       % velocity preshock [m/s] - laboratory fixed
    str1.v_shock = u1; % velocity preshock [m/s] - shock fixed
    str1.theta = theta; % detached angle [deg]
    if length(x) > 4
        str2 = x{5};
    else
        str2 = [];
    end
end

function f = fun2solve(x, u1, theta)
    f = [x(1) - theta - atan(x(2) / (u1 * cos(x(1))))
         cos(x(1)) - x(2) / u1 * cot(x(1) - theta)];
end