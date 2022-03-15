function [str1, str2] = overdriven_detonation(varargin)
    % Solve Overdriven detonation 

    % Unpack input data
    [self, str1, str2] = unpack(varargin);
    % Compute CJ speed
    if isempty(str1.cj_speed)
        [str1_cj, ~] = cj_detonation(self, str1);
        str1.cj_speed = str1_cj.u;
    end
    % Solve overdriven detonation
    [str1, str2] = shock_incident(self, str1, str1.cj_speed * str1.overdriven, str2);
    % Assign CJ speed
    str2.cj_speed = str1.cj_speed;
end

% SUB-PASS FUNCTIONS
function [self, str1, str2] = unpack(x)
    % Unpack input data
    self = x{1};
    str1 = x{2};
    overdriven = x{3};
    str1.overdriven = overdriven;
    if length(x) > 3
        str2 = x{4};
    else
        str2 = [];
    end
    if isempty(str2)
        str1.cj_speed = [];
    else
        str1.cj_speed = str2.cj_speed;
    end
end