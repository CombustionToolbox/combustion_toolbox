function [mix1, mix2] = overdriven_detonation(varargin)
    % Solve Overdriven detonation 

    % Unpack input data
    [self, mix1, mix2] = unpack(varargin);
    % Compute CJ speed
    if isempty(mix1.cj_speed)
        [str1_cj, ~] = cj_detonation(self, mix1);
        mix1.cj_speed = str1_cj.u;
    end
    % Solve overdriven detonation
    [mix1, mix2] = shock_incident(self, mix1, mix1.cj_speed * mix1.overdriven, mix2);
    % Assign CJ speed
    mix2.cj_speed = mix1.cj_speed;
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