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
function [self, mix1, mix2] = unpack(x)
    % Unpack input data
    self = x{1};
    mix1 = x{2};
    overdriven = x{3};
    mix1.overdriven = overdriven;
    if length(x) > 3
        mix2 = x{4};
    else
        mix2 = [];
    end
    if isempty(mix2)
        mix1.cj_speed = [];
    else
        mix1.cj_speed = mix2.cj_speed;
    end
end