function [str1, str2] = overdriven_detonation(varargin)
% Unpack input data
[self, str1, ~] = unpack(varargin);

[str1, ~] = cj_detonation(self, str1);
[str1, str2] = shock_incident(self, str1, str1.u * str1.overdriven);


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
end