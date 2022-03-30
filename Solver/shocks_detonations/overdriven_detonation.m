function [mix1, mix2] = overdriven_detonation(self, mix1, overdriven, varargin)
    % Compute pre-shock and post-shock states of an overdriven planar detonation
    %
    % Args:
    %     self (struct):      Data of the mixture, conditions, and databases
    %     mix1 (struct):      Properties of the mixture in the pre-shock state
    %     overdriven (float): Overdriven factor [-]
    %
    % Optional Args:
    %     mix2 (struct):      Properties of the mixture in the post-shock state (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     - mix1 (struct):      Properties of the mixture in the pre-shock state
    %     - mix2 (struct):      Properties of the mixture in the post-shock state

    % Unpack input data
    [self, mix1, mix2] = unpack(self, mix1, overdriven, varargin);
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
function [self, mix1, mix2] = unpack(self, mix1, overdriven, x)
    % Unpack input data
    mix1.overdriven = overdriven;
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