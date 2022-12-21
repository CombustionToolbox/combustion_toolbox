function [mix1, mix2] = det_overdriven(self, mix1, drive_factor, varargin)
    % Compute pre-shock and post-shock states of an overdriven planar detonation
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the mixture in the pre-shock state
    %     drive_factor (float): Overdriven factor [-]
    %
    % Optional Args:
    %     mix2 (struct): Properties of the mixture in the post-shock state (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     - mix1 (struct): Properties of the mixture in the pre-shock state
    %     - mix2 (struct): Properties of the mixture in the post-shock state

    % Unpack input data
    [self, mix1, mix2] = unpack(self, mix1, drive_factor, varargin);
    % Compute CJ speed and initial guess
    if isempty(mix1.cj_speed)
        % Compute CJ speed
        [str1_cj, ~] = det_cj(self, mix1);
        mix1.cj_speed = str1_cj.u;
        % The initial guess is computed as for an incident shock
    end

    % Solve overdriven detonation
    [mix1, mix2] = shock_incident(self, mix1, mix1.cj_speed * mix1.drive_factor, mix2);
    % Assign CJ speed
    mix2.cj_speed = mix1.cj_speed;
end

% SUB-PASS FUNCTIONS
function [self, mix1, mix2] = unpack(self, mix1, drive_factor, x)
    % Unpack input data
    mix1.drive_factor = drive_factor;

    if ~isempty(x)
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
