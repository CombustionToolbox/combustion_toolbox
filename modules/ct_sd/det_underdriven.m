function [mix1, mix2] = det_underdriven(self, mix1, drive_factor, varargin)
    % Compute pre-shock and post-shock states of an overdriven planar detonation
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the mixture in the pre-shock state
    %     drive_factor (float): Underdriven factor [-]
    %
    % Optional Args:
    %     mix2 (struct): Properties of the mixture in the post-shock state (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (struct): Properties of the mixture in the pre-shock state
    %     * mix2 (struct): Properties of the mixture in the post-shock state
    
    % Definitions
    zeta = 0.1;
    % Unpack input data
    [self, mix1, mix2] = unpack(self, mix1, drive_factor, varargin);
    if isempty(mix1.cj_speed)
        % Compute CJ state
        [str1_cj, str2_cj] = det_cj(self, mix1);
        mix1.cj_speed = str1_cj.u;
        % Adjust guess for underdriven detonation
        u1 = mix1.cj_speed * drive_factor; % [m/s]
        p1 = convert_bar_to_Pa(mix1.p); % [Pa]
        mix2.rho = 1 / (zeta / str2_cj.rho + (1 - zeta) / mix1.rho);
        mix2.p = convert_Pa_to_bar(p1 - mix1.rho * u1^2 * (mix1.rho / mix2.rho - 1));
        mix2.Xi = str2_cj.Xi;
        mix2.N = str2_cj.N;
        mix2.T = mix2.p / (mix2.rho * self.C.R0 / (str2_cj.W)) * 1e2; % assuming ideal EoS
    end
    
    % Solve detonation
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
