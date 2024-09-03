function [mix1, mix2] = detonationUnderdriven(obj, mix1, driveFactor, varargin)
    % Compute pre-shock and post-shock states of an overdriven planar detonation
    %
    % Args:
    %     obj (DetonationSolver): DetonationSolver object
    %     mix1 (Mixture): Properties of the mixture in the pre-shock state
    %     driveFactor (float): Underdriven factor [-]
    %
    % Optional Args:
    %     mix2 (Mixture): Properties of the mixture in the post-shock state (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (Mixture): Properties of the mixture in the pre-shock state
    %     * mix2 (Mixture): Properties of the mixture in the post-shock state
    %
    % Examples:
    %     * [mix1, mix2] = detonationUnderdriven(DetonationSolver(), mix1, 1.5)
    %     * [mix1, mix2] = detonationUnderdriven(DetonationSolver(), mix1, 1.5, mix2)
    
    % Definitions
    R0 = combustiontoolbox.common.Constants.R0;
    zeta = 0.1;
    
    % Unpack input data
    [mix1, mix2] = unpack(mix1, driveFactor, varargin);

    if isempty(mix1.cjSpeed)
        % Compute CJ state
        [mix1_cj, mix2_cj] = detonationCJ(obj, mix1);
        mix1.cjSpeed = mix1_cj.u;
        % Adjust guess for underdriven detonation
        u1 = mix1.cjSpeed * driveFactor; % [m/s]
        p1 = convert_bar_to_Pa(mix1.p); % [Pa]
        mix2.rho = 1 / (zeta / mix2_cj.rho + (1 - zeta) / mix1.rho);
        mix2.p = convert_Pa_to_bar(p1 - mix1.rho * u1^2 * (mix1.rho / mix2.rho - 1));
        mix2.Xi = mix2_cj.Xi;
        mix2.N = mix2_cj.N;
        mix2.T = mix2.p / (mix2.rho * R0 / (mix2_cj.W)) * 1e5; % assuming ideal EoS
    end
    
    % Solve detonation
    [mix1, mix2] = shockIncident(obj.shockSolver, mix1, mix1.cjSpeed * mix1.driveFactor, mix2);

    % Assign CJ speed and driveFactor
    mix2.cjSpeed = mix1.cjSpeed;
    mix2.driveFactor = mix1.driveFactor;
end

% SUB-PASS FUNCTIONS
function [mix1, mix2] = unpack(mix1, driveFactor, x)
    % Unpack input data
    mix1.driveFactor = driveFactor;

    if ~isempty(x)
        mix2 = x{1};
    else
        mix2 = mix1.copy;
    end

    if isempty(mix2)
        mix1.cjSpeed = [];
    else
        mix1.cjSpeed = mix2.cjSpeed;
    end

end
