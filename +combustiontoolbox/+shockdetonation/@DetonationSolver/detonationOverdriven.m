function [mix1, mix2] = detonationOverdriven(obj, mix1, driveFactor, varargin)
    % Compute pre-shock and post-shock states of an overdriven planar detonation
    %
    % Args:
    %     obj (DetonationSolver): DetonationSolver object
    %     mix1 (Mixture): Properties of the mixture in the pre-shock state
    %     driveFactor (float): Overdriven factor [-]
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
    %     * [mix1, mix2] = detonationOverdriven(DetonationSolver(), mix1, 1.5)
    %     * [mix1, mix2] = detonationOverdriven(DetonationSolver(), mix1, 1.5, mix2)
    
    % Import packages
    import combustiontoolbox.shockdetonation.ShockSolver

    % Unpack input data
    [mix1, mix2] = unpack(mix1, driveFactor, varargin);
    
    % Compute CJ speed and initial guess
    if isempty(mix1.cjSpeed)
        % Compute CJ speed
        [mix1_cj, ~] = detonationCJ(obj, mix1);
        mix1.cjSpeed = mix1_cj.u;
        % The initial guess is computed as for an incident shock
    end

    % Solve overdriven detonation
    [mix1, mix2] = shockIncident(ShockSolver(), mix1, mix1.cjSpeed * mix1.driveFactor, mix2);

    % Assign CJ speed
    mix2.cjSpeed = mix1.cjSpeed;
end

% SUB-PASS FUNCTIONS
function [mix1, mix2] = unpack(mix1, driveFactor, x)
    % Unpack input data
    mix1.driveFactor = driveFactor;

    if ~isempty(x)
        mix2 = x{1};
    else
        mix2 = [];
    end

    if isempty(mix2)
        mix1.cjSpeed = [];
    else
        mix1.cjSpeed = mix2.cjSpeed;
    end
    
end
