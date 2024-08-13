function mix4 = rocketExit(obj, mix2, mix3, mix4_guess, areaRatio, varargin)
    % Compute thermochemical composition for a given areaRatio
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     obj (RocketSolver): RocketSolver object
    %     mix2 (Mixture): Properties of the mixture at the outlet of the chamber
    %     mix3 (Mixture): Properties of the mixture at the throat
    %     mix4 (Mixture): Properties of the mixture at the exit (previous calculation)
    %     areaRatio (float): Ratio area_exit / area_throat
    %
    % Optional Args:
    %     mix2_inj (Mixture): Properties of the mixture at the injector of the chamber
    %
    % Returns:
    %     mix3 (Mixture): Properties of the mixture at the throat
    %
    % Examples:
    %     * mix4 = rocketExit(RocketSolver(), mix2, mix3, mix4)
    %     * mix4 = rocketExit(RocketSolver(), mix2, mix3, mix4, mix2_inj)

    % Check if areaRatio is provided
    if isempty(areaRatio)
        mix4 = [];
        return
    end
    
    % Initialize mixture
    mix4 = copy(mix3);

    % Check extra inputs
    if nargin > 5
        mix2_inj = varargin{1};
    else
        mix2_inj = mix2;
    end

    % Definitions
    obj.equilibriumSolver.problemType = 'SP';

    % Compute pressure guess [bar] for Infinite-Area-Chamber (IAC)
    logP = obj.rocketGuessExitIAC(mix2, mix3, areaRatio, obj.FLAG_SUBSONIC);

    % Initialization
    STOP = 1; it = 0; logP0 = logP;
    
    % Loop
    while STOP > obj.tol0 && it < obj.itMax
        % Update iteration
        it = it + 1;

        % Extract pressure [bar]
        mix4.p = extractPressure(logP, mix2.p);

        % Solve chemical equilibrium (SP)
        solve(obj.equilibriumSolver, mix4, mix4_guess);

        % Compute velocity at the exit point
        mix4.u = computeVelocity(mix2_inj, mix4);
        mix4.mach = mix4.u / mix4.sound;

        % Compute new estimate
        logP = computeLogPressureRatio(mix3, mix4, logP, areaRatio);

        % Compute error
        STOP = abs(logP - logP0);

        % Update guess
        logP0 = logP;
        mix4_guess = mix4;

        % Debug
        % aux_STOP(it) = STOP;
    end
    
    % Debug
    % debug_plot_error(it, aux_STOP, 1);

    % Assign values
    mix4.p = extractPressure(logP, mix2.p); % [bar]
    mix4.uShock = mix4.u; % [m/s]
    mix4.areaRatio = areaRatio; % [-]
    mix4.problemType = mix3.problemType;
end

% SUB-PASS FUNCTIONS
function logP = computeLogPressureRatio(mix3, mix4, logP, areaRatio)
    % Compute pressure [bar]
    Aratio_guess = computeAreaRatio(mix3, mix4);
    dlogP = computeDerivative(mix4);
    logP = logP + dlogP * (log(areaRatio) - log(Aratio_guess));
end

function velocity = computeVelocity(mix2, mix3)
    % Compute velocity
    velocity = sqrt(2 * (enthalpy_mass(mix2) - enthalpy_mass(mix3)) * 1e3); % [m/s]
end

function value = computeDerivative(mix4)
    % Compute derivative log(Pratio) / log(areaRatio)
    value = (mix4.gamma_s * mix4.u^2) / (mix4.u^2 - mix4.sound^2);
end

function value = computeAreaRatio(mix3, mix4)
    % Compute area ratio A_exit / A_throat
    value = (mix3.rho * mix3.u) / (mix4.rho * mix4.u);
end

function pressure = extractPressure(logP, pressure_inf)
    % Extract pressure from log pressure ratio
    pressure = pressure_inf / exp(logP);
end
