function [mix1, mix2_inj, mix2_c, mix3, mix4] = rocketFAC(obj, mix1, varargin)
    % Compute chemical equilibria at the injector, outlet of the chamber and at the throat
    % using the Finite-Area-Chamber (FAC) model
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     obj (RocketSolver): RocketSolver object
    %     mix1 (Mixture): Properties of the initial mixture
    %
    % Optional Args:
    %     * mix2_inj (Mixture): Properties of the mixture at the injector of the chamber (previous calculation)
    %     * mix2_c (Mixture): Properties of the mixture at the outlet of the chamber (previous calculation)
    %     * mix3 (Mixture): Properties of the mixture at the throat (previous calculation)
    %     * mix4 (Mixture): Properties of the mixture at the given exit points (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix2_inj (Mixture): Properties of the mixture at the injector of the chamber
    %     * mix2_c (Mixture): Properties of the mixture at the outlet of the chamber
    %     * mix3 (Mixture): Properties of the mixture at the throat
    %     * mix4 (Mixture): Array of mixtures objects at the given exit points
    %
    % Example:
    %     [mix2_inj, mix2_c, mix3] = rocketFAC(obj, mix1, mix2_inj, mix2_c, mix3)
    
    % Definitions
    mix2_inj = copy(mix1);
    areaRatio = mix1.areaRatio;
    areaRatioChamber = mix1.areaRatioChamber;
    obj.FLAG_SUBSONIC = true;

    % Unpack additional input parameters
    switch nargin
        case 3
            mix2_inj_guess = varargin{1};
        case 4
            mix2_inj_guess = varargin{1};
            mix2_c_guess = varargin{2};
        case 5
            mix2_inj_guess = varargin{1};
            mix2_c_guess = varargin{2};
            mix3_guess = varargin{3};
        case 6
            mix2_inj_guess = varargin{1};
            mix2_c_guess = varargin{2};
            mix3_guess = varargin{3};
            mix4_guess = varargin{4};
        otherwise
            mix2_inj_guess = [];
            mix2_c_guess = [];
            mix3_guess = [];
            mix4_guess = [];
    end
    
    % Obtain mixture compostion and properties at the injector of the chamber (equivalent to the outlet of the chamber using IAC model)
    mix2_inj = rocketChamberIAC(obj, mix2_inj, mix2_inj_guess);

    % Set areaRatio = areaChamber / areaThroat
    mix2_inj.areaRatio = areaRatioChamber;

    % Get results
    pressure_inj = mix2_inj.p; % [bar]

    % Compute guess pressure_inf
    pressure_inf = obj.rocketGuessInjectorFAC(pressure_inj, areaRatioChamber); % [bar]
    
    % Create a temporal shallow copy of mix1 to avoid overwrite results
    temp_mix1 = mix1.copy;

    % Initialization
    STOP = 1; it = 0;
    pressure_inj = convert_bar_to_Pa(pressure_inj); % [Pa]
    mix2_inf_guess = mix2_inj;
    temp_mix1.areaRatio = areaRatioChamber;
    
    % Loop
    while STOP > obj.tol0 && it < obj.itMax
        % Update iteration
        it = it + 1;

        % Get pressure guess [bar]
        temp_mix1.p = pressure_inf;

        % Obtain mixture compostions
        [~, mix2_inf, mix3, mix2_c] = rocketIAC(obj, temp_mix1, mix2_inf_guess, mix3_guess, mix2_c_guess, false);

        % Get results
        pressure_c = convert_bar_to_Pa(mix2_c.p); % [Pa]
        density_c = mix2_c.rho; % [kg/m3]
        velocity_c = mix2_c.u; % [m/s]

        % Compute pressuse_inj_a
        pressure_inj_a = (pressure_c + density_c * velocity_c^2); % [Pa]

        % Check convergence
        STOP = abs(pressure_inj_a - pressure_inj) / pressure_inj_a;

        % Update guess
        pressure_inf = pressure_inf * pressure_inj / pressure_inj_a; % [bar]
        mix2_inf.p = pressure_inf;
        mix2_inf_guess = mix2_inf;
        mix3_guess = mix3;
        mix2_c_guess = mix2_c;

        % Debug
        % aux_lambda(it) = pressure_inj / pressure_inj_a;
        % aux_STOP(it) = STOP;
    end
    
    % Debug
    % debug_plot_error(it, aux_STOP, aux_lambda);
    
    % Initial velocity of the gas
    mix1.u = 0; mix1.uShock = 0; mix1.mach = 0;
    mix1.I_sp = 0; mix1.I_vac = 0;

    % Velocity at the injector of the chamber
    mix2_inj.u = 0; mix2_inj.uShock = 0; mix2_inj.mach = 0;
    mix2_inj.I_sp = 0; mix2_inj.I_vac = 0;

    % Assign values
    mix2_c.areaRatioChamber = areaRatioChamber; % [-]

    % Compute chemical equilibria at the given exit points
    obj.FLAG_SUBSONIC = false;
    mix4 = rocketExit(obj, mix2_c, mix3, mix4_guess, areaRatio, mix2_inj);

    % Compute rocket parameters
    [mix3, mix2_c, mix4] = obj.rocketParameters(mix2_inj, mix3, mix2_c, mix4);
end