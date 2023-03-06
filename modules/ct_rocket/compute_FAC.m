function [mix2_inj, mix2_c, mix3] = compute_FAC(self, mix1, mix2_inj, mix2_c, mix3)
    % Compute chemical equilibria at the injector, outlet of the chamber and at the throat
    % using the Finite-Area-Chamber (FAC) model
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the initial mixture
    %     mix2_inj (struct): Properties of the mixture at the injector of the chamber (previous calculation)
    %     mix2_c (struct): Properties of the mixture at the outlet of the chamber (previous calculation)
    %     mix3 (struct): Properties of the mixture at the throat (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix2_inj (struct): Properties of the mixture at the injector of the chamber
    %     * mix2_c (struct): Properties of the mixture at the outlet of the chamber
    %     * mix3 (struct): Properties of the mixture at the throat

    % Abbreviations
    TN = self.TN;
    % Definitions
    Aratio_chamber = self.PD.Aratio_c.value;
    % Obtain mixture compostion and properties at the injector of the chamber (equivalent to the outlet of the chamber using IAC model)
    mix2_inj = compute_chamber_IAC(self, mix1, mix2_inj);
    % Set A_chamber/A_throat
    mix2_inj.Aratio = Aratio_chamber;
    % Get results
    pressure_inj = pressure(mix2_inj); % [bar]
    % Compute guess pressure_inf
    pressure_inf = compute_initial_guess_pressure(pressure_inj, Aratio_chamber); % [bar]
    % Initialization
    STOP = 1; it = 0;
    mix2_inf = [];
    pressure_inj = convert_bar_to_Pa(pressure_inj); % [Pa]
    % Loop
    while STOP > TN.tol_rocket && it < TN.it_rocket
        % Update iteration
        it = it + 1;
        % Get guess
        mix1.p = pressure_inf; % [bar]
        % Obtain mixture compostions
        [mix2_inf, mix3, mix2_c] = compute_IAC_model(self, mix1, mix2_inf, mix3, mix2_c, Aratio_chamber);
        % Get results
        pressure_c = convert_bar_to_Pa(pressure(mix2_c)); % [Pa]
        density_c = density(mix2_c); % [kg/m3]
        velocity_c = velocity_relative(mix2_c); % [m/s]
        % Compute pressuse_inj_a
        pressure_inj_a = (pressure_c + density_c * velocity_c^2); % [Pa]
        % Check convergence
        STOP = abs(pressure_inj_a - pressure_inj) / pressure_inj_a;
        % Update guess
        pressure_inf = pressure_inf * pressure_inj / pressure_inj_a; % [bar]
        mix2_inf.p = pressure_inf;
        % Debug
        % aux_lambda(it) = pressure_inj / pressure_inj_a;
        % aux_STOP(it) = STOP;
    end

    % debug_plot_error(it, aux_STOP, aux_lambda);
    % Assign values
    mix2_c.Aratio = Aratio_chamber; % [-]
end

% SUB-PASS FUNCTIONS
function pressure_inf = compute_initial_guess_pressure(pressure_inj, Aratio_chamber)
    % Compute initial guess of the pressure assuming an Infinite-Area-Chamber
    % (IAC) for the Finite-Area-Chamber model (FAC)
    pressure_inf = pressure_inj * (1.0257 - 1.2318 * Aratio_chamber) / (1 - 1.26505 * Aratio_chamber);
end

function [mix2_inf, mix3, mix2_c] = compute_IAC_model(self, mix1, mix2_inf, mix3, mix2_c, Aratio)
    % Solve IAC model
    self.PD.FLAG_IAC = true;
    self.PD.FLAG_SUBSONIC = true;
    mix2_inj = [];
    [~, mix2_inf, mix3, mix2_c] = solve_model_rocket(self, mix1, mix2_inj, mix2_inf, mix3, mix2_c, Aratio);
end
