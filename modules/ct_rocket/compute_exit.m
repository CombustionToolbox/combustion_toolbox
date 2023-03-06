function mix4 = compute_exit(self, mix2, mix3, mix4, Aratio, varargin)
    % Compute thermochemical composition for a given Aratio
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix2 (struct): Properties of the mixture at the outlet of the chamber
    %     mix3 (struct): Properties of the mixture at the throat
    %     mix4 (struct): Properties of the mixture at the exit (previous calculation)
    %     Aratio (struct): Ratio area_exit / area_throat
    %
    % Returns:
    %     mix3 (struct): Properties of the mixture at the throat

    % Check Aratio
    if isempty(Aratio)
        mix4 = [];
        return
    end

    % Check extra inputs
    if nargin > 5
        mix2_in = varargin{1};
    else
        mix2_in = mix2;
    end

    % Definitions
    self.PD.ProblemType = 'SP';
    % Compute pressure guess [bar] for Infinite-Area-Chamber (IAC)
    logP = guess_pressure_exit_IAC(mix2, mix3, Aratio, self.PD.FLAG_SUBSONIC);
    % Initialization
    STOP = 1; it = 0; logP0 = logP;
    % Loop
    while STOP > self.TN.tol_rocket && it < self.TN.it_rocket
        it = it + 1;
        % Extract pressure [bar]
        pressure = extract_pressure(logP, mix2.p);
        % Solve chemical equilibrium (SP)
        mix4 = equilibrate(self, mix2, pressure, mix4);
        % Compute velocity at the exit point
        mix4.u = compute_velocity(mix2_in, mix4);
        % Compute new estimate
        logP = compute_log_pressure_ratio(mix3, mix4, logP, Aratio);
        % Compute error
        STOP = abs(logP - logP0);
        % Update guess
        logP0 = logP;
        % Debug
        % aux_STOP(it) = STOP;
    end

    % debug_plot_error(it, aux_STOP, 1);
    % Assign values
    mix4.p = extract_pressure(logP, mix2.p); % [bar]
    mix4.v_shock = mix4.u; % [m/s]
    mix4.Aratio = Aratio; % [-]
end

% SUB-PASS FUNCTIONS
function logP = compute_log_pressure_ratio(mix3, mix4, logP, Aratio)
    % Compute pressure [bar]
    Aratio_guess = compute_Aratio(mix3, mix4);
    dlogP = compute_derivative(mix4);
    logP = logP + dlogP * (log(Aratio) - log(Aratio_guess));
end

function velocity = compute_velocity(mix2, mix3)
    % Compute velocity
    velocity = sqrt(2 * (enthalpy_mass(mix2) - enthalpy_mass(mix3)) * 1e3); % [m/s]
end

function value = compute_derivative(mix4)
    % Compute derivative log(Pratio) / log(Aratio)
    value = (mix4.gamma_s * mix4.u^2) / (mix4.u^2 - mix4.sound^2);
end

function value = compute_Aratio(mix3, mix4)
    % Compute area ratio A_exit / A_throat
    value = (mix3.rho * mix3.u) / (mix4.rho * mix4.u);
end

function pressure = extract_pressure(logP, pressure_inf)
    % Extract pressure from log pressure ratio
    pressure = pressure_inf / exp(logP);
end
