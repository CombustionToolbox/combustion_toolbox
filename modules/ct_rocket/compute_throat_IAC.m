function mix3 = compute_throat_IAC(self, mix2, mix3)
    % Compute thermochemical composition for the Infinite-Area-Chamber (IAC) model
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix2 (struct): Properties of the mixture at the outlet of the chamber
    %     mix3 (struct): Properties of the mixture at the throat (previous calculation)
    %
    % Returns:
    %     mix3 (struct): Properties of the mixture at the throat

    % Definitions
    self.PD.ProblemType = 'SP';
    % Compute pressure guess [bar] for Infinite-Area-Chamber (IAC)
    pressure = guess_pressure_IAC_model(mix2);
    % Initialization
    STOP = 1; it = 0;
    % Loop
    while STOP > self.TN.tol_rocket && it < self.TN.it_rocket
        it = it + 1;
        mix3 = equilibrate(self, mix2, pressure, mix3);
        mix3.u = compute_velocity(mix2, mix3);
        pressure = compute_pressure(mix3);
        STOP = compute_STOP(mix3);
    end

    % Assign values
    mix3.p = pressure; % [bar]
    mix3.v_shock = mix3.u; % [m/s]
    mix3.Aratio = 1; % [-]
end

function pressure = compute_pressure(mix3)
    % Compute pressure [bar]
    Mach = velocity_relative(mix3) / soundspeed(mix3); % [-]
    pressure = mix3.p * (1 + mix3.gamma_s * Mach^2) / (1 + mix3.gamma_s); % [bar]
end

function velocity = compute_velocity(mix2, mix3)
    % Compute velocity
    velocity = sqrt(2 * (enthalpy_mass(mix2) - enthalpy_mass(mix3)) * 1e3); % [m/s]
end

function STOP = compute_STOP(mix3)
    % Compute STOP criteria
    STOP = abs((mix3.u^2 - mix3.sound^2) / mix3.u^2); % [-]
end
