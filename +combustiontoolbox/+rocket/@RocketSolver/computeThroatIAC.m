function mix3 = computeThroatIAC(obj, mix2, mix3)
    % Compute thermochemical composition for the Infinite-Area-Chamber (IAC) model
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     obj (RocketSolver): RocketSolver object
    %     mix2 (Mixture): Properties of the mixture at the outlet of the chamber
    %     mix3 (Mixture): Properties of the mixture at the throat (previous calculation)
    %
    % Returns:
    %     mix3 (Mixture): Properties of the mixture at the throat
    %
    % Example:
    %     mix3 = computeThroatIAC(RocketSolver(), mix2, mix3)

    % Definitions
    obj.equilibriumSolver.problemType = 'SP';
    
    % Check input
    if isempty(mix3)
        mix3 = copy(mix2);
    end

    % Compute pressure guess [bar] for Infinite-Area-Chamber (IAC)
    pressure = guess_pressure_IAC_model(mix2);
    mix3.p = pressure;

    % Initialization
    STOP = 1; it = 0;
    
    % Loop
    while STOP > obj.tol0 && it < obj.itMax
        it = it + 1;
        solve(obj.equilibriumSolver, mix3, mix3);
        mix3.u = compute_velocity(mix2, mix3);
        mix3.mach = mix3.u / mix3.sound;
        pressure = compute_pressure(mix3);
        mix3.p = pressure;
        STOP = compute_STOP(mix3);
    end

    % Assign values
    mix3.p = pressure; % [bar]
    mix3.uShock = mix3.u; % [m/s]
    mix3.mach = mix3.u / mix3.sound; % [-]
    mix3.areaRatio = 1; % [-]
end

function pressure = compute_pressure(mix3)
    % Compute pressure [bar]
    pressure = mix3.p * (1 + mix3.gamma_s * mix3.mach^2) / (1 + mix3.gamma_s); % [bar]
end

function velocity = compute_velocity(mix2, mix3)
    % Compute velocity
    velocity = sqrt(2 * (enthalpy_mass(mix2) - enthalpy_mass(mix3)) * 1e3); % [m/s]
end

function STOP = compute_STOP(mix3)
    % Compute STOP criteria
    STOP = abs((mix3.u^2 - mix3.sound^2) / mix3.u^2); % [-]
end

function pressure = guess_pressure_IAC_model(mix)
    % Compute pressure guess [bar] at the throat considering an Infinite-Area-Chamber (IAC)
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     mix (Mixture): Properties of the mixture
    %
    % Returns:
    %     pressure (float): Pressure at the throat [bar]
    %
    % Example:
    %     pressure = guess_pressure_IAC_model(mix)

    pressure = mix.p / ((mix.gamma_s + 1) / 2)^(mix.gamma_s / (mix.gamma_s - 1));
end
