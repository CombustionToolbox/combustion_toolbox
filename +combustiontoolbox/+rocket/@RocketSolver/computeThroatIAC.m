function mix3 = computeThroatIAC(obj, mix2, mix3_guess)
    % Compute thermochemical composition for the Infinite-Area-Chamber (IAC) model
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     obj (RocketSolver): RocketSolver object
    %     mix2 (Mixture): Properties of the mixture at the outlet of the chamber
    %     mix3_guess (Mixture): Properties of the mixture at the throat (previous calculation)
    %
    % Returns:
    %     mix3 (Mixture): Properties of the mixture at the throat
    %
    % Examples:
    %     * mix3 = computeThroatIAC(RocketSolver(), mix2, mix3_guess)
    %     * mix3 = computeThroatIAC(RocketSolver(), mix2, [])

    % Definitions
    obj.equilibriumSolver.problemType = 'SP';

    % Initialize mixture
    mix3 = copy(mix2);

    % Compute pressure guess [bar] for Infinite-Area-Chamber (IAC)
    mix3.p = guess_pressure_IAC_model(mix2);

    % Initialization
    STOP = 1; it = 0;
    
    % Loop
    while STOP > obj.tol0 && it < obj.itMax
        % Update iteration
        it = it + 1;

        % Solve chemical equilibrium (SP)
        solve(obj.equilibriumSolver, mix3, mix3_guess);

        % Compute velocity at the throat
        mix3.u = compute_velocity(mix2, mix3);
        mix3.mach = mix3.u / mix3.sound;

        % Compute new estimate
        mix3.p = compute_pressure(mix3);

        % Compute STOP criteria
        STOP = compute_STOP(mix3);

        % Update guess
        mix3_guess = mix3;
    end

    % Assign values
    mix3.uShock = mix3.u; % [m/s]
    mix3.areaRatio = 1; % [-]
    mix3.problemType = mix2.problemType;
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
