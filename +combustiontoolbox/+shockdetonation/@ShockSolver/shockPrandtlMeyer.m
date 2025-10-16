function [mix1, mix2] = shockPrandtlMeyer(obj, mix1, u1, theta2, varargin)
    % Computes the downstream state associated with an isentropic expansion
    % through a given deflection angle theta2 by marching along the SV (defined
    % entropy and specific volume) expansion path. The objective function
    % thetaGuess - theta2 = 0 is solved using a Newton-Raphson method
    %
    % Args:
    %     obj (ShockSolver): ShockSolver object
    %     mix1 (Mixture): Properties of the mixture in the pre-expansion state
    %     u1 (float): Pre-expansion velocity [m/s]
    %     theta2 (float): Deflection angle [deg]
    %
    % Optional Args:
    %     * mix2 (Mixture): Properties of the mixture in the post-expansion state (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (Mixture): Properties of the initial state
    %     * mix2 (Mixture): Properties of the final state after expansion
    %
    % Notes:
    %     * mix2.polar stores the path history (states along the expansion)
    %
    % Example:
    %     [mix1, mix2] = shockPrandtlMeyer(ShockSolver(), mix1, 15.0)

    % Definitions
    problemType = 'SHOCK_PRANDTL_MEYER';

    % Unpack input data
    [mix1, mix2] = unpack(mix1, u1, theta2, varargin{:});

    % Check supersonic flow
    if mix1.mach < 1
        error('flow is subsonic with Mach %.4f', mix1.mach);
    end

    % Check caloric model
    if ~obj.equilibriumSolver.FLAG_TCHEM_FROZEN || ~obj.equilibriumSolver.FLAG_FROZEN
        obj.equilibriumSolver.problemType = 'TP';
        mach1 = mix1.mach;
        mix1 = solve(obj.equilibriumSolver, mix1);
        mix1.mach = mach1;
        mix1.u = mix1.mach * mix1.sound;
        mix1.uShock = mix1.u;
        mix1.problemType = problemType;
    end

    % Check theta2 value
    if theta2 == 0
        return
    end

    % Definitions
    obj.equilibriumSolver.problemType = 'SV'; % Process at defined entropy and volume
    theta2 = theta2 * pi / 180;               % Deflection angle [rad]
    h0 = mix1.h / mix1.mi + 0.5 * mix1.u^2;   % Stagnation enthalpy [J/kg]
    
    % Get initial guess considering the calorically perfect gas approximation
    [mach2Guess, rho2Guess] = getGuess(mix1, mix2, theta2);

    % Solve problem assuming a calorically perfect gaseous mixture
    if obj.equilibriumSolver.FLAG_TCHEM_FROZEN
        % Definitions
        gamma = mix1.gamma;

        % Isentropic relations for calorically perfect gas (T and p)
        T2 = mix1.T * ( (1 + 0.5 * (gamma - 1) * mix1.mach^2) / (1 + 0.5 * (gamma - 1) * mach2Guess^2) );
        p2 = mix1.p * ( T2 / mix1.T )^(gamma / (gamma - 1));

        % Set state
        mix2.p = p2; mix2.T = T2;
        mix2 = equilibrateT(obj.equilibriumSolver, mix1, mix2, T2);

        % Final state
        mix2.mach = mach2Guess;          % [-]
        mix2.u = mix2.mach * mix2.sound; % [m/s]
        mix2.uShock = mix2.u;            % [m/s]
        mix2.theta = theta2 * 180 / pi;  % [deg]
        mix2.problemType = problemType;
        return
    end

    % Solve problem assuming a thermally perfect gaseous mixture
    if obj.equilibriumSolver.FLAG_FROZEN
        mix2 = setProperties(mix2, 'entropySpecific', mix1.sSpecific, 'volume', 1/rho2Guess, 'mach', mach2Guess);

        % Final state
        mix2.uShock = mix2.u;           % [m/s]
        mix2.theta = theta2 * 180 / pi; % [deg]
        mix2.problemType = problemType;
        return
    end

    % Solve problem assuming a calorically imperfect gaseous mixture at chemical equilibrium

    % Initialization
    it = 0; itMax = obj.itMax; STOP = 1;
    
    % Batch defined entropy-volume (SV) integration from rho1 to rho2Guess to build initial path
    [mixArray2, rhoArray2, machArray2, integrandArray] = tracePrandtlMeyerExpansion(obj, mix1, rho2Guess, h0, theta2);

    % Get guess from final state
    mix2 = mixArray2(end);
    rho2Guess = rhoArray2(end);

    while STOP > obj.tol0 && it < itMax
        % Update iteration
        it = it + 1;

        % Compute new entropy-volume state
        mix2 = stateSV(obj.equilibriumSolver, mix2, mix1.sSpecific, rho2Guess, mix2);

        % Update flow velocity and Mach number
        mix2.u = sqrt( 2 * (h0 - mix2.h / mix2.mi) );
        mix2.uShock = mix2.u;
        machArray2(end + 1) = mix2.u ./ mix2.sound;

        % Compute integrand Prandt-Meyer function
        integrandArray(end + 1) = getIntegrandPrandtMeyer(rhoArray2(end), machArray2(end));
        
        % Sort data (descending density) to preserve monotonic integration
        [rhoArray2, index] = sort(rhoArray2, 'descend');
        machArray2 = machArray2(index);
        integrandArray = integrandArray(index);

        % Select domain up to the current density guess
        FLAG_INTEGRATION = rhoArray2 >= rho2Guess;

        % Current deflection angle from trapezoidal integration
        theta = trapz( rhoArray2(FLAG_INTEGRATION), integrandArray(FLAG_INTEGRATION) );

        % Compute residual and derivative
        f0 = theta - theta2;
        df = integrandArray( find(FLAG_INTEGRATION, 1, 'last') );

        % Apply correction
        rho2 = rho2Guess - f0 / df;
        
        % Compute error
        STOP = max(abs( (rho2 - rho2Guess) / rho2), abs(f0));

        % Update guess
        rho2Guess = rho2;

        % Append state
        mix2.theta = theta * 180 / pi; % [deg]
        rhoArray2(end + 1) = rho2;  % [kg/m3]
        mixArray2(end + 1) = mix2;
    end
    
    % Final state
    mix2.polar = mixArray2;
    mix2.problemType = problemType;
end

% SUB-PASS FUNCTIONS
function [mix1, mix2] = unpack(mix1, u1, theta, varargin)
    % Unpack input data
    mix1.u = u1;                     % Pre-shock velocity [m/s] - laboratory fixed
    mix1.uShock = u1;                % Pre-shock velocity [m/s] - shock fixed
    mix1.mach = mix1.u / mix1.sound; % Pre-shock Mach number [-]
    mix1.theta = theta;              % Wave angle  [deg]

    if nargin > 3
        mix2 = varargin{1}.copy();
    else
        mix2 = mix1.copy();
    end
end

function [mach2Guess, rho2Guess] = getGuess(mix1, mix2, theta2)
    % Get initial guess for the downstream state based on the previous 
    % calculation or the calorically perfect gas approximation
    %
    % Args:
    %     mix1 (Mixture): Properties of the mixture in the pre-expansion state
    %     mix2 (Mixture): Properties of the mixture in the post-expansion state (previous calculation)
    %     theta2 (float): Deflection angle [rad]
    %
    % Returns:
    %     Tuple containing
    %
    %     * mach2Guess (float): Initial guess for the Mach number of the post-expansion state [-]
    %     * rho2Guess (float): Initial guess for the density of the post-expansion state [kg/m3]

    % Get initial guess from previous calculation
    % if mix2.mach ~= mach1
    %     mach2Guess = mix2.mach;
    %     rho2Guess = mix2.rho;
    %     return
    % end

    % Definitions
    rho1 = mix1.rho;    % Pre-expansion Density [kg/m3]
    mach1 = mix1.mach;  % Pre-expansion Mach number
    gamma = mix1.gamma; % Pre-expansion adiabatic index [-]

    % Calorically perfect gas approximation
    nu1 = getPrandtlMeyerPerfect(mach1, gamma); % [rad]
    nu2 = nu1 + theta2; % [rad]
    mach2Guess = getInversePrandtlMeyerPerfect(nu2, gamma);
    rho2Guess = rho1 * ( (1 + 0.5 * (gamma - 1) * mach1^2) / ...
                         (1 + 0.5 * (gamma - 1) * mach2Guess^2) )^(1 / (gamma - 1));
end

function [mixArray2, rhoArray2, machArray2, integrandArray] = tracePrandtlMeyerExpansion(obj, mix1, rho2, h0, theta2)
    % Trace Prandtl–Meyer expansion path along SV (defined entropy and volume)
    %
    % Generates a discretized expansion path from the initial density rho1 down
    % to a target density guess rho2 by marching along the SV equilibrium curve.
    % Each intermediate state is computed by enforcing (entropy, specific volume)
    % and recovering flow velocity and Mach number from the stagnation enthalpy 
    % constraint
    %
    % Args:
    %     obj (ShockSolver): ShockSolver object
    %     mix1 (Mixture): Properties of the mixture in the pre-shock
    %     rho2 (float): Final density [kg/m3]
    %     h0 (float): Stagnation enthalpy [J/kg]
    %     theta2 (float): Deflection angle [rad]
    %
    % Returns:
    %     Tuple containing
    %     
    %     * mixArray2 (Mixture): Array of Mixture objects with the states along the expansion
    %     * rhoArray2 (float): Density vector [kg/m3]
    %     * machArray2 (float): Mach number vector [-]
    %     * integrandArray (float): Integrand of the Prandt-Meyer function [rad-m3/kg]

    % Definitions
    numCases = obj.numPointsPrandtlMeyer; % Number of points for integration
    decay = 2;                            % Controls clustering intensity (1 = uniform, >1 = more clustering near rho1)
    xi = linspace(0, 1, numCases);        % Normalized coordinate
    rhoArray2 = mix1.rho - (mix1.rho - rho2) * (1 - exp(-decay * xi)) / (1 - exp(-decay)); % Density vector [kg/m3]

    % Preallocation
    integrandArray = zeros( size(rhoArray2) );
    machArray2 = mix1.mach * ones( size(rhoArray2) );
    mixArray2 = mix1.empty(0, numCases);
    mixArray2(1) = mix1.copy();

    % Compute states
    for i = 2:numCases
        mixArray2(i) = stateSV(obj.equilibriumSolver, mixArray2(i - 1), mix1.sSpecific, rhoArray2(i), mixArray2(i - 1));
        u = sqrt( 2 * (h0 - mixArray2(i).h / mixArray2(i).mi) );
        machArray2(i) = u ./ mixArray2(i).sound;

        % Assign properties
        mixArray2(i).u = u;
        mixArray2(i).uShock = u;
        mixArray2(i).mach = machArray2(i);
    end

    % Compute integrand Prandt-Meyer function
    integrandArray(2:end) = getIntegrandPrandtMeyer(rhoArray2(2:end), machArray2(2:end));

    % Trapezoidal integration
    theta = cumtrapz(rhoArray2, integrandArray);

    % Get residual and derivative
    f0 = theta(end) - theta2;
    df = integrandArray(end);

    % Apply correction
    rhoArray2(end + 1) = rho2 - f0 / df;

    % Assign deflection angles
    for i = 1:numCases
        mixArray2(i).theta = theta(i) * 180 / pi;
    end

end

function mix = stateSV(equilibriumSolver, mix, entropySpecific, density, mixGuess)
    % Compute thermochemical equilibrium state at given entropy and density
    %
    % Args:
    %     equilibriumSolver (EquilibriumSolver): EquilibriumSolver object
    %     mix (Mixture): Mixture object with initial guess
    %     entropySpecific (float): entropySpecific [J/kg-K]
    %     density (float): Density [kg/m3]
    %     mixGuess (Mixture): Previous mixture state (for better convergence)
    %
    % Returns:
    %     mix (Mixture): Mixture object with updated properties

    % Definitions
    problemType = mix.problemType;

    % Set state
    mix = setProperties(mix, 'entropySpecific', entropySpecific, 'volume', 1 / density);
    
    % Solve equilibrium state
    solve(equilibriumSolver, mix, mixGuess);

    % Restore problem type
    mix.problemType = problemType;
end

function integrand = getIntegrandPrandtMeyer(rho, mach)
    % Compute integrand of the Prandt-Meyer function
    %
    % Args:
    %     rho (float): Density [kg/m3]
    %     mach (float): Mach number [-]
    %
    % Returns:
    %     integrand (float): Integrand of the Prandt-Meyer function [rad-m3/kg]
    
    integrand = - sqrt(mach.^2 - 1) ./ (mach.^2 .* rho);
end

function nu = getPrandtlMeyerPerfect(mach, gamma)
    % Prandtl-Meyer function for calorically perfect gas
    %
    % Args:
    %     mach (float): Mach number [-]
    %     gamma (float): Adiabatic index [-]
    %
    % Returns:
    %     nu (float): Prandtl-Meyer function value [rad]

    nu = sqrt( (gamma + 1) / (gamma - 1) ) * atan( sqrt( (gamma - 1) / (gamma + 1) * (mach.^2 - 1) ) ) - atan( sqrt(mach.^2 - 1) );
end

function mach = getInversePrandtlMeyerPerfect(nu, gamma)
    % Solve inverse of Prandtl-Meyer function for calorically perfect gas
    % using Newton's method
    %
    % Args:
    %     nu (float): Prandtl-Meyer function value [rad]
    %     gamma (float): Adiabatic index [-]
    %
    % Returns:
    %     mach (float): Mach number [-]

    % Definitions
    tol0 = 1e-5; itMax = 100;

    % Initialization
    machGuess = 2; STOP = 1; it = 0;

    % Newton-Raphson iteration
    while STOP > tol0 && it < itMax
        % Update iteration
        it = it + 1;
        
        % Get residual and derivative
        f = getPrandtlMeyerPerfect(machGuess, gamma) - nu;
        df = sqrt(machGuess^2 - 1) / (machGuess * (1 + 0.5 * (gamma - 1) * machGuess^2));

        % Apply correction
        mach = machGuess - f/df;

        % Compute error
        STOP = abs( (mach - machGuess) / mach );

        % Update guess
        machGuess = mach;
    end

    if STOP > tol0 
        warning('getInversePrandtlMeyerPerfect: No convergence after %d iterations', itMax);
    end

end