function [P, T, STOP] = detonationGuess(obj, mix1)
    % Obtain guess of the jump conditions for a Chapman-Jouguet detonation
    % as in NASA's CEA code (see Sec. 8.3 of [1])
    %
    % [1] Gordon, S., & McBride, B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     obj (DetonationSolver): DetonationSolver object
    %     mix1 (Mixture): Properties of the mixture in the pre-shock state
    %
    % Returns:
    %     Tuple containing
    %
    %     * P (float): Pressure ratio [-]
    %     * T (float): Temperature ratio [-]
    %     * STOP (float): Relative error [-]
    %
    % Example:
    %     [P, T, STOP] = detonationGuess(obj, mix1)

    % Definitions
    R0 = combustiontoolbox.common.Constants.R0; % [J/mol-K]
    equilibriumSolver = obj.equilibriumSolver;
    equilibriumSolver.problemType = 'HP';

    % Parameters
    T1 = mix1.T; % [K]
    MW1 = mix1.MW; % [kg/mol]
    h1 = mix1.h / mix1.mi; % [J/kg]
    m1 = mix1.mi; % [kg]
    itMax = obj.itGuess;
    p2p1_0 = 15;

    % Initialization
    mix2 = mix1.copy;

    % Estimate post-shock state considering h2 = h1 + 3/4 * R * T1 / W1 * (p2 / p1)_0
    mix2.h = (h1 + 3/4 * R0 * T1 / MW1 * p2p1_0) * m1; % [J]
    mix2.p = p2p1_0 * mix1.p;

    % Equilibrate HP
    mix2 = equilibriumSolver.equilibrate(mix2);

    % Get properties
    cp2 = mix2.cp ./ mix2.mi; % [J/kg-K]
    gamma2_s = mix2.gamma_s; % [-]
    MW2 = mix2.MW; % [kg]
    
    % Initialization guess (loop)
    P_0 = p2p1_0; % [-]
    T_0 = mix2.T / T1; % [-]
    P = P_0;
    T = T_0;

    % Debug
    % fprintf('T_guess,0 = %.4g\n', T)
    % fprintf('\nP_guess,0 = %.4g\n', P)
    
    % Miscellaneous
    it = 0; STOP = 1;

    % Lopp
    while STOP > obj.tol0 && it < itMax
        % Update iteration
        it = it + 1;
        % Update ratios
        P_old = P;
        T_old = T;
        % Compute ratios
        alpha = 1 / T * MW2 / MW1;
        P = (1 + gamma2_s) / (2 * gamma2_s * alpha) * (1 + sqrt(1 - 4 * gamma2_s * alpha / (1 + gamma2_s)^2));
        r = alpha * P;
        T = T_0 - 3/4 * R0 / (MW1 * cp2) * P_0 + R0 * gamma2_s / (2 * MW1 * cp2) * (r^2 - 1) / r * P;
        % Compute STOP criteria
        STOP = norm([P - P_old, T - T_old]);
    end

    % Debug
    % fprintf('T_guess,%d = %.4g\n', it, T)
    % fprintf('\nP_guess,%d = %.4g\n', it, P)
end