function [P, T, STOP] = det_compute_guess_CEA(self, mix1)
    % Obtain guess of the jump conditions for a Chapman-Jouguet detonation
    % as in NASA's CEA code (see Sec. 8.3 of [1])
    %
    % [1] Gordon, S., & McBride, B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the mixture in the pre-shock state
    %
    % Returns:
    %     Tuple containing
    %
    %     * P (float): Pressure ratio [-]
    %     * T (float): Temperature ratio [-]
    %     * STOP (float): Relative error [-]

    % Abbreviations
    R0 = self.C.R0; % [J/mol-K]
    % Parameters
    T1 = temperature(mix1); % [K]
    W1 = meanMolecularWeight(mix1) * 1e-3; % [kg/mol]
    h1 = enthalpy_mass(mix1) * 1e3; % [J/kg]
    m1 = mass(mix1); % [kg]
    itMax = self.TN.it_guess_det;
    p2p1_0 = 15;
    % Estimate post-shock state considering h2 = h1 + 3/4 * R * T1 / W1 * (p2 / p1)_0
    mix1.h = (h1 + 3/4 * R0 * T1 / W1 * p2p1_0) * m1 * 1e-3; % [kJ]
    self.PD.ProblemType = 'HP';
    mix2 = equilibrate(self, mix1, p2p1_0 * mix1.p);
    cp2 = cp_mass(mix2) * 1e3; % [J/kg-K]
    gamma2_s = adiabaticIndex_sound(mix2); % [-]
    W2 = meanMolecularWeight(mix2) * 1e-3; % [kg]
    % Initialization
    P_0 = p2p1_0; % [-]
    T_0 = temperature(mix2) / T1; % [-]
    P = P_0;
    T = T_0;

    % Debug
    % fprintf('T_guess,0 = %.4g\n', T)
    % fprintf('\nP_guess,0 = %.4g\n', P)

    it = 0; STOP = 1;

    while STOP > self.TN.tol_shocks && it < itMax
        % Update iteration
        it = it + 1;
        % Update ratios
        P_old = P;
        T_old = T;
        % Compute ratios
        alpha = 1 / T * W2 / W1;
        P = (1 + gamma2_s) / (2 * gamma2_s * alpha) * (1 + sqrt(1 - 4 * gamma2_s * alpha / (1 + gamma2_s)^2));
        r = alpha * P;
        T = T_0 - 3/4 * R0 / (W1 * cp2) * P_0 + R0 * gamma2_s / (2 * W1 * cp2) * (r^2 - 1) / r * P;
        % Compute STOP criteria
        STOP = norm([P - P_old, T - T_old]);
    end

    % Debug
    % fprintf('T_guess,%d = %.4g\n', it, T)
    % fprintf('\nP_guess,%d = %.4g\n', it, P)
end
