function [P, T, M1, R, Q, STOP] = det_compute_guess(self, mix1, phi, overdriven)
    % Obtain guess of the jump conditions for a Chapman-Jouguet detonation.
    % Only valid if the mixture have CHON. It computes the guess assuming
    % first a complete combustion, next it recomputes assuming an incomplete
    % combustion from the composition obtained in the previous step.
    %
    % Args:
    %     self (struct):      Data of the mixture, conditions, and databases
    %     mix1 (struct):      Properties of the mixture in the pre-shock state
    %     phi (float):        Equivalence ratio [-]
    %     overdriven (float): Overdriven ratio [-] respect to the sound velocity of the mixture
    %
    % Returns:
    %     Tuple containing
    %
    %     * P (float):        Pressure ratio [-]
    %     * T (float):        Temperature ratio [-]
    %     * M1 (float):       Pre-shock Mach number [-]
    %     * R (float):        Density ratio [-]
    %     * Q (float):        Dimensionless Heat release []
    %     * STOP (float):     Relative error [-]

    % Parameters
    gamma1 = mix1.gamma;
    a1 = mix1.sound;
    DeltaQ = 0;
    itMax = self.TN.it_guess_det;
    % Compute moles considering COMPLETE combustion
    [N_2_cc, LS] = complete_combustion(self, mix1, phi);
    [hfi_1, N_1, Yi_fuel, W_fuel] = compute_hfi_1_molar(self);
    [P, T, M1, M2, R, Q] = body_guess_cj(self, N_2_cc, LS, gamma1, DeltaQ, overdriven, a1, hfi_1, N_1, Yi_fuel, W_fuel);
    % Compute moles considering INCOMPLETE combustion
    [~, ind_a, ind_b] = intersect(self.S.LS, LS);
    N_2 = zeros(1, length(self.S.LS));
    N_2(ind_a) = N_2_cc(ind_b);

    LS = self.S.LS;
    W1 = mix1.W;
    lambda = 2 - sqrt(1 + M2^2);

    %     % Debug
    %     fprintf('T_guess,0 = %.4g\n', T)
    %     fprintf('\nP_guess,0 = %.4g\n', P)

    it = 0; STOP = 1;

    while STOP > self.TN.tol_shocks && it < itMax
        % Update iteration
        it = it + 1;
        % Assign initial values
        P_0 = P; T_0 = T; N_2_0 = N_2;
        % Compute post-shock pressure and temperature
        p2 = P_0 * mix1.p;
        T2 = T_0 * mix1.T;
        % Compute post-shock compostion (moles)
        N_2 = equilibrium_gibbs(self, p2, T2, mix1, []); N_2 = N_2(:, 1)';
        % Apply correction
        N_2 = N_2_0 + lambda .* (N_2 - N_2_0);
        % Compute Mean molecular weight (MW)
        W2 = compute_W(N_2, LS, self.DB);
        % Compute adiabatic index
        gamma2 = compute_gamma(N_2, T * mix1.T, LS, self.DB);
        gamma2 = gamma1 + lambda * (gamma2 - gamma1);
        % Compute ratios and dimensionless heat release
        [P, T, M1, ~, R, Q] = body_guess_cj(self, N_2, LS, gamma1, DeltaQ, overdriven, a1, hfi_1, N_1, Yi_fuel, W_fuel);
        T = T * W2 / W1;
        % Update ratios
        P = P_0 + lambda * (P - P_0);
        T = T_0 + lambda * (T - T_0);
        % Compute correction dimensionless heat release
        DeltaQ = (gamma1 + 1) / (2 * gamma1 * (gamma1 - 1)) * T * (gamma2 - gamma1);
        % Compute STOP criteria
        STOP = norm([P - P_0, T - T_0]);
    end

    % Debug
    % fprintf('T_guess,%d = %.4g\n', it, T)
    % fprintf('\nP_guess,%d = %.4g\n', it, P)
end

% SUB-PASS FUNCTIONS
function [P, T, M1, M2, R, Q] = body_guess_cj(self, N_2, LS, gamma, DeltaQ, overdriven, a1, hfi_1, N_1, Yi_fuel, W_fuel)
    % Get enthalpy of formation [J/mol] and some parameters
    hfi_2 = compute_hfi_molar(N_2, LS, self.DB);
    % Compute heat release - molar base [J/mol]
    q_molar =- dot(N_2, hfi_2) + dot(N_1, hfi_1);
    % Compute heat release - molar base [J/kg_mixture]
    q = q_molar * Yi_fuel / W_fuel;
    % Compute dimensionless heat release | with correction
    Q = abs((gamma^2 - 1) / (2 * a1^2) * q + DeltaQ);
    % Compute minimum upstream Mach number (CJ condition)
    Mcj = sqrt(1 + Q) + sqrt(Q);
    % Compute jump relations and upstream Mach number (M1 >= Mcj)
    M1 = 1.001 * overdriven * Mcj;
    R = compute_R(gamma, Q, M1);
    P = compute_P(gamma, Q, M1);
    M2 = compute_M2(gamma, Q, M1);
    T = P / R;
end

function [hfi_1, ni, Yi_fuel, W_fuel] = compute_hfi_1_molar(self)
    % Compute properties initial mixture

    % Compute properties matrix
    properties_matrix = self.PD.R_Fuel + self.PD.R_Oxidizer + self.PD.R_Inert;
    properties_matrix(:, 1:self.C.M0.ind_ni - 1) = self.C.M0.value(:, 1:self.C.M0.ind_ni - 1);
    % Get moles
    ni = properties_matrix(:, self.C.M0.ind_ni);
    % Compute enthalpy of formation [J/mol]
    hfi_1 = properties_matrix(:, self.C.M0.ind_hfi) * 1e3;
    % Compute molecular weight Fuel mixture [kg/mol]
    W_fuel = MolecularWeight(self.PS.strR_Fuel) * 1e-3;
    % Compute mass mixture [kg]
    mi = dot(properties_matrix(:, self.C.M0.ind_W) * 1e-3, ni);
    % Compute mass fractions
    Yi = (ni .* properties_matrix(:, self.C.M0.ind_W) * 1e-3) ./ mi;
    Yi_fuel = sum(Yi(self.PD.R_Fuel(:, self.C.M0.ind_ni) > 0));
end

function hfi = compute_hfi_molar(N, LS, DB)
    % Compute enthalpy of formation in molar basis
    for i = length(N):-1:1
        hfi(i) = DB.(LS{i}).hf; % [J/mol]
    end

    hfi(isnan(hfi)) = 0;
    hfi(isinf(hfi)) = 0;
end

function gamma = compute_gamma(N, T, LS, DB)
    % Compute adiabatic index
    for i = length(N):-1:1
        cP(i) = species_cP(LS{i}, T, DB); % [J/kg-K]
        cV(i) = species_cV(LS{i}, T, DB); % [J/kg-K]
    end

    gamma = dot(N, cP) / dot(N, cV);
end

function W = compute_W(N, LS, DB)
    % Compute Mean molecular Weight (MW)
    for i = length(N):-1:1
        W(i) = DB.(LS{i}).mm; % [g/mol]
    end

    N_total = sum(N);
    W = dot(N / N_total, W);
end

function [P, R, M1, M2] = check_CJ_condition(self, P, R, gamma, Q, M1, M2, Mcj)
    % Check Chapman-Jouguet condition
    if M1 < Mcj || strcmpi(self.PD.ProblemType, 'DET')
        M1 = 1.001 * Mcj;
        P = compute_P(gamma, Q, M1);
        R = compute_R(gamma, Q, M1);
        M2 = compute_M2(gamma, Q, M1);
    end

end

function P = compute_P(gamma, Q, M1)
    % Compute pressure ratio (p2 / p1)
    P = real((gamma * M1^2 + 1 + gamma * ((M1^2 - 1)^2 - 4 * Q * M1^2)^(1/2)) / (gamma + 1));
end

function R = compute_R(gamma, Q, M1)
    % Compute density ratio (rho2 / rho1)
    R = real(((gamma + 1) * M1^2) / (gamma * M1^2 + 1 - ((M1^2 - 1)^2 - 4 * Q * M1^2)^(1/2)));
end

function M2 = compute_M2(gamma, Q, M1)
    % Compute post-shock Mach number
    M2 = ((gamma * M1^2 + 1 - ((M1^2 - 1)^2 - 4 * Q * M1^2)^(1/2)) / (gamma * M1^2 + 1 + gamma * ((M1^2 - 1)^2 - 4 * Q * M1^2)^(1/2)))^(1/2);
end

function x = compute_M1_R1(P, gamma, Q, Mcj)
    % Compute pre-shock Mach number
    Rcj = compute_R(gamma, Q, Mcj);
    options = optimset('Display', 'off');
    x = real(fsolve(@(x) root2d(x, P, gamma, Q), [Mcj, Rcj], options));
end

function F = root2d(x, P, gamma, Q)
    % Rankine-Hugoniot relations as a function of pre-shock Mach number (M1) and denstiy ratio (R) 
    % F == lhs - rhs = 0
    % x(1) == M1; x(2) == R
    F(1) = x(2) - ((gamma + 1) * x(1)^2) / (gamma * x(1)^2 + 1 - ((x(1)^2 - 1)^2 - 4 * Q * x(1)^2)^(1/2));
    F(2) = x(1) - (P - 1) * (1 - x(2)^(-1)) / gamma;
end

function F = root2d_test(x, P, gamma, R)
    % Rankine-Hugoniot relations as a function of pre-shock Mach number (M1) and dimensionless heat release (Q) 
    % F == lhs - rhs = 0
    % x(1) == M1; x(2) == Q
    F(1) = R - ((gamma + 1) * x(1)^2) / (gamma * x(1)^2 + 1 - ((x(1)^2 - 1)^2 - 4 * x(2) * x(1)^2)^(1/2));
    F(2) = x(1) - (P - 1) * (1 - R^(-1)) / gamma;
end
