function [P, T, M1, R, Q, STOP] = compute_guess_det(self, mix1, phi, overdriven)
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
    %     - P (float):          Pressure ratio [-]
    %     - T (float):          Temperature ratio [-]
    %     - M1 (float):         Pre-shock Mach number [-]
    %     - R (float):          Density ratio [-]
    %     - Q (float):          Dimensionless Heat release []
    %     - STOP (float):       Relative error [-] 

    % Paramenters
    gamma1 = mix1.gamma;
    a1     = mix1.sound;
    DeltaQ = 0;
    itMax  = 5;
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
    P = P/3; T = T/3; % complete reaction gives too high values
    it = 0; STOP = 1;
    while STOP > self.TN.tol_shocks && it < itMax
        it = it + 1;
        
        P_0 = P; T_0 = T; N_2_0 = N_2;
        
        p2 = P_0 * mix1.p;
        T2 = T_0 * mix1.T;
        
        N_2 = equilibrium(self, p2, T2, mix1); N_2 = N_2(:, 1)';
        N_2 = N_2_0  + lambda .* (N_2 - N_2_0);
        
        W2 = compute_W(N_2, LS, self.DB);

        gamma2 = compute_gamma(N_2, T * mix1.T, LS, self.DB);
        gamma2 = gamma1 + lambda * (gamma2 - gamma1);

        [P, T, M1, ~, R, Q] = body_guess_cj(self, N_2, LS, gamma1, DeltaQ, overdriven, a1, hfi_1, N_1, Yi_fuel, W_fuel);
        T = T * W2 / W1;
        
        P  = P_0  + lambda * (P - P_0);
        T  = T_0  + lambda * (T - T_0);
        DeltaQ = (gamma1 + 1) / (2*gamma1 * (gamma1 - 1)) * T * (gamma2 - gamma1);

        STOP = norm([P - P_0, T - T_0]);
    end
end

% SUB-PASS FUNCTIONS
function [P, T, M1, M2, R, Q] = body_guess_cj(self, N_2, LS, gamma, DeltaQ, overdriven, a1, hfi_1, N_1, Yi_fuel, W_fuel)
    % Get enthalpy of formation [J/mol] and some parameters
    hfi_2 = compute_hfi_molar(N_2, LS, self.DB);
    % Compute heat release - molar base [J/mol]
    q_molar = - sum(N_2 .* hfi_2) + sum(N_1 .* hfi_1);
    % Compute heat release - molar base [J/kg_mixture]
    q = q_molar * Yi_fuel / W_fuel;
    % Compute dimensionless heat release | with correction
    Q = abs((gamma^2 - 1) / (2*a1^2) * q + DeltaQ); 
    % Compute minimum upstream Mach number (CJ condition)
    Mcj = sqrt(1 + Q) + sqrt(Q);
    % Compute jump relations and upstream Mach number (M1 >= Mcj)
    M1 = 1.001 * overdriven * Mcj;
    R  = compute_R(gamma,  Q, M1);
    P  = compute_P(gamma,  Q, M1);
    M2 = compute_M2(gamma, Q, M1);
    T  = P / R;
end

function [hfi_1, N_1, Yi_fuel, W_fuel] = compute_hfi_1_molar(self)
    LS_1 = [self.PD.S_Fuel, self.PD.S_Oxidizer, self.PD.S_Inert];
    N_1  = [self.PD.N_Fuel, self.PD.N_Oxidizer, self.PD.N_Inert];
    hfi_1 = compute_hfi_molar(N_1, LS_1, self.DB);
    R = self.PD.R_Fuel + self.PD.R_Oxidizer + self.PD.R_Inert;
    mi = sum(R(:,11)) * 1e-3;
    Yi = R(:,11) ./ mi * 1e-3;
    Yi_fuel = sum(Yi(self.PD.R_Fuel(:, 1) > 0));
    W_fuel  = self.PS.strR_Fuel.W * 1e-3; % [kg/mol]
end

function hfi = compute_hfi_molar(N, LS, strThProp)
    for i = length(N):-1:1
        hfi(i) = strThProp.(LS{i}).hf; % [J/mol]
    end
    hfi(isnan(hfi)) = 0;
    hfi(isinf(hfi)) = 0;
end

function gamma = compute_gamma(N, T, LS, strThProp)
    for i = length(N):-1:1
        cP(i) = species_cP(LS{i}, T, strThProp); % [J/kg-K]
        cV(i) = species_cV(LS{i}, T, strThProp); % [J/kg-K]
    end
    gamma = sum(N .* cP) / sum(N .* cV);
end

function W = compute_W(N, LS, strThProp)
    for i = length(N):-1:1
        W(i) = strThProp.(LS{i}).mm; % [g/mol]
    end
    N_total = sum(N);
    W = sum(N/N_total .* W);
end

function [P, R, M1, M2] = check_CJ_condition(self, P, R, gamma, Q, M1, M2, Mcj)
    if M1 < Mcj || strcmpi(self.PD.ProblemType, 'DET')
        M1 = 1.001*Mcj;
        P  = compute_P(gamma, Q, M1);
        R  = compute_R(gamma, Q, M1);
        M2 = compute_M2(gamma, Q, M1);
    end
end

function P = compute_P(gamma, Q, M1)
    P = real((gamma * M1^2 + 1 + gamma * ((M1^2 - 1)^2 - 4*Q*M1^2)^(1/2)) / (gamma + 1));
end

function R = compute_R(gamma, Q, M1)
    R = real(((gamma + 1) * M1^2) / (gamma * M1^2 + 1 - ((M1^2 - 1)^2 - 4*Q*M1^2)^(1/2)));
end

function M2 = compute_M2(gamma, Q, M1)
    M2 = ((gamma * M1^2 + 1 - ((M1^2 - 1)^2 - 4*Q * M1^2)^(1/2)) / (gamma * M1^2 + 1 + gamma * ((M1^2 - 1)^2 - 4*Q * M1^2)^(1/2)))^(1/2);
end

function x = compute_M1_R1(P, gamma, Q, Mcj)
    Rcj = compute_R(gamma, Q, Mcj);   
    options = optimset('Display','off');
    x = real(fsolve(@(x) root2d(x, P, gamma, Q), [Mcj, Rcj], options));
end

function F = root2d(x, P, gamma, Q)
    %    F == lhs - rhs = 0
    % x(1) == M1; x(2) == R
    F(1) = x(2) - ((gamma + 1) * x(1)^2) / (gamma * x(1)^2 + 1 - ((x(1)^2 - 1)^2 - 4*Q*x(1)^2)^(1/2));
    F(2) = x(1) - (P - 1) * (1 - x(2)^(-1)) / gamma;
end

function F = root2d_test(x, P, gamma, R)
    %    F == lhs - rhs = 0
    % x(1) == M1; x(2) == Q
    F(1) = R - ((gamma + 1) * x(1)^2) / (gamma * x(1)^2 + 1 - ((x(1)^2 - 1)^2 - 4*x(2)*x(1)^2)^(1/2));
    F(2) = x(1) - (P - 1) * (1 - R^(-1)) / gamma;
end