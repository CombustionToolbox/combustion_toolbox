function [P, T, M1, R, Q] = compute_guess_cj(self, str1, phi)
    % Paramenters
    gamma  = str1.gamma;
    a1     = str1.sound;
    P      = 15;
    lambda = 0.2;
    % Compute moles considering COMPLETE combustion
    [N_2, LS] = complete_combustion(self, str1, phi);
    [P, T, M1, R, Q] = body_guess_cj(self, P, N_2, LS, gamma, a1);
    % Compute moles considering INCOMPLETE combustion
    it = 0; itMax = 2; STOP = 1;
    while STOP > 1e-1 && it < itMax
        it = it + 1;
        
        P_0 = P; T_0 = T; M1_0 = M1; R_0 = R; Q_0 = Q;
        
        p2 = P_0 * str1.p;
        T2 = T_0 * str1.T;
        N_2 = equilibrium(self, p2, T2, str1); N_2 = N_2(:, 1)';
        LS = self.S.LS;
        [P, T, M1, R, Q] = body_guess_cj(self, P_0, N_2, LS, gamma, a1);
        P  = P_0  + lambda * (P - P_0);
        T  = T_0  + lambda * (T - T_0);
        M1 = M1_0 + lambda * (real(M1) - M1_0);
        R  = R_0  + lambda * (R - R_0);
        Q  = Q_0  + lambda * (Q - Q_0);
        
        STOP = norm([P - P_0, T - T_0, M1 - M1_0, R - R_0, Q - Q_0]);
    end
end

% NESTED FUNCTIONS
function [P, T, M1, R, Q] = body_guess_cj(self, P, N_2, LS, gamma, a1)
    % Get enthalpy of formation [J/mol] and some parameters
    [hfi_1, hfi_2, N_1, Yi_fuel, W_fuel] = get_hf0_molar(self, N_2, LS);
    % Compute heat release - molar base [J/mol]
    q_molar = - sum(N_2 .* hfi_2) + sum(N_1 .* hfi_1);
    % Compute heat release - molar base [J/kg_mixture]
    q = q_molar * Yi_fuel / W_fuel;
    % Compute dimensionless heat release
    Q = (gamma^2 - 1) / (2*a1^2) * q;
    % Compute minimum upstream Mach number (CJ condition)
    Mcj = sqrt(1 + Q) + sqrt(Q);
    % Compute jump relations and upstream Mach number (M1 >= Mcj)
    x = compute_M1_R1(P, gamma, Q, Mcj);
    M1 = x(1); R = x(2);
    % Check M1 >= Mcj
    [P, R, M1] = check_CJ_condition(self, P, R, gamma, Q, M1, Mcj);
    T = P / R;
end

function [hfi_1, hfi_2, N_1, Yi_fuel, W_fuel] = get_hf0_molar(self, N_2, LS_2)
    LS_1 = [self.PD.S_Fuel, self.PD.S_Oxidizer, self.PD.S_Inert];
    N_1  = [self.PD.N_Fuel, self.PD.N_Oxidizer, self.PD.N_Inert];
    
    hfi_1 = loop_hf0_molar(N_1, LS_1, self.strThProp);
    hfi_2 = loop_hf0_molar(N_2, LS_2, self.strThProp);
    
    R = self.PD.R_Fuel + self.PD.R_Oxidizer + self.PD.R_Inert;
    mi = sum(R(:,11)) * 1e-3;
    Yi = R(:,11) ./ mi * 1e-3;
    Yi_fuel = sum(Yi(self.PD.R_Fuel(:, 1) > 0));
    W_fuel  = self.PS.strR_Fuel.W * 1e-3; % [kg/mol]
end

function hfi = loop_hf0_molar(N, LS, strThProp)
    for i = length(N):-1:1
        hfi(i) = strThProp.(LS{i}).hf; % [J/mol]
    end
    hfi(isnan(hfi)) = 0;
    hfi(isinf(hfi)) = 0;
end

function [P, R, M1] = check_CJ_condition(self, P, R, gamma, Q, M1, Mcj)
    if M1 < Mcj || strcmpi(self.PD.ProblemType, 'DET')
        M1 = Mcj;
        P  = compute_P(gamma, Q, M1);
        R  = compute_R(gamma, Q, M1);
    end
end

function P = compute_P(gamma, Q, M1)
    P = real((gamma * M1^2 + 1 + gamma * ((M1^2 - 1)^2 - 4*Q*M1^2)^(1/2)) / (gamma + 1));
end

function R = compute_R(gamma, Q, M1)
    R = real(((gamma + 1) * M1^2) / (gamma * M1^2 + 1 - ((M1^2 - 1)^2 - 4*Q*M1^2)^(1/2)));
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