function [R, P] = therm_effects_diatomic_sym(varargin)
    % Get parameters
    [T, alpha, Br, theta_r, theta_v, theta_d] = get_parameters(varargin);
    % Definitions general case (dissociation)
    epsilon_v = @(T) (theta_v/T) / (exp(theta_v/T) - 1);
    alpha_R   = @(T) -((alpha(T) .* (1 - alpha(T)))/(2 - alpha(T)));
    alpha_T   = @(T) -alpha_R(T) .* (1/2 + (theta_d./T) .* ((1 - (1 + theta_v./theta_d) .* exp(-theta_v./T)) ./ (1 - exp(-theta_v./T))));
    
    R  = @(T) ((1 - alpha(T)) ./ alpha(T).^2) .* Br .* sqrt(T) .* (1 - exp(-theta_v ./ T)) .* exp(-theta_d ./ T);
    P  = @(T) (1 + alpha(T)) .* R(T) .* T;
    c2qc1sq = @(T) 5/7 * T * (1 + alpha(T) + alpha_R(T) + ((1 + alpha(T) + alpha_T(T)) * (2*(1 + alpha(T)) - alpha_R(T) * (1 - 2*epsilon_v(T) + 2*theta_d/T)))/(5 + alpha(T) + 2*(1 - alpha(T)) * epsilon_v(T)^2 * exp(theta_v/T) + alpha_T(T) * (1 - 2*epsilon_v(T) + 2*theta_d/T)));
    M1 = f_M1(R, P); % f(T)
    M2 = f_M2(M1, R, c2qc1sq); % f(T)
    P0 = f_P0(R); % f(T, gamma)
    Gammas = f_Gammas(M1, R, P);
    thetaOfzeta = f_thetaOfzeta(R, M2);
    
    sigma_a = f_sigma_a(R, M2, Gammas);
    sigma_b = f_sigma_b(M2, Gammas);
    sigma_c = f_sigma_c(R, M2, Gammas);
    
    pi_l1 = f_pi_l1(R, M2, sigma_b, sigma_c);
    pi_l2 = f_pi_l2(R, M2, sigma_b, sigma_c);
    pi_s  = f_pi_s(R, M2, sigma_b, sigma_c);
    
    ka = f_ka(M2);
    wa = f_wa(M2);
    
    Delta_ua = f_Delta_ua(ka, wa, pi_s);
    Delta_va = f_Delta_va(wa, pi_s);
    
    Omega_1 = f_Omega_1(R, M2);
    Omega_2 = f_Omega_2(R, M2, Gammas);
    Delta_Omega_l1 = f_Delta_Omega_l1(pi_l1, Omega_1, Omega_2);
    Delta_Omega_l2 = f_Delta_Omega_l2(pi_l2, Omega_2);
    Delta_Omega_s  = f_Delta_Omega_s(pi_s, Omega_1, Omega_2);
    
    Delta = f_Delta(M2);
    
    Delta_u_l1 = f_Delta_u_l1(Delta_Omega_l1, Delta);
    Delta_u_l2 = f_Delta_u_l2(Delta_Omega_l2, Delta);
    Delta_u_s  = f_Delta_u_s(Delta_Omega_s  , Delta);
    Delta_u    = @(zeta) f_Delta_u(Delta_u_l1, Delta_u_l2, Delta_u_s, zeta);
    
    Delta_v_l1 = f_Delta_v_l1(M2, Delta_Omega_l1, Delta);
    Delta_v_l2 = f_Delta_v_l2(M2, Delta_Omega_l2, Delta);
    Delta_v_s  = f_Delta_v_s(M2, Delta_Omega_s, Delta);
    Delta_v    = @(zeta) f_Delta_v(Delta_v_l1, Delta_v_l2, Delta_v_s, zeta);
    
    PDF_3D = f_PDF_3D(R, M2);
    
    L3Drl = f_L3Drl(Delta_u_l1, Delta_u_l2, PDF_3D);
    % Definitions ideal case (no dissociation / frozen chemistry)
    R_ideal  = @(T) 3*(1 - T^(-1)) * (1 + sqrt((1 + T^(-1)) / (9*(1 - T^(-1))^2)));
    P_ideal  = @(T) R_ideal(T) * T;
    c2qc1sq_ideal = @(T) T;
    M1_ideal = f_M1(R_ideal, P_ideal);
    M2_ideal = f_M2(M1_ideal, R_ideal, c2qc1sq_ideal);
    P0_ideal = f_P0(R_ideal);
    Gammas_ideal = f_Gammas(M1_ideal, R_ideal, P_ideal);
    Gammas_ideal2 = @(T) 1 / M1_ideal(T)^2;
    thetaOfzeta_ideal = @(T, zeta) arctan((M2_ideal(T) * R_ideal(T)) / (sqrt(1 - M2_ideal(T)^2) * zeta));
    
    % Definitions strong limit (complete dissociation)
    R_strong = @(T) (-6 + 8*T + 2*theta_d + sqrt(8*T + (-6 + 8*T + 2*theta_d)^2)) / (4*T);
    P_strong = @(T) 2*R_strong(T) * T;
    c2qc1sq_strong = @(T) 50/21 * T;
    M1_ideal = f_M1(R_strong, P_strong);
    M2_strong = f_M2(M1_strong, R_strong, c2qc1sq_strong);
    Gammas_strong = f_Gammas(M1_strong, R_strong, P_strong);
end

% NESTED FUNCTIONS
function P0 = f_P0(R)
    P0 = @(T, gamma) ((gamma + 1) - (gamma - 1) * R(T)^(-1)) / ((gamma + 1) * R(T)^(-1) - (gamma - 1));
end

function M1 = f_M1(R, P)
    % Incident Mach number
    M1 = @(T) sqrt(5/7 * (P(T) - 1) / (1 - R(T)^(-1)));
end

function M2 = f_M2(M1, R, c2qc1sq)
    % Incident Mach number
    M2 = @(T) (M1(T)/R(T)) / sqrt(c2qc1sq(T));
end

function Gammas = f_Gammas(M1, R, P)
    Gammas = @(T) 7/5 * (M1(T)^2 / R(T)^2) * ((P(T) - P(T - 1/1000))/(R(T) - R(T - 1/1000)))^(-1);
end

function thetaOfzeta = f_thetaOfzeta(R, M2)
    thetaOfzeta = @(T, zeta) atan((M2(T) * R(T)) / (sqrt(1 - M2(T)^2) * zeta));
end

function sigma_a = f_sigma_a(R, M2, Gammas)
    sigma_a = @(T) (R(T) / (R(T) - 1)) * ((1 - Gammas(T)) / (2*M2(T)));
end

function sigma_b = f_sigma_b(M2, Gammas)
    sigma_b = @(T) (1 + Gammas(T)) / (2*M2(T));
end

function sigma_c = f_sigma_c(R, M2, Gammas)
    sigma_c = @(T) (((M2(T)^2 * R(T))) / (1 - M2(T)^2)) * ((1 - Gammas(T)) / 2);
end

function pi_l1 = f_pi_l1(R, M2, sigma_b, sigma_c)
    pi_l1 = @(T, zeta) ((-(1 - R(T).^-1) .* (sigma_b(T) .* zeta.^2 - sigma_c(T))) ./ (zeta.^2 * (1 - zeta.^2) + (sigma_b(T) .* zeta.^2 - sigma_c(T)).^2)) .* (zeta.^2 - (R(T) .* M2(T).^2) ./ (1 - M2(T).^2));
end

function pi_l2 = f_pi_l2(R, M2, sigma_b, sigma_c)
    pi_l2 = @(T, zeta) (((1 - R(T)^-1) * zeta * sqrt((1 - zeta^2))) / (zeta^2 * (1 - zeta^2) + (sigma_b(T) * zeta^2 - sigma_c(T))^2)) * (zeta^2 - (R(T) * M2(T)^2) / (1 - M2(T)^2));
end

function pi_s = f_pi_s(R, M2, sigma_b, sigma_c)
    pi_s = @(T, zeta) ((-(1 - R(T)^-1)) / (zeta * sqrt(zeta^2 - 1) + sigma_b(T) * zeta^2 - sigma_c(T))) * (zeta^2 - (R(T) * M2(T)^2) / (1 - M2(T)^2));
end

function ka = f_ka(M2)
    ka = @(T, zeta) (zeta * M2(T) - sqrt(zeta^2 - 1)) / (sqrt(1 - M2(T)^2));
end

function wa = f_wa(M2)
    wa = @(T, zeta) (zeta - M2(T) * sqrt(zeta^2 - 1)) / (sqrt(1 - M2(T)^2));
end

function Delta_ua = f_Delta_ua(ka, wa, pi_s)
    Delta_ua = @(T, zeta) (ka(T, zeta) / wa(T, zeta)) * pi_s(T, zeta);
end

function Delta_va = f_Delta_va(wa, pi_s)
    Delta_va = @(T, zeta) (1 / wa(T, zeta)) * pi_s(T, zeta);
end

function Omega_1 = f_Omega_1(R, M2)
    Omega_1 = @(T, zeta) R(T) * (1 + zeta^2 * ((1 - M2(T)^2) / (R(T)^2 * M2(T)^2)));
end

function Omega_2 = f_Omega_2(R, M2, Gammas)
    Omega_2 = @(T, zeta) ((R(T) - 1) * (1 - Gammas(T))) / (2*M2(T));
end

function Delta_Omega_l1 = f_Delta_Omega_l1(pi_l1, Omega_1, Omega_2)
    Delta_Omega_l1 = @(T, zeta) Omega_2(T, zeta) * pi_l1(T, zeta) + Omega_1(T, zeta);
end

function Delta_Omega_l2 = f_Delta_Omega_l2(pi_l2, Omega_2)
    Delta_Omega_l2 = @(T, zeta) Omega_2(T, zeta) * pi_l2(T, zeta);
end

function Delta_Omega_s = f_Delta_Omega_s(pi_s, Omega_1, Omega_2)
    Delta_Omega_s = @(T, zeta) Omega_2(T, zeta) * pi_s(T, zeta) + Omega_1(T, zeta);
end

function Delta = f_Delta(M2)
    Delta = @(T, zeta) 1 + ((1 - M2(T)^2) / M2(T)^2) * zeta^2;
end

function Delta_u_l1 = f_Delta_u_l1(Delta_Omega_l1, Delta)
    Delta_u_l1 = @(T, zeta) Delta_Omega_l1(T, zeta) / Delta(T, zeta);
end

function Delta_u_l2 = f_Delta_u_l2(Delta_Omega_l2, Delta)
    Delta_u_l2 = @(T, zeta) Delta_Omega_l2(T, zeta) / Delta(T, zeta);
end

function Delta_u_s = f_Delta_u_s(Delta_Omega_s, Delta)
    Delta_u_s = @(T, zeta) Delta_Omega_s(T, zeta) / Delta(T, zeta);
end

function Delta_u = f_Delta_u(Delta_u_l1, Delta_u_l2, Delta_u_s, zeta)
    if zeta > 1
        Delta_u = @(T, zeta) Delta_u_s(T, zeta);
    else
        Delta_u = @(T, zeta) sqrt((Delta_u_l1(T, zeta)^2 + Delta_u_l2(T, zeta)^2));
    end
end

%%%
function Delta_v_l1 = f_Delta_v_l1(M2, Delta_Omega_l1, Delta)
    Delta_v_l1 = @(T, zeta) zeta * ((sqrt(1 - M2(T)^2)) / M2(T)) * (Delta_Omega_l1(T, zeta) / Delta(T, zeta));
end

function Delta_v_l2 = f_Delta_v_l2(M2, Delta_Omega_l2, Delta)
    Delta_v_l2 = @(T, zeta) zeta * ((sqrt(1 - M2(T)^2)) / M2(T)) * (Delta_Omega_l2(T, zeta) / Delta(T, zeta));
end

function Delta_v_s = f_Delta_v_s(M2, Delta_Omega_s, Delta)
    Delta_v_s = @(T, zeta) zeta * ((sqrt(1 - M2(T)^2)) / M2(T)) * (Delta_Omega_s(T, zeta) / Delta(T, zeta));
end

function Delta_v = f_Delta_v(Delta_v_l1, Delta_v_l2, Delta_v_s, zeta)
    if zeta > 1
        Delta_v = @(T, zeta) Delta_v_s(T, zeta);
    else
        Delta_v = @(T, zeta) sqrt((Delta_v_l1(T, zeta)^2 + Delta_v_l2(T, zeta)^2));
    end
end

function PDF_3D = f_PDF_3D(R, M2)
    PDF_3D = @(T, zeta) 3/2 * ((M2(T)^4 * R(T)^4 * sqrt(1 - M2(T)^2)) / (M2(T)^2 * R(T)^2 + zeta^2 * (1 - M2(T)^2))^(5/2));
end
%%%

function L3Drl = f_L3Drl(Delta_u_l1, Delta_u_l2, PDF_3D)
    fun = @(T, zeta) (Delta_u_l1(T, zeta)^2 + Delta_u_l2(T, zeta)^2) * PDF_3D(T, zeta);
    L3Drl = @(T) integral(@(zeta) fun(T, zeta), 0, 1);
end

function [Tr, Tv, Td, G, m] = get_initial_state(gas)
    Tr = vpa(gas.Tr);
    Tv = vpa(gas.Tv);
    Td = vpa(gas.Td);
    G  = vpa(gas.G);
    m  = vpa(gas.m);
end

function [T, alpha, Br, theta_r, theta_v, theta_d] = get_parameters(varargin)

    gas = varargin{1}{1};
    T   = varargin{1}{2};
    if nargin > 2, T0 = varargin{1}{3}; else, T0 = 300;    end % [K]
    if nargin > 3, p0 = varargin{1}{3}; else, p0 = 101325; end % [Pa]
    
    syms x positive
    digits(250);
    % Get initial state
    [Tr, Tv, Td, G, m] = get_initial_state(gas);
    % Constants and some functions
    T0   = vpa(T0);
    p0   = vpa(p0);
    kb   = vpa(1.38064852  * 1e-23); % Boltzmann's constant
    h    = vpa(6.62607004  * 1e-34); % Planck's constant
    hbar = vpa(1.054571817 * 1e-34); % Reduced Planck's constant
    Rg   = vpa(kb / (2*m));
    r1   = vpa(p0 / (Rg * T0));
    Br   = vpa(m * ((pi * m * kb) / hbar^2)^(3/2) * Tr * G * (sqrt(T0)/r1));
    theta_r = Tr / T0;
    theta_v = Tv / T0;
    theta_d = Td / T0;
    
    Den = @(T) vpa(((1 - x) / x^2) * Br * sqrt(T) * (1 - exp(-theta_v / T)) * exp(-theta_d / T));
    Tem = @(T, R) vpa((6 - R^(-1) - 2*x * theta_d - 2*(1 - x) * theta_v / (exp(theta_v/T) - vpa(1))) / (2*(x + 3) - R * (1 + x)));
    
    eqs = @(T) Tem(T, Den(T)) == T;
    f_alpha = @(T, x0) abs(vpasolve(eqs(T), x, vpa(x0)));
    
    alpha = @(T) double(f_alpha(vpa(T), [vpa(0), vpa(1)]));
    
    Br = double(Br);
    theta_r = double(theta_r);
    theta_v = double(theta_v);
    theta_d = double(theta_d);
end
