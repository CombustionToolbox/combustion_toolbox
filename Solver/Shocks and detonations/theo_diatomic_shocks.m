function [R, P, M1, M2, K3D] = theo_diatomic_shocks(varargin)
    global Infty steps
    Infty = 1e6; steps = 1e6;
    % Get parameters
    [T, alpha, Br, theta_r, theta_v, theta_d] = get_parameters(varargin);
    
    NT = length(T);
    % Definitions general case (dissociation)
    epsilon_v = f_epsilon_v(T, theta_v);
    alpha_R   = f_alpha_R(alpha);
    alpha_T   = f_alpha_T(alpha_R, T, theta_v, theta_d);
    
    R = f_R(alpha, T, Br, theta_v, theta_d);
    P = f_P(alpha, T, R);
    c2qc1sq = 5/7 .* T .* (1 + alpha + alpha_R + ((1 + alpha + alpha_T) .* (2.*(1 + alpha) - alpha_R .* (1 - 2.*epsilon_v + 2.*theta_d./T))) ./ (5 + alpha + 2.*(1 - alpha) .* epsilon_v.^2 .* exp(theta_v./T) + alpha_T .* (1 - 2.*epsilon_v + 2.*theta_d./T)));
    M1 = f_M1(R, P); 
    M2 = f_M2(M1, R, c2qc1sq); 
    P0 = @(gamma) f_P0(R); 
    
    R_per = f_R(alpha, T - 1/1000, Br, theta_v, theta_d);
    P_per = f_P(alpha, T - 1/1000, R_per);
    Gammas = f_Gammas(M1, R, P, R_per, P_per);
    
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
    
    L3Drl = f_L3Drl(Delta_u_l1, Delta_u_l2, PDF_3D, NT);
    L3Drs = f_L3Drs(Delta_u_s, PDF_3D, NT);
    L3Dr  = f_L3Dr(L3Drl, L3Drs);
    L3Da  = f_L3Da(Delta_ua, PDF_3D, NT);
    L3D   = f_L3D(L3Dr, L3Da);
    
    T3Drl = f_T3Drl(Delta_v_l1, Delta_v_l2, PDF_3D, NT);
    T3Drs = f_T3Drs(Delta_v_s, PDF_3D, NT);
    T3Dr  = f_T3Dr(T3Drl, T3Drs);
    T3Da  = f_T3Da(Delta_va, PDF_3D, NT);
    T3D   = f_T3D(T3Dr, T3Da);
    
    K3D = f_K3D(L3D, T3D);
    
%     % Definitions ideal case (no dissociation / frozen chemistry)
    R_ideal  = @(T) 3*(1 - T^(-1)) * (1 + sqrt((1 + T^(-1)) / (9*(1 - T^(-1))^2)));
%     P_ideal  = @(T) R_ideal(T) * T;
%     c2qc1sq_ideal = @(T) T;
%     M1_ideal = f_M1(R_ideal, P_ideal);
%     M2_ideal = f_M2(M1_ideal, R_ideal, c2qc1sq_ideal);
%     P0_ideal = f_P0(R_ideal);
%     Gammas_ideal = f_Gammas(M1_ideal, R_ideal, P_ideal);
%     Gammas_ideal2 = @(T) 1 / M1_ideal(T)^2;
%     thetaOfzeta_ideal = @(T, zeta) arctan((M2_ideal(T) * R_ideal(T)) / (sqrt(1 - M2_ideal(T)^2) * zeta));
%     
%     % Definitions strong limit (complete dissociation)
%     R_strong = @(T) (-6 + 8*T + 2*theta_d + sqrt(8*T + (-6 + 8*T + 2*theta_d)^2)) / (4*T);
%     P_strong = @(T) 2*R_strong(T) * T;
%     c2qc1sq_strong = @(T) 50/21 * T;
%     M1_ideal = f_M1(R_strong, P_strong);
%     M2_strong = f_M2(M1_strong, R_strong, c2qc1sq_strong);
%     Gammas_strong = f_Gammas(M1_strong, R_strong, P_strong);
end

% NESTED FUNCTIONS
function epsilon_v = f_epsilon_v(T, theta_v)
    epsilon_v = (theta_v./T) ./ (exp(theta_v./T) - 1);
end

function alpha_R = f_alpha_R(alpha)
    alpha_R = -((alpha .* (1 - alpha))/(2 - alpha));
end

function alpha_T = f_alpha_T(alpha_R, T, theta_v, theta_d)
    alpha_T = -alpha_R .* (1/2 + (theta_d./T) .* ((1 - (1 + theta_v./theta_d) .* exp(-theta_v./T)) ./ (1 - exp(-theta_v./T))));
end

function R = f_R(alpha, T, Br, theta_v, theta_d)
    R = ((1 - alpha) ./ alpha.^2) .* Br .* sqrt(T) .* (1 - exp(-theta_v ./ T)) .* exp(-theta_d ./ T);
end

function P = f_P(alpha, T, R)
    P = (1 + alpha) .* R .* T;
end

function P0 = f_P0(R)
    P0 = ((gamma + 1) - (gamma - 1) * R^(-1)) / ((gamma + 1) * R^(-1) - (gamma - 1));
end

function M1 = f_M1(R, P)
    % Incident Mach number
    M1 = sqrt(5/7 * (P - 1) ./ (1 - R.^(-1)));
end

function M2 = f_M2(M1, R, c2qc1sq)
    % Incident Mach number
    M2 = (M1 ./ R) ./ sqrt(c2qc1sq);
end

function Gammas = f_Gammas(M1, R, P, R_per, P_per)
    Gammas = (7/5) .* (M1.^2 ./ R.^2) .* ((P - P_per)./(R - R_per)).^(-1);
end

function thetaOfzeta = f_thetaOfzeta(R, M2)
    thetaOfzeta = @(zeta) atan((M2 .* R) ./ (sqrt(1 - M2.^2) .* zeta));
end

function sigma_a = f_sigma_a(R, M2, Gammas)
    sigma_a = (R ./ (R - 1)) .* ((1 - Gammas) ./ (2*M2));
end

function sigma_b = f_sigma_b(M2, Gammas)
    sigma_b = (1 + Gammas) ./ (2*M2);
end

function sigma_c = f_sigma_c(R, M2, Gammas)
    sigma_c = (((M2.^2 .* R)) ./ (1 - M2.^2)) .* ((1 - Gammas) ./ 2);
end

function pi_l1 = f_pi_l1(R, M2, sigma_b, sigma_c)
    pi_l1 = @(zeta) ((-(1 - R.^(-1)) .* (sigma_b .* zeta.^2 - sigma_c)) ./ (zeta.^2 .* (1 - zeta.^2) + (sigma_b .* zeta.^2 - sigma_c).^2)) .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2));
end

function pi_l2 = f_pi_l2(R, M2, sigma_b, sigma_c)
    pi_l2 = @(zeta) (((1 - R.^(-1)) .* zeta .* sqrt((1 - zeta.^2))) ./ (zeta.^2 .* (1 - zeta.^2) + (sigma_b .* zeta.^2 - sigma_c).^2)) .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2));
end

function pi_s = f_pi_s(R, M2, sigma_b, sigma_c)
    pi_s = @(zeta) ((-(1 - R.^(-1))) ./ (zeta .* sqrt(zeta.^2 - 1) + sigma_b .* zeta.^2 - sigma_c)) .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2));
end

function ka = f_ka(M2)
    ka = @(zeta) (zeta .* M2 - sqrt(zeta.^2 - 1)) ./ (sqrt(1 - M2.^2));
end

function wa = f_wa(M2)
    wa = @(zeta) (zeta - M2 .* sqrt(zeta.^2 - 1)) ./ (sqrt(1 - M2.^2));
end

function Delta_ua = f_Delta_ua(ka, wa, pi_s)
    Delta_ua = @(zeta) (ka(zeta) ./ wa(zeta)) .* pi_s(zeta);
end

function Delta_va = f_Delta_va(wa, pi_s)
    Delta_va = @(zeta) (1 ./ wa(zeta)) .* pi_s(zeta);
end

function Omega_1 = f_Omega_1(R, M2)
    Omega_1 = @(zeta) R .* (1 + zeta.^2 .* ((1 - M2.^2) ./ (R.^2 .* M2.^2)));
end

function Omega_2 = f_Omega_2(R, M2, Gammas)
    Omega_2 = @(zeta) ((R - 1) .* (1 - Gammas)) ./ (2*M2);
end

function Delta_Omega_l1 = f_Delta_Omega_l1(pi_l1, Omega_1, Omega_2)
    Delta_Omega_l1 = @(zeta) Omega_2(zeta) .* pi_l1(zeta) + Omega_1(zeta);
end

function Delta_Omega_l2 = f_Delta_Omega_l2(pi_l2, Omega_2)
    Delta_Omega_l2 = @(zeta) Omega_2(zeta) .* pi_l2(zeta);
end

function Delta_Omega_s = f_Delta_Omega_s(pi_s, Omega_1, Omega_2)
    Delta_Omega_s = @(zeta) Omega_2(zeta) .* pi_s(zeta) + Omega_1(zeta);
end

function Delta = f_Delta(M2)
    Delta = @(zeta) 1 + ((1 - M2.^2) ./ M2.^2) .* zeta.^2;
end

function Delta_u_l1 = f_Delta_u_l1(Delta_Omega_l1, Delta)
    Delta_u_l1 = @(zeta) Delta_Omega_l1(zeta) ./ Delta(zeta);
end

function Delta_u_l2 = f_Delta_u_l2(Delta_Omega_l2, Delta)
    Delta_u_l2 = @(zeta) Delta_Omega_l2(zeta) ./ Delta(zeta);
end

function Delta_u_s = f_Delta_u_s(Delta_Omega_s, Delta)
    Delta_u_s = @(zeta) Delta_Omega_s(zeta) ./ Delta(zeta);
end

function Delta_u = f_Delta_u(Delta_u_l1, Delta_u_l2, Delta_u_s, zeta)
    if zeta > 1
        Delta_u = Delta_u_s(zeta);
    else
        Delta_u = sqrt((Delta_u_l1(zeta).^2 + Delta_u_l2(zeta).^2));
    end
end

function Delta_v_l1 = f_Delta_v_l1(M2, Delta_Omega_l1, Delta)
    Delta_v_l1 = @(zeta) zeta .* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_l1(zeta) ./ Delta(zeta));
end

function Delta_v_l2 = f_Delta_v_l2(M2, Delta_Omega_l2, Delta)
    Delta_v_l2 = @(zeta) zeta .* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_l2(zeta) ./ Delta(zeta));
end

function Delta_v_s = f_Delta_v_s(M2, Delta_Omega_s, Delta)
    Delta_v_s = @(zeta) zeta .* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_s(zeta) ./ Delta(zeta));
end

function Delta_v = f_Delta_v(Delta_v_l1, Delta_v_l2, Delta_v_s, zeta)
    if zeta > 1
        Delta_v = @(zeta) Delta_v_s(zeta);
    else
        Delta_v = @(zeta) sqrt((Delta_v_l1(zeta).^2 + Delta_v_l2(zeta).^2));
    end
end
% PDF
function PDF_3D = f_PDF_3D(R, M2)
    PDF_3D = @(zeta) 3/2 .* ((M2.^4 .* R.^4 .* sqrt(1 - M2.^2)) ./ (M2.^2 .* R.^2 + zeta.^2 .* (1 - M2.^2)).^(5/2));
end
% Aux function for integration
function value = aux_int(fun, z0, zf, NT)
    global steps
    zeta = linspace(z0, zf, steps);
    for i=length(zeta):-1:1
        int_fun(i, :) = fun(zeta(i));
    end
    for i=NT:-1:1
        value(i) = trapz(zeta, int_fun(:, i));
    end
end
% Longitudinal
function L3Drl = f_L3Drl(Delta_u_l1, Delta_u_l2, PDF_3D, NT)
    fun = @(zeta) (Delta_u_l1(zeta).^2 + Delta_u_l2(zeta).^2) .* PDF_3D(zeta);
    L3Drl = aux_int(fun, 0, 1, NT);
end

function L3Drs = f_L3Drs(Delta_u_s, PDF_3D, NT)
    global Infty
    fun = @(zeta) Delta_u_s(zeta).^2 .* PDF_3D(zeta);
    L3Drs = aux_int(fun, 1, Infty, NT);
end

function L3Dr = f_L3Dr(L3Drl, L3Drs)
    L3Dr = L3Drl + L3Drs;
end

function L3Da = f_L3Da(Delta_ua, PDF_3D, NT)
    global Infty
    fun = @(zeta) Delta_ua(zeta).^2 .* PDF_3D(zeta);
    L3Da = aux_int(fun, 1, Infty, NT);
end

function L3D = f_L3D(L3Dr, L3Da)
    L3D = L3Dr + L3Da;
end
% Transversal
function T3Drl = f_T3Drl(Delta_v_l1, Delta_v_l2, PDF_3D, NT)
    fun = @(zeta) 1/2 * (Delta_v_l1(zeta).^2 + Delta_v_l2(zeta).^2) .* PDF_3D(zeta);
    T3Drl = aux_int(fun, 0, 1, NT);
end

function T3Drs = f_T3Drs(Delta_v_s, PDF_3D, NT)
    global Infty
    fun = @(zeta) 1/2 * Delta_v_s(zeta).^2 .* PDF_3D(zeta);
    T3Drs = aux_int(fun, 1, Infty, NT);
end

function T3Dr = f_T3Dr(T3Drl, T3Drs)
    T3Dr = T3Drl + T3Drs;
end

function T3Da = f_T3Da(Delta_va, PDF_3D, NT)
    global Infty
    fun = @(zeta) 1/2 * Delta_va(zeta).^2 .* PDF_3D(zeta);
    T3Da = aux_int(fun, 1, Infty, NT);
end

function T3D = f_T3D(T3Dr, T3Da)
    T3D = T3Dr + T3Da;
end

% Kinetic Energy (TKE)
function K3D = f_K3D(L3D, T3D)
    K3D = 1/3 * (L3D + 2*T3D);
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
    
    alpha = zeros(1, length(T));
    
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
    
%     Br = vpa(1e6); theta_v = vpa(10); theta_d = vpa(100); % Generic case
    
    Den = @(T) vpa(((1 - x) / x^2) * Br * sqrt(T) * (1 - exp(-theta_v / T)) * exp(-theta_d / T));
    Tem = @(T, R) vpa((6 - R^(-1) - 2*x * theta_d - 2*(1 - x) * theta_v / (exp(theta_v/T) - vpa(1))) / (2*(x + 3) - R * (1 + x)));
    
    eqs = @(T) Tem(T, Den(T)) == T;
    f_alpha = @(T, x0) abs(vpasolve(eqs(T), x, vpa(x0)));
    
    T = vpa(T);
    
    value = f_alpha(T(end), [vpa(0), vpa(1)]);
    alpha(end) = value(1);
    for i = length(T)-1:-1:1
        value = f_alpha(T(i), [vpa(0), vpa(1)]);
        alpha(i) = value(1);
    end
    
    T  = double(T);
    Br = double(Br);
    theta_r = double(theta_r);
    theta_v = double(theta_v);
    theta_d = double(theta_d);
end
