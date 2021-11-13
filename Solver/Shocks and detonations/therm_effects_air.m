function r = therm_effects_air(varargin)
%     load test_O2.mat
%     load test.mat
    % Get RH jump conditions
    self = compute_shock(varargin);
     [r.R, r.P, r.T, r.M1, r.M2, r.Gammas] = get_parameters(self);
%    [R, P, T, M1, M2, Gammas] = get_parameters(self);
%     save test_O2.mat  

    Nl = 1e4;
    Ns = 1e4;
    zeta_l = linspace(0, 1, Nl)';
    zeta_s = linspace(1, 1e5, Ns)'; 
%     zeta = [zeta_l; zeta_s];

%     thetaOfzeta = f_thetaOfzeta(r.R, r.M2, zeta);
    
%     sigma_a = f_sigma_a(r.R, r.M2, r.Gammas);
    sigma_b = f_sigma_b(r.M2, r.Gammas);
    sigma_c = f_sigma_c(r.R, r.M2, r.Gammas);
    
    pi_l1 = f_pi_l1(r.R, r.M2, sigma_b, sigma_c, zeta_l);
    pi_l2 = f_pi_l2(r.R, r.M2, sigma_b, sigma_c, zeta_l);
    pi_s  = f_pi_s(r.R, r.M2, sigma_b, sigma_c, zeta_s);
    
    ka = f_ka(r.M2, zeta_s);
    wa = f_wa(r.M2, zeta_s);
    
    Delta_ua = f_Delta_ua(ka, wa, pi_s);
    Delta_va = f_Delta_va(wa, pi_s);
    
    Omega_1_l = f_Omega_1(r.R, r.M2, zeta_l);
    Omega_1_s = f_Omega_1(r.R, r.M2, zeta_s);
    Omega_2 = f_Omega_2(r.R, r.M2, r.Gammas);
    Delta_Omega_l1 = f_Delta_Omega_l1(pi_l1, Omega_1_l, Omega_2);
    Delta_Omega_l2 = f_Delta_Omega_l2(pi_l2, Omega_1_l);
    Delta_Omega_s  = f_Delta_Omega_s(pi_s, Omega_1_s, Omega_2);
    
    Delta_l = f_Delta(r.M2, zeta_l);
    Delta_s = f_Delta(r.M2, zeta_s);
    
    Delta_u_l1 = f_Delta_u_l1(Delta_Omega_l1, Delta_l);
    Delta_u_l2 = f_Delta_u_l2(Delta_Omega_l2, Delta_l);
    Delta_u_s  = f_Delta_u_s(Delta_Omega_s  , Delta_s);
%     Delta_u    = f_Delta_u(Delta_u_l1, Delta_u_l2, Delta_u_s, zeta);
    
    Delta_v_l1 = f_Delta_v_l1(r.M2, Delta_Omega_l1, Delta_l, zeta_l);
    Delta_v_l2 = f_Delta_v_l2(r.M2, Delta_Omega_l2, Delta_l, zeta_l);
    Delta_v_s  = f_Delta_v_s(r.M2, Delta_Omega_s, Delta_s, zeta_s);
%     Delta_v    = f_Delta_v(Delta_v_l1, Delta_v_l2, Delta_v_s, zeta);
    
    PDF_3D_l = f_PDF_3D(r.R, r.M2, zeta_l);
    PDF_3D_s = f_PDF_3D(r.R, r.M2, zeta_s);
%     PDF_3D = f_PDF_3D(r.R, r.M2, zeta);

    L3Drl = f_L3Drl(Delta_u_l1, Delta_u_l2, PDF_3D_l);
    L3Drs = f_L3Drs(Delta_u_s, PDF_3D_s);
    r.L3Dr = f_L3Dr(L3Drl, L3Drs);
    r.L3Da = f_L3Da(Delta_ua, PDF_3D_s);
    r.L3D = f_L3D(r.L3Dr, r.L3Da);

    T3Drl = f_T3Drl(Delta_v_l1, Delta_v_l2, PDF_3D_l);
    T3Drs = f_T3Drs(Delta_v_s, PDF_3D_s);
    r.T3Dr = f_T3Dr(T3Drl, T3Drs);
    r.T3Da = f_T3Da(Delta_va, PDF_3D_s);
    r.T3D = f_T3D(r.T3Dr, r.T3Da);

    r.K3Dr = 1/3 * (r.L3Dr + 2*r.T3Dr);
    r.K3Da = 1/3 * (r.L3Da + 2*r.T3Da);

    r.K3D = 1/3 * (r.L3D + 2*r.T3D);
end

% NESTED FUNCTIONS
function P0 = f_P0(R)
    P0 = @(gamma) ((gamma + 1) - (gamma - 1) .* R.^(-1)) ./ ((gamma + 1).* R.^(-1) - (gamma - 1));
end

function Gammas = f_Gammas(M1, R, P)
    Gammas =  7/5 * (M1(1:end-1).^2 ./ R(1:end-1).^2) .* ((P(2:end)  - P(1:end-1)) ./ (R(2:end)  - R(1:end-1))).^(-1);
end

function thetaOfzeta = f_thetaOfzeta(R, M2, zeta)
    thetaOfzeta = atan((M2 .* R) ./ (sqrt(1 - M2.^2) .* zeta(1:end)));
end

function sigma_a = f_sigma_a(R, M2, Gammas)
    sigma_a =  (R  ./ (R  - 1)) .* ((1 - Gammas) ./ (2*M2));
end

function sigma_b = f_sigma_b(M2, Gammas)
    sigma_b =  (1 + Gammas) / (2*M2);
end

function sigma_c = f_sigma_c(R, M2, Gammas)
    sigma_c =  (((M2.^2 .* R)) ./ (1 - M2.^2)) .* ((1 - Gammas) / 2);
end

function pi_l1 = f_pi_l1(R, M2, sigma_b, sigma_c, zeta)
    pi_l1 = ((-(1 - R.^-1) .* (sigma_b  .* zeta.^2 - sigma_c)) ./ (zeta.^2.* (1 - zeta.^2) + (sigma_b .* zeta.^2 - sigma_c).^2)) .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2));
end

function pi_l2 = f_pi_l2(R, M2, sigma_b, sigma_c, zeta)
    pi_l2 = (((1 - R.^-1) .* zeta .* sqrt((1 - zeta.^2))) ./ (zeta.^2 .* (1 - zeta.^2) + (sigma_b .* zeta.^2 - sigma_c).^2)) .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2));
end

function pi_s = f_pi_s(R, M2, sigma_b, sigma_c, zeta)
    pi_s = ((-(1 - R.^-1)) ./ (zeta .* sqrt(zeta.^2 - 1) + sigma_b .* zeta.^2 - sigma_c)) .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2));
end

function ka = f_ka(M2, zeta)
    ka = (zeta .* M2  - sqrt(zeta.^2 - 1)) ./ (sqrt(1 - M2.^2));
end

function wa = f_wa(M2, zeta)
    wa = (zeta - M2  .* sqrt(zeta.^2 - 1)) ./ (sqrt(1 - M2.^2));
end

function Delta_ua = f_Delta_ua(ka, wa, pi_s)
    Delta_ua = (ka ./ wa) .* pi_s;
end

function Delta_va = f_Delta_va(wa, pi_s)
    Delta_va = (1 ./ wa) .* pi_s;
end

function Omega_1 = f_Omega_1(R, M2, zeta)
    Omega_1 = R .* (1 + zeta.^2 .* ((1 - M2.^2) ./ (R.^2 .* M2.^2)));
end

function Omega_2 = f_Omega_2(R, M2, Gammas)
    Omega_2 = ((R - 1) .* (1 - Gammas)) ./ (2*M2);
end

function Delta_Omega_l1 = f_Delta_Omega_l1(pi_l1, Omega_1, Omega_2)
    Delta_Omega_l1 = Omega_2 .* pi_l1 + Omega_1;
end

function Delta_Omega_l2 = f_Delta_Omega_l2(pi_l2, Omega_2)
    Delta_Omega_l2 = Omega_2 .* pi_l2;
end

function Delta_Omega_s = f_Delta_Omega_s(pi_s, Omega_1, Omega_2)
    Delta_Omega_s = Omega_2 .* pi_s + Omega_1;
end

function Delta = f_Delta(M2, zeta)
    Delta = 1 + ((1 - M2.^2) ./ M2.^2) .* zeta.^2;
end

function Delta_u_l1 = f_Delta_u_l1(Delta_Omega_l1, Delta)
    Delta_u_l1 = Delta_Omega_l1 ./ Delta;
end

function Delta_u_l2 = f_Delta_u_l2(Delta_Omega_l2, Delta)
    Delta_u_l2 = Delta_Omega_l2 ./ Delta;
end

function Delta_u_s = f_Delta_u_s(Delta_Omega_s, Delta)
    Delta_u_s = Delta_Omega_s ./ Delta;
end

function Delta_u = f_Delta_u(Delta_u_l1, Delta_u_l2, Delta_u_s, zeta)
    j = 0;
    count = 0;
    for i=length(zeta):-1:1
        if zeta(i) > 1
            Delta_u(i, :) = Delta_u_s(end-j, :);
            j = j + 1;
        elseif zeta(i) == 1 && count < 1
            count = count + 1;
            Delta_u(i, :) = Delta_u_s(end-j, :);
        elseif zeta(i) == 1 && count < 2
            count = count + 1;
            Delta_u(i, :) = sqrt((Delta_u_l1(i, :).^2 + Delta_u_l2(i, :).^2));
        else
            Delta_u(i, :) = sqrt((Delta_u_l1(i, :).^2 + Delta_u_l2(i, :).^2));
        end
    end
end

%%%
function Delta_v_l1 = f_Delta_v_l1(M2, Delta_Omega_l1, Delta, zeta)
    Delta_v_l1 = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_l1 ./ Delta);
end

function Delta_v_l2 = f_Delta_v_l2(M2, Delta_Omega_l2, Delta, zeta)
    Delta_v_l2 = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_l2 ./ Delta);
end

function Delta_v_s = f_Delta_v_s(M2, Delta_Omega_s, Delta, zeta)
    Delta_v_s = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_s ./ Delta);
end

function Delta_v = f_Delta_v(Delta_v_l1, Delta_v_l2, Delta_v_s, zeta)
    j = 0;
    count = 0;
    for i=length(zeta):-1:1
        if zeta(i) > 1
            Delta_v(i, :) = Delta_v_s(end-j, :);
            j = j + 1;
        elseif zeta(i) == 1 && count < 1
            count = count + 1;
            Delta_v(i, :) = Delta_v_s(end-j, :);
        elseif zeta(i) == 1 && count < 2
            count = count + 1;
            Delta_v(i, :) = sqrt((Delta_v_l1(i, :).^2 + Delta_v_l2(i, :).^2));
        else
            Delta_v(i, :) = sqrt((Delta_v_l1(i, :).^2 + Delta_v_l2(i, :).^2));
        end
    end
end

function PDF_3D = f_PDF_3D(R, M2, zeta)
    PDF_3D = 3/2 * ((M2.^4 .* R.^4 .* sqrt(1 - M2.^2)) ./ (M2.^2 .* R.^2 + zeta.^2 .* (1 - M2.^2)).^(5/2));
end
%%%

function L3Drl = f_L3Drl(Delta_u_l1, Delta_u_l2, PDF_3D)
    fun = (Delta_u_l1.^2 + Delta_u_l2.^2) .* PDF_3D;
    for i=length(PDF_3D(1, :)):-1:1
        L3Drl(i) = trapz(fun(:, i));
    end
end

function L3Drs = f_L3Drs(Delta_u_s, PDF_3D)
    fun = Delta_u_s.^2 .* PDF_3D;
    for i=length(PDF_3D(1, :)):-1:1
        L3Drs(i) = trapz(fun(:, i));
    end
end

function L3Dr = f_L3Dr(L3Drl, L3Drs)
    L3Dr = L3Drl + L3Drs;
end

function L3Da = f_L3Da(Delta_ua, PDF_3D)
    fun = Delta_ua.^2 .* PDF_3D;
    for i=length(PDF_3D(1, :)):-1:1
        L3Da(i) = trapz(fun(:, i));
    end
end

function L3D = f_L3D(L3Dr, L3Da)
    L3D = L3Dr + L3Da;
end




function T3Drl = f_T3Drl(Delta_v_l1, Delta_v_l2, PDF_3D)
    fun = (Delta_v_l1.^2 + Delta_v_l2.^2 + 3/2) .* PDF_3D;
    for i=length(PDF_3D(1, :)):-1:1
        T3Drl(i) = 0.5 * trapz(fun(:, i));
    end
end

function T3Drs = f_T3Drs(Delta_v_s, PDF_3D)
    fun = (Delta_v_s.^2 + 3/2).* PDF_3D;
    for i=length(PDF_3D(1, :)):-1:1
        T3Drs(i) = 0.5 * trapz(fun(:, i));
    end
end

function T3Dr = f_T3Dr(T3Drl, T3Drs)
    T3Dr = T3Drl + T3Drs;
end

function T3Da = f_T3Da(Delta_va, PDF_3D)
    fun = Delta_va.^2 .* PDF_3D;
    for i=length(PDF_3D(1, :)):-1:1
        T3Da(i) = 0.5 * trapz(fun(:, i));
    end
end

function T3D = f_T3D(T3Dr, T3Da)
    T3D = T3Dr + T3Da;
end


function [R, P, T, M1, M2, Gammas] = get_parameters(self)
    R = cell2vector(self.PS.strP, 'rho') ./ cell2vector(self.PS.strR, 'rho');
    P = cell2vector(self.PS.strP, 'p') ./ cell2vector(self.PS.strR, 'p');
    T = cell2vector(self.PS.strP, 'T') ./ cell2vector(self.PS.strR, 'T');
    M1 = cell2vector(self.PS.strR, 'u') ./ cell2vector(self.PS.strR, 'sound');
    M2 = cell2vector(self.PS.strP, 'v_shock') ./ cell2vector(self.PS.strP, 'sound');
    Gammas = f_Gammas(M1, R, P);
    R = R(1:end-1);
    P = P(1:end-1);
    T = T(1:end-1);
    M1 = M1(1:end-1);
    M2 = M2(1:end-1);
end

function self = compute_shock(varargin)
    % Read list of species considered as products
    if ~isempty(varargin{nargin})
        ListProducts = varargin{1}{1};
    else 
        ListProducts = 'Air_ions';
    end
    if nargin > 2, T = varargin{1}{2}; else, T = 300; end % [K]
    if nargin > 3, p = varargin{1}{3}; else, p = 1; end % [bar]
    self = App(ListProducts);
    self.Misc.FLAG_RESULTS = false;
    % Initial conditions
    self = set_prop(self, 'TR', T, 'pR', 1.01325 * p);
    self.PD.S_Oxidizer = {'O2'};
    self.PD.S_Inert    = {'N2', 'Ar', 'CO2'};
%     self.PD.proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;
    % Additional inputs
%     u1 = logspace(2, 5, 500); u1 = u1(u1<20000); u1 = u1(u1>=360);
    u1 = logspace(2, 5, 2e4); u1 = u1(u1<15000); u1 = u1(u1>=360);
%     u1 = [10000, 12000];
    self = set_prop(self, 'u1', u1, 'phi', self.PD.phi.value(1) * ones(1, length(u1)));
    % Solve problem
    self = SolveProblem(self, 'SHOCK_I');
end