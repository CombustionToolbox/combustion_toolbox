function r = therm_effects_air(varargin)
%     load test_O2.mat
%     load CombustionToolboxAir.mat R P T M1 M2 Gammas
%       load thermo_num_air.mat
%     r.R = R; r.P = P; r.T = T; r.M1 = M1; r.M2 = M2; r.Gammas = Gammas;
%     clearvars R P T M1 M2 Gammas;

    % Get RH jump conditions
    self = compute_shock(varargin);
%     [r.R, r.P, r.T, r.M1, r.M2, r.Gammas, r.alpha] = get_parameters(self);
    
    [R, P, T, M1, M2, Gammas] = get_parameters(self);
%     clearvar self varargin
    save thermo_num_air.mat  
%
%     Nl = 1e4;
%     Ns = 1e4;
%     zeta_l = linspace(0, 1, Nl)';
%     zeta_s = linspace(1, 1e5, Ns)'; 
%     zeta = [zeta_l; zeta_s];

%     thetaOfzeta = thetaOfzeta(r.R, r.M2);
    Ndata = length(r.R);
    for i=Ndata:-1:2
        R = r.R(i); M2 = r.M2(i); Gammas = r.Gammas(i);

        r.L3Dr(i) = L3Dr(R, M2, Gammas);
        r.L3Da(i) = L3Da(R, M2, Gammas);
    
        r.T3Dr(i) = T3Dr(R, M2, Gammas);
        r.T3Da(i) = T3Da(R, M2, Gammas);
    end
    
    r.R = r.R(2:end);
    r.P = r.P(2:end);
    r.M1 = r.M1(2:end);
    r.M2 = r.M2(2:end);
    r.Gammas = r.Gammas(2:end);
    r.T = r.T(2:end);
    r.L3Dr = r.L3Dr(2:end); 
    r.L3Da = r.L3Da(2:end);
    r.T3Dr = r.T3Dr(2:end);
    r.T3Da = r.T3Da(2:end);

    r.L3D = r.L3Dr + r.L3Da;
    r.T3D = r.T3Dr + r.T3Da;
    r.K3D = 1/3 * (r.L3D + 2*r.T3D);
    r.Rturb3D = r.R .* sqrt(r.K3D);
    r.RRe3D = (sqrt(r.K3D) .* sqrt((1 + 2 * r.R.^2) / 3)) ./ (r.T.^(0.7));
    r.Anisotropy = 1 - (4 * r.L3D) ./ (3 * r.K3D + r.L3D);
%     r.Anisotropy = 1 - (2 * r.L3D) ./ (r.L3D + r.T3D);
%     % CHECKS
%     j = 6000;
%     zeta = 0.5;
%     R = r.R(j); M2 = r.M2(j); Gammas = r.Gammas(j);
%     check_solution(R, M2, Gammas, zeta)
end

% SUB PASS FUNCTIONS
function check_solution(R, M2, Gammas, zeta)
    fprintf('\nCHECK SOLUTION\n\n');
    fprintf('sigma_a        = %.6f\n', sigma_a(R, M2, Gammas));
    fprintf('sigma_b        = %.6f\n', sigma_b(M2, Gammas));
    fprintf('sigma_c        = %.6f\n', sigma_c(R, M2, Gammas));
    fprintf('\n')
    fprintf('pi_l1          = %.6f\n', pi_l1(R, M2, Gammas, zeta));
    fprintf('pi_l2          = %.6f\n', pi_l2(R, M2, Gammas, zeta));
    fprintf('pi_s           = %.6f\n', pi_s(R, M2, Gammas, zeta));
    fprintf('ka             = %.6f\n', ka(M2, zeta));
    fprintf('wa             = %.6f\n', wa(M2, zeta));
    fprintf('\n')
    fprintf('Delta_ua       = %.6f\n', Delta_ua(R, M2, Gammas, zeta));
    fprintf('Delta_va       = %.6f\n', Delta_va(R, M2, Gammas, zeta));
    fprintf('\n')
    fprintf('Omega_1        = %.6f\n', Omega_1(R, M2, zeta));
    fprintf('Omega_2        = %.6f\n', Omega_2(R, M2, Gammas));
    fprintf('Delta_Omega_l1 = %.6f\n', Delta_Omega_l1(R, M2, Gammas, zeta));
    fprintf('Delta_Omega_l2 = %.6f\n', Delta_Omega_l2(R, M2, Gammas, zeta));
    fprintf('Delta_Omega_s  = %.6f\n', Delta_Omega_s(R, M2, Gammas, zeta));
    fprintf('\n')
    fprintf('Delta          = %.6f\n', Delta(M2, zeta));
    fprintf('\n')   
    fprintf('Delta_u_l1     = %.6f\n', Delta_u_l1(R, M2, Gammas, zeta));
    fprintf('Delta_u_l2     = %.6f\n', Delta_u_l2(R, M2, Gammas, zeta));
    fprintf('Delta_u_s      = %.6f\n', Delta_u_s(R, M2, Gammas, zeta));
    fprintf('Delta_u        = %.6f\n', Delta_u(R, M2, Gammas, zeta));
    fprintf('\n')
    fprintf('Delta_v_l1     = %.6f\n', Delta_v_l1(R, M2, Gammas, zeta));
    fprintf('Delta_v_l2     = %.6f\n', Delta_v_l2(R, M2, Gammas, zeta));
    fprintf('Delta_v_s      = %.6f\n', Delta_v_s(R, M2, Gammas, zeta));
    fprintf('Delta_v        = %.6f\n', Delta_v(R, M2, Gammas, zeta));
    fprintf('\n')
    fprintf('L3Drl          = %.6f\n', L3Drl(R, M2, Gammas));
    fprintf('L3Drs          = %.6f\n', L3Drs(R, M2, Gammas));
    fprintf('L3Dr           = %.6f\n', L3Dr(R, M2, Gammas));
    fprintf('L3Da           = %.6f\n', L3Da(R, M2, Gammas));
    fprintf('L3D            = %.6f\n', L3D(R, M2, Gammas));
    fprintf('\n')
    fprintf('T3Drl          = %.6f\n', T3Drl(R, M2, Gammas));
    fprintf('T3Drs          = %.6f\n', T3Drs(R, M2, Gammas));
    fprintf('T3Dr           = %.6f\n', T3Dr(R, M2, Gammas));
    fprintf('T3Da           = %.6f\n', T3Da(R, M2, Gammas));
    fprintf('T3D            = %.6f\n', T3D(R, M2, Gammas));
end

function value = P0(R)
    value = @(gamma) ((gamma + 1) - (gamma - 1) .* R.^(-1)) ./ ((gamma + 1).* R.^(-1) - (gamma - 1)); 
end

function value = compute_Gammas(M1, R, P)
    value =  7/5 * (M1(1:end-1).^2 ./ R(1:end-1).^2) .* ((P(2:end)  - P(1:end-1)) ./ (R(2:end)  - R(1:end-1))).^(-1);
end

function value = thetaOfzeta(R, M2, zeta)
    value = atan((M2 .* R) ./ (sqrt(1 - M2.^2) .* zeta(1:end)));
end

function value = sigma_a(R, M2, Gammas)
    value =  (R  ./ (R  - 1)) .* ((1 - Gammas) ./ (2*M2));
end

function value = sigma_b(M2, Gammas)
    value =  (1 + Gammas) ./ (2*M2);
end

function value = sigma_c(R, M2, Gammas)
    value =  (((M2.^2 .* (R - 1))) ./ (1 - M2.^2)) .* sigma_a(R, M2, Gammas);
end

function value = pi_l1(R, M2, Gammas, zeta)
    value = ((-(1 - R.^-1) .* (sigma_b(M2, Gammas) .* zeta.^2 - sigma_c(R, M2, Gammas))) ./ (zeta.^2.* (1 - zeta.^2) + (sigma_b(M2, Gammas) .* zeta.^2 - sigma_c(R, M2, Gammas)).^2)) .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2));
end

function value = pi_l2(R, M2, Gammas, zeta)
    value = (((1 - R.^-1) .* zeta .* sqrt((1 - zeta.^2))) ./ (zeta.^2 .* (1 - zeta.^2) + (sigma_b(M2, Gammas) .* zeta.^2 - sigma_c(R, M2, Gammas)).^2)) .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2));
end

function value = pi_s(R, M2, Gammas, zeta)
    value = ((-(1 - R.^-1)) ./ (zeta .* sqrt(zeta.^2 - 1) + sigma_b(M2, Gammas) .* zeta.^2 - sigma_c(R, M2, Gammas))) .* (zeta.^2 - (R .* M2.^2) ./ (1 - M2.^2));
end

function value = ka(M2, zeta)
    value = (zeta .* M2  - sqrt(zeta.^2 - 1)) ./ (sqrt(1 - M2.^2));
end

function value = wa(M2, zeta)
    value = (zeta - M2  .* sqrt(zeta.^2 - 1)) ./ (sqrt(1 - M2.^2));
end

function value = Delta_ua(R, M2, Gammas, zeta)
    value = (ka(M2, zeta) ./ wa(M2, zeta)) .* pi_s(R, M2, Gammas, zeta);
end

function value = Delta_va(R, M2, Gammas, zeta)
    value = (1 ./ wa(M2, zeta)) .* pi_s(R, M2, Gammas, zeta);
end

function value = Omega_1(R, M2, zeta)
    value = R .* (1 + zeta.^2 .* ((1 - M2.^2) ./ (R.^2 .* M2.^2)));
end

function value = Omega_2(R, M2, Gammas)
    value = ((R - 1) .* (1 - Gammas)) ./ (2*M2);
end

function value = Delta_Omega_l1(R, M2, Gammas, zeta)
    value = Omega_2(R, M2, Gammas) .* pi_l1(R, M2, Gammas, zeta) + Omega_1(R, M2, zeta);
end

function value = Delta_Omega_l2(R, M2, Gammas, zeta)
    value = Omega_2(R, M2, Gammas) .* pi_l2(R, M2, Gammas, zeta);
end

function value = Delta_Omega_s(R, M2, Gammas, zeta)
    value = Omega_2(R, M2, Gammas) .* pi_s(R, M2, Gammas, zeta) + Omega_1(R, M2, zeta);
end

function value = Delta(M2, zeta)
    value = 1 + ((1 - M2.^2) ./ M2.^2) .* zeta.^2;
end

function value = Delta_u_l1(R, M2, Gammas, zeta)
    value = Delta_Omega_l1(R, M2, Gammas, zeta) ./ Delta(M2, zeta);
end

function value = Delta_u_l2(R, M2, Gammas, zeta)
    value = Delta_Omega_l2(R, M2, Gammas, zeta) ./ Delta(M2, zeta);
end

function value = Delta_u_s(R, M2, Gammas, zeta)
    value = Delta_Omega_s(R, M2, Gammas, zeta) ./ Delta(M2, zeta);
end

function value = Delta_u(R, M2, Gammas, zeta)
    count = 0;   
    if zeta > 1
        value = Delta_u_s(R, M2, Gammas, zeta);
    elseif zeta == 1 && count < 1
        count = count + 1;
        value = Delta_u_s(R, M2, Gammas, zeta);
    elseif zeta == 1 && count < 2
        count = count + 1;
        value = sqrt((Delta_u_l1(R, M2, Gammas, zeta).^2 + Delta_u_l2(R, M2, Gammas, zeta).^2));
    else
        value = sqrt((Delta_u_l1(R, M2, Gammas, zeta).^2 + Delta_u_l2(R, M2, Gammas, zeta).^2));
    end
end

function value = Delta_v_l1(R, M2, Gammas, zeta)
    value = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_l1(R, M2, Gammas, zeta) ./ Delta(M2, zeta));
end

function value = Delta_v_l2(R, M2, Gammas, zeta)
    value = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_l2(R, M2, Gammas, zeta) ./ Delta(M2, zeta));
end

function value = Delta_v_s(R, M2, Gammas, zeta)
    value = zeta.* ((sqrt(1 - M2.^2)) ./ M2) .* (Delta_Omega_s(R, M2, Gammas, zeta) ./ Delta(M2, zeta));
end

function value = Delta_v(R, M2, Gammas, zeta)
    count = 0;   
    if zeta > 1
        value = Delta_v_s(R, M2, Gammas, zeta);
    elseif zeta == 1 && count < 1
        count = count + 1;
        value = Delta_v_s(R, M2, Gammas, zeta);
    elseif zeta == 1 && count < 2
        count = count + 1;
        value = sqrt((Delta_v_l1(R, M2, Gammas, zeta).^2 + Delta_v_l2(R, M2, Gammas, zeta).^2));
    else
        value = sqrt((Delta_v_l1(R, M2, Gammas, zeta).^2 + Delta_v_l2(R, M2, Gammas, zeta).^2));
    end
end

function value = pdf3D(R, M2, zeta)
    value = 3/2 * ((M2.^4 .* R.^4 .* sqrt(1 - M2.^2)) ./ (M2.^2 .* R.^2 + zeta.^2 .* (1 - M2.^2)).^(5/2));
end

function value = L3Drl(R, M2, Gammas)
    fun = @(zeta) (Delta_u_l1(R, M2, Gammas, zeta).^2 + Delta_u_l2(R, M2, Gammas, zeta).^2) .* pdf3D(R, M2, zeta); 
    value = integral(fun, 0, 1, 'Waypoints', linspace(0, 1, 2e3));
end

function value = L3Drs(R, M2, Gammas)
    fun = @(zeta) Delta_u_s(R, M2, Gammas, zeta).^2 .* pdf3D(R, M2, zeta);
    value = integral(fun, 1.005, Inf, 'Waypoints', [linspace(1.005, 100, 1e4), linspace(100.01, 1e5, 1e4)]);
end

function value = L3Dr(R, M2, Gammas)
    value = L3Drl(R, M2, Gammas) + L3Drs(R, M2, Gammas);
end

function value = L3Da(R, M2, Gammas)
    fun = @(zeta) Delta_ua(R, M2, Gammas, zeta).^2 .* pdf3D(R, M2, zeta);
    value = integral(fun, 1.005, Inf, 'Waypoints', [linspace(1.005, 100, 1e4), linspace(100.01, 1e5, 1e4)]);
end

function value = L3D(R, M2, Gammas)
    value = L3Dr(R, M2, Gammas) + L3Da(R, M2, Gammas);
end

function value = T3Drl(R, M2, Gammas)
    fun = @(zeta) (Delta_v_l1(R, M2, Gammas, zeta).^2 + Delta_v_l2(R, M2, Gammas, zeta).^2 + 3/2) .* pdf3D(R, M2, zeta);
    value =  0.5 * integral(fun, 0, 1, 'Waypoints', linspace(0, 1, 2e3));
end

function value = T3Drs(R, M2, Gammas)
    fun = @(zeta) (Delta_v_s(R, M2, Gammas, zeta).^2 + 3/2).* pdf3D(R, M2, zeta);
    value =   0.5 * integral(fun, 1.01, Inf, 'Waypoints', [linspace(1, 100, 1e4), linspace(100.01, 1e5, 1e5)]);
end

function value = T3Dr(R, M2, Gammas)
    value = T3Drl(R, M2, Gammas) + T3Drs(R, M2, Gammas);
end

function value = T3Da(R, M2, Gammas)
    fun = @(zeta) Delta_va(R, M2, Gammas, zeta).^2 .* pdf3D(R, M2, zeta);
    value =  0.5 * integral(fun, 1.002, Inf, 'Waypoints', [linspace(1.002, 100, 1e4), linspace(100.01, 1e10, 1e4)]);
end

function value = T3D(R, M2, Gammas)
    value = T3Dr(R, M2, Gammas) + T3Da(R, M2, Gammas);
end


function [R, P, T, M1, M2, Gammas, alpha] = get_parameters(self)
    R = cell2vector(self.PS.strP, 'rho') ./ cell2vector(self.PS.strR, 'rho');
    P = cell2vector(self.PS.strP, 'p') ./ cell2vector(self.PS.strR, 'p');
    T = cell2vector(self.PS.strP, 'T') ./ cell2vector(self.PS.strR, 'T');
    M1 = cell2vector(self.PS.strR, 'u') ./ cell2vector(self.PS.strR, 'sound');
    M2 = cell2vector(self.PS.strP, 'v_shock') ./ cell2vector(self.PS.strP, 'sound');
    Gammas = compute_Gammas(M1, R, P);

    Yi = cell2vector(self.PS.strP, 'Yi');
    alpha = Yi(2, :);

    R = R(1:end-1);
    P = P(1:end-1);
    T = T(1:end-1);
    M1 = M1(1:end-1);
    M2 = M2(1:end-1);
    alpha = alpha(1:end-1);
end

function self = compute_shock(varargin)
    % Read list of species considered as products
    if ~isempty(varargin{nargin})
        ListProducts = varargin{1}{1};
    else 
        ListProducts = 'Air';
    end
    if nargin > 2, T = varargin{1}{2}; else, T = 300; end % [K]
    if nargin > 3, p = varargin{1}{3}; else, p = 1; end % [bar]
    self = App(ListProducts);
    self.Misc.FLAG_RESULTS = false;
    % Initial conditions
    self = set_prop(self, 'TR', T, 'pR', 1.01325 * p);
    self.PD.S_Oxidizer = {'O2'};
    self.PD.S_Inert    = {'N2', 'Ar', 'CO2'};
    self.PD.proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;
    % Additional inputs
%     initial_velocity_sound = 3.529546069689621e+02; % N2
%     initial_velocity_sound = 329.4216; % O2
    initial_velocity_sound = 347.1035; % Air

%     u1 = linspace(initial_velocity_sound, 500, 200);
%     u1 = [u1, linspace(500.1, 3000, 1000)];
%     u1 = [u1, linspace(3000.1, 5000, 500)];
%     u1 = [u1, linspace(5000.1, 12000, 8000)];
%     u1 = [u1, linspace(12000.1, 15000, 500)];
    u1 = linspace(initial_velocity_sound,  10000, 1000);
%     u1 = logspace(2, 5, 2e3);
%     u1 = u1(u1 < 15000); u1 = u1(u1 >= initial_velocity_sound);

    self = set_prop(self, 'u1', u1, 'phi', self.PD.phi.value(1) * ones(1, length(u1)));
    % Solve problem
    self = SolveProblem(self, 'SHOCK_I');
end