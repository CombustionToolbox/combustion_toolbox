function results = compute_LIA(R, M2, Gammas)
    % Routine to characterizes the turbulence amplification of a shock
    % interacting with weak turbulence using linear theory.
    %
    % These calculations are based on our previous theoretical work [1]
    % and have been extended to multi-component mixtures [2] using the
    % Combustion Toolbox [3].
    %
    % Note:
    %     The LIA results from [1] and [2] are computed using Mathematica
    %     and there can be some numerical differences with the results
    %     obtained with this MATLAB routine.
    %
    % References:
    %     [1] Huete, C., Cuadra, A., Vera, M., Urzay, & J. (2021). Thermochemical
    %         effects on hypersonic shock waves interacting with weak turbulence.
    %         Physics of Fluids 33, 086111 (featured article). DOI: 10.1063/5.0059948.
    %
    %     [2] Cuadra, A., Vera, M., Di Renzo, M., & Huete, C. (2023). Linear Theory
    %         of Hypersonic Shocks Interacting with Turbulence in Air. In 2023 AIAA
    %         SciTech Forum, National Harbor, USA. DOI: 10.2514/6.2023-0075.
    %
    %     [3] Cuadra, A., Huete, C., Vera, M., (2022). Combustion Toolbox:
    %         A MATLAB-GUI based open-source tool for solving gaseous
    %         combustion problems. Zenodo. DOI: 10.5281/zenodo.5554911.
    %
    % Args:
    %     R (float): Density ratio rho_2 / rho_1
    %     M2 (float): Post-shock Mach number
    %     Gammas (float): Inverse normalized Hugonitor slope, see Eq. (22) in [1] (Eq. (4) in [2])
    %
    % Returns:
    %     results (struct): Struct with the results
    
    % Definitions
    N = length(R);
    
    % Compute acoustic and vortical modes of the longitudinal and
    % transverse components of the turbulent kinetic energy (TKE)
    % amplification
    for i = N:-1:1
        results.L3Dr(i) = L3Dr(R(i), M2(i), Gammas(i));
        results.L3Da(i) = L3Da(R(i), M2(i), Gammas(i));
        results.T3Dr(i) = T3Dr(R(i), M2(i), Gammas(i));
        results.T3Da(i) = T3Da(R(i), M2(i), Gammas(i));
    end

    % Compute longitudinal contribution of the turbulent kinetic energy (TKE) amplification
    results.L3D = results.L3Dr + results.L3Da;

    % Compute transverse contribution of the turbulent kinetic energy (TKE) amplification
    results.T3D = results.T3Dr + results.T3Da;

    % Compute amplification of the turbulent kinetic energy (TKE)
    results.K3D = 1/3 * (results.L3D + 2*results.T3D);

    % Compute amplification of turbulence intensity
    results.I3D = results.R .* sqrt(results.K3D);

    % Compute amplification of turbulent Reynolds number
    results.RRe3D = (sqrt(results.K3D) .* sqrt((1 + 2 * results.R.^2) / 3)) ./ (results.T.^(0.7));

    % Compute anisotropy
    results.Anisotropy = 1 - (4 * results.L3D) ./ (3 * results.K3D + results.L3D);

end

% SUB PASS FUNCTIONS
function value = P0(R)
    value = @(gamma) ((gamma + 1) - (gamma - 1) .* R.^(-1)) ./ ((gamma + 1).* R.^(-1) - (gamma - 1)); 
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
    select = 0;
    if zeta > 1
        value = Delta_u_s(R, M2, Gammas, zeta);
    elseif zeta == 1 && select == 1
        value = Delta_u_s(R, M2, Gammas, zeta);
    elseif zeta == 1 && select == 2
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
    select = 0;   
    if zeta > 1
        value = Delta_v_s(R, M2, Gammas, zeta);
    elseif zeta == 1 && select == 1
        value = Delta_v_s(R, M2, Gammas, zeta);
    elseif zeta == 1 && select == 2
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
    value = integral(fun, 1, Inf, 'Waypoints', logspace(0, 5, 1e5));
end

function value = L3Dr(R, M2, Gammas)
    value = L3Drl(R, M2, Gammas) + L3Drs(R, M2, Gammas);
end

function value = L3Da(R, M2, Gammas)
    fun = @(zeta) Delta_ua(R, M2, Gammas, zeta).^2 .* pdf3D(R, M2, zeta);
    value = integral(fun, 1, Inf, 'Waypoints', logspace(0, 5, 1e5));
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
    value =   0.5 * integral(fun, 1, Inf, 'Waypoints', logspace(0, 5, 1e5));
end

function value = T3Dr(R, M2, Gammas)
    value = T3Drl(R, M2, Gammas) + T3Drs(R, M2, Gammas);
end

function value = T3Da(R, M2, Gammas)
    fun = @(zeta) Delta_va(R, M2, Gammas, zeta).^2 .* pdf3D(R, M2, zeta);
    value =  0.5 * integral(fun, 1, Inf, 'Waypoints', logspace(0, 5, 1e5));
end

function value = T3D(R, M2, Gammas)
    value = T3Dr(R, M2, Gammas) + T3Da(R, M2, Gammas);
end

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