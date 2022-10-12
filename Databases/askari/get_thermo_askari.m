function cP = get_thermo_askari(LS, T, p, DB)
    % Function that computes the vector of specific heats at constant
    % pressure for the given set of species [J/(mol-K)]
    %
    % Args:
    %     LS (cell): List of species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     cP (float): Specific heat at constant pressure [J/(mol-K)]

    for i=length(LS):-1:1
        species = LS{i};
        coeff = compute_coefficents(species, p, DB);
        W(i) = compute_W(T, coeff); % [kg/mol]
        cP(i) = compute_cp(T, coeff) * W; % [J/(mol-K)]
    end
end

% SUB-PASS FUNCTIONS
function W = compute_W(T, coeff)
    W = coeff.lambda(4) - sum(coeff.a_m ./ (1 + (coeff.b_m / T).^coeff.c_m)) * 1e-3; % [kg/mol]
end

function [cp, cp_r, cp_f] = compute_cp(T, coeff)
    cp_f = coeff.lambda(1) + sum(coeff.a_f ./ (1 + (coeff.b_f / T).^coeff.c_f)); % [kJ/kg-K]
    cp_r = sum(coeff.a_r .* exp(-(log(T ./ coeff.b_r) / coeff.c_r).^2)); % [kJ/kg-K]
    cp_f = cp_f * 1e3;  % [J/kg-K]
    cp_r = cp_r * 1e3;  % [J/kg-K]
    cp   = cp_f + cp_r; % [J/kg-K]
end

function polynomial = compute_polynomial(species, p, DB, variable)
    log_p_pow_j = log(p).^(0:1:DB.(species).m);
    for i = length(DB.(species).(variable)(:, 2)):-1:1
        polynomial(i) = dot(DB.(species).(variable)(i, :), log_p_pow_j);
    end
end

function coeff = compute_coefficents(species, p, DB)
    % Frozen coefficients [kJ/kg-K]
    coeff.a_f = exp(compute_polynomial(species, p, DB, 'a_f'));
    coeff.b_f = exp(compute_polynomial(species, p, DB, 'b_f'));
    coeff.c_f = exp(compute_polynomial(species, p, DB, 'c_f'));
    % Reaction coefficients [kJ/kg-K]
    coeff.a_r = exp(compute_polynomial(species, p, DB, 'a_r'));
    coeff.b_r = exp(compute_polynomial(species, p, DB, 'b_r'));
    coeff.c_r = exp(compute_polynomial(species, p, DB, 'c_r'));
    % Mean molar mass coefficients [kg/kgmol]
    coeff.a_m = exp(compute_polynomial(species, p, DB, 'a_m'));
    coeff.b_m = exp(compute_polynomial(species, p, DB, 'b_m'));
    coeff.c_m = exp(compute_polynomial(species, p, DB, 'c_m'));
    % lambda
    coeff.lambda = compute_polynomial(species, p, DB, 'lambda');
end