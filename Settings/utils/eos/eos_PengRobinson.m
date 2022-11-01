function [V, a_mix, b_mix, Vi, a, b] = eos_PengRobinson(self, T, p, species, Xi)
    % Compute molar volume of the mixture considering Peng-Robinson
    % Equation of State (EoS)
    % 
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     species (cell): List of the species in the mixture
    %     Xi (float): molar fractions of the mixture
    %     T (float): Temperature [K]
    %     p (float): Pressure [bar]
    % 
    % Returns:
    %     Tuple containing
    %
    %     - V (float): molar volume of the mixture [m3/mol]
    %     - Vi (float): molar volume of the components [m3/mol]

    % Definitions
    R0 = self.C.R0; % Universal gas constant [J/mol-K]
    T_c = set_prop_DB(species, 'critical_temperature', self.DB); % [K]
    p_c = set_prop_DB(species, 'critical_pressure', self.DB); % [Pa]
    w = set_prop_DB(species, 'acentric_factor', self.DB); % [-]
    % Change of units
    p = convert_bar_to_Pa(p); % [Pa]
    % Compute individual coefficients aij and bj
    [a, b] = compute_coefficients_ab(T, T_c, p_c, w, R0); 
    % Compute mixture coefficients
    [a_mix, b_mix] = apply_mixing_rule(Xi, a, b);
    % Compute guess molar volume assuming ideal EoS
    V_guess = compute_v_molar_guess(T, p, R0);
    % Compute molar volume assuming non ideal EoS
    Vi = compute_v_molar_roots_components(self, species, a, b, T, p, R0);
    V = compute_v_molar_roots(a_mix, b_mix, T, p, R0);
    [V_newton, STOP, it] = compute_v_molar_newton(self, V_guess, a_mix, b_mix, T, p, R0);
end

% SUB-PASS FUNCTIONS
function [a, b] = compute_coefficients_ab(T, T_c, p_c, w, R0)
    % Compute individual coefficients (a, b)
    alpha = get_alpha(T, T_c, w);
    a = 0.45724 * alpha .* R0^2 .* T_c.^2 ./ p_c;
    b = 0.0778 * R0 * T_c ./ p_c;

    a(isnan(a)) = 0;
    b(isnan(b)) = 0;
end

function [a_mix, b_mix] = apply_mixing_rule(Xi, a, b)
    % Compute coefficients of the mixture applying the mixing rule ...
    kij = 0; % Binary Interaction Parameter (BIP)
    aij = sqrt(kron(a, a')) .* (1 - kij);
    a_mix = dot(Xi' * aij, Xi);
    b_mix = dot(Xi, b);
end

function value = get_alpha(T, T_c, w)
    % Compute parameter alpha
    T_r = T ./ T_c; % Reduced temperature
    beta = 0.37464 + 1.54226 * w - 0.26992 * w.^2;
    value = (1 + beta .* (1 - sqrt(T_r))).^2;
end

function V_guess = compute_v_molar_guess(T, p, R)
    % Compute molar volume considering ideal EoS
    V_guess = R * T / p;
end

function value = get_f(v, a, b, T, p, R0)
    % Function f(v) = 0
    value = v^3 + (b - R0 * T / p) * v^2 - (3*b^2 + 2*R0 * T * b / p - a / p) * v + (b^3 + R0 * T * b^2 - a * b) / p;
end

function value = get_fprime(v, a, b, T, p, R0)
    % Derivative fprime = df(v)
    value = 3*v^2 + 2*(b - R0 * T / p) * v - (3*b^2 + 2*R0 * T * b / p - a / p);
end

function [f, fprime] = compute_newton_functions(v, a, b, T, p, R0)
    % Evaluate function f and its derivative fprime
    f = get_f(v, a, b, T, p, R0);
    fprime = get_fprime(v, a, b, T, p, R0);
end

function [x, STOP, it] = compute_v_molar_newton(self, x0, a, b, T, p, R0)
    % Obtain molar volume [m3/mol] using Newton-Raphson method
    it = 0; STOP = 1.0;
    while STOP > self.TN.tol_eos && it < self.TN.it_eos
        it = it + 1;
        
        [f0, fprime0] = compute_newton_functions(x0, a, b, T, p, R0);

        x = abs(x0 - f0 / fprime0);

        STOP = abs((x - x0) / x);
        x0 = x;
    end
end

function V = compute_v_molar_roots(a, b, T, p, R0)
    % Obtain molar volume [m3/mol] solving the cubic polynomial
    coefficients = [1, (b - R0 * T / p), - (3*b^2 + 2*R0 * T * b / p - a / p), (b^3 + R0 * T * b^2 - a * b) / p];
    r = roots(coefficients);
    r = r(imag(r) == 0);
    V = max(r);
end

function Vi = compute_v_molar_roots_components(self, species, a, b, T, p, R0)
    % Obtain molar volume [m3/mol] solving the cubic polynomial
    for i = length(species):-1:1
        r = compute_v_molar_roots(a(i), b(i), T, p, R0);
        if self.DB.(species{i}).phase
            Vi(i) = min(r);
        else
            Vi(i) = max(r);
        end
    end
end