function [R_j, P_j, T_j] = theo_diatomic_shocks(gas, T, T0, p1)
    % Get parameters
    [alpha, Br, theta_r, theta_v, theta_d] = get_parameters(gas, T, T0, p1);
    % Definitions
    R_j = ((1 - alpha) ./ alpha.^2) .* Br .* sqrt(T) .* (1 - exp(-theta_v ./ T)) .* exp(-theta_d ./ T);
    T_j = (6 - R_j.^(-1) - 2*alpha .* theta_d - 2.*(1 - alpha) * theta_v ./ (exp(theta_v./T) - 1)) ./ (2.*(alpha + 3) - R_j .* (1 + alpha));
    P_j = (1 + alpha) .* R_j .* T;
    
    alpha_R = -((alpha .* (1 - alpha))/(2 - alpha));
    alpha_T = -alpha_R .* (1/2 + (theta_d./T) .* ((1 - (1 + theta_v./theta_d) .* exp(-theta_v./T)) ./ (1 - exp(-theta_v./T))));
end

% NESTED FUNCTIONS
function [Tr, Tv, Td, G, m] = get_initial_state(gas)
    Tr = vpa(gas.Tr);
    Tv = vpa(gas.Tv);
    Td = vpa(gas.Td);
    G  = vpa(gas.G);
    m  = vpa(gas.m);
end

function [alpha, Br, theta_r, theta_v, theta_d] = get_parameters(gas, T_vector, T0, p0)
    syms x positive
    digits(250);
    alpha = zeros(1, length(T_vector));
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
    
    T_vector = vpa(T_vector);
    
    value = f_alpha(T_vector(end), [vpa(0), vpa(1)]);
    alpha(end) = value(1);
    for i = length(T_vector)-1:-1:1
        value = f_alpha(T_vector(i), [vpa(0), vpa(1)]);
        alpha(i) = value(1);
    end
    
    Br = double(Br);
    theta_r = double(theta_r);
    theta_v = double(theta_v);
    theta_d = double(theta_d);
end
