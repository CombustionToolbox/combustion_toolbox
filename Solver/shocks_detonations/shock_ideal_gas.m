function [R, P, T] = shock_ideal_gas(gamma, M1)
    % Compute jump conditions for a given adiabatic index (gamma) and
    % a array of incident Mach numbers (M1).
    R = (gamma + 1) .* M1.^2 ./ ((gamma - 1) .* M1.^2 + 2);
    P = (2 * gamma * M1.^2 - (gamma - 1)) ./ (gamma + 1);
    T = (((gamma - 1) .* M1.^2 + 2) .* (2 * gamma * M1.^2 - (gamma - 1))) ./ ((gamma + 1).^2 .* M1.^2);
end