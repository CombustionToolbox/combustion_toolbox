function [R, P, T] = shock_ideal_gas(gamma, M1)
    % Compute jump conditions assuming a thermochemically frozen gas
    %
    % Args:
    %     gamma (float): Adiabatic index [-]
    %     M1 (float):    Pre-shock Mach number [-]
    %
    % Returns:
    %     Tuple containing
    %
    %     - R (float):     Density ratio [-]
    %     - P (float):     Pressure ratio [-]
    %     - T (float):     Temperature ratio [-]

    R = (gamma + 1) .* M1.^2 ./ ((gamma - 1) .* M1.^2 + 2);
    P = (2 * gamma * M1.^2 - (gamma - 1)) ./ (gamma + 1);
    T = (((gamma - 1) .* M1.^2 + 2) .* (2 * gamma * M1.^2 - (gamma - 1))) ./ ((gamma + 1).^2 .* M1.^2);
end