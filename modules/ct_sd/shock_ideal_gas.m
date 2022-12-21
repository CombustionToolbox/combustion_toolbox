function [R, P, T, Gammas, M1] = shock_ideal_gas(gamma, M1)
    % Compute jump conditions assuming a thermochemically frozen gas
    %
    % Args:
    %     gamma (float): Adiabatic index [-]
    %     M1 (float):    Pre-shock Mach number [-]
    %
    % Returns:
    %     Tuple containing
    %
    %     - R (float):      Density ratio [-]
    %     - P (float):      Pressure ratio [-]
    %     - T (float):      Temperature ratio [-]
    %     - Gammas (float): Rankine-Hugoniot slope parameter [-]

    R = (gamma + 1) .* M1.^2 ./ ((gamma - 1) .* M1.^2 + 2);
    P = (2 * gamma * M1.^2 - (gamma - 1)) ./ (gamma + 1);
    T = (((gamma - 1) .* M1.^2 + 2) .* (2 * gamma * M1.^2 - (gamma - 1))) ./ ((gamma + 1).^2 .* M1.^2);

    if length(M1) > 1
        % Compute Gamma
        % Gammas = compute_Gammas_frozen(M1, R, P);
        Gammas = compute_Gammas(u2, rho2, p2);
        % Remove last point of the dataset to match dimensions with Gammas
        R = R(1:end-1);
        P = P(1:end-1);
        T = T(1:end-1);
        M1 = M1(1:end-1);
    else
        Gammas = [];
    end
end
