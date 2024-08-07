function [R, P, T, Gammas, M1] = shockIncidentIdeal(obj, gamma, M1)
    % Compute jump conditions assuming a thermochemically frozen gas (calorically perfect gas)
    %
    % Args:
    %     obj (ShockSolver): ShockSolver object
    %     gamma (float): Adiabatic index [-]
    %     M1 (float): Pre-shock Mach number [-]
    %
    % Returns:
    %     Tuple containing
    %
    %     * R (float): Density ratio [-]
    %     * P (float): Pressure ratio [-]
    %     * T (float): Temperature ratio [-]
    %     * Gammas (float): Rankine-Hugoniot slope parameter [-]
    %     * M1 (float): Pre-shock Mach number [-]
    %
    % Example:
    %     [R, P, T, Gammas, M1] = shockIncidentIdeal(ShockSolver(), 1.4, 2.0)

    R = (gamma + 1) .* M1.^2 ./ ((gamma - 1) .* M1.^2 + 2);
    P = (2 * gamma * M1.^2 - (gamma - 1)) ./ (gamma + 1);
    T = (((gamma - 1) .* M1.^2 + 2) .* (2 * gamma * M1.^2 - (gamma - 1))) ./ ((gamma + 1).^2 .* M1.^2);

    if isscalar(M1)
        Gammas = [];
        return
    end

    % Compute Gamma numerically
    Gammas = obj.computeGammasFrozen(M1, R, P);

    % Remove last point of the dataset to match dimensions with Gammas
    R = R(1:end-1);
    P = P(1:end-1);
    T = T(1:end-1);
    M1 = M1(1:end-1);
end
