function Gammas = compute_Gammas_frozen(M1, R, P)
    % Compute slope of the Hugoniot curve for thermochemically frozen air
    %
    % Args:
    %     M1 (float): Pre-shock Mach number [-]
    %     R (float): Density ratio [-]
    %     P (float): Pressure ratio [-]
    %
    % Returns:
    %     Gammas (float): Slope of the Hugoniot curve [-]

    Gammas =  7/5 * (M1(1:end-1).^2 ./ R(1:end-1).^2) .* compute_first_derivative(P, R).^(-1);
end