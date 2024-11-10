function gamma = getAdiabaticIndex(obj, T)
    % Compute adiabatic index of the species [-] at the given temperature
    % [K] using piecewise cubic Hermite interpolating polynomials and
    % linear extrapolation
    %
    % Args:
    %     obj (Species): Species object
    %     T (float): Temperature [K]
    %
    % Returns:
    %     gamma (float): Adiabatic index [-]
    %
    % Example:
    %     gamma = getAdiabaticIndex(obj, 300)

    gamma = getHeatCapacityPressure(obj, T) ./ getHeatCapacityVolume(obj, T);

    assert(any(~isnan(gamma)), 'Adibatic index equal NaN');
end