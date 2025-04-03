function DeT = getThermalInternalEnergy(obj, T)
    % Compute thermal internal energy [J/mol] of the species at the given
    % temperature [K] using piecewise cubic Hermite interpolating
    % polynomials and linear extrapolation
    %
    % Args:
    %     obj (Species): Species object
    %     T (float): Temperature [K]
    %
    % Returns:
    %     DeT (float): Thermal internal energy in molar basis [J/mol]
    %
    % Example:
    %     DeT = getThermalInternalEnergy(obj, 300)

    % Compute internal energy [J/mol]
    e0 = getInternalEnergy(obj, T);

    % Compute thermal internal energy [J/mol]
    DeT = e0 - obj.ef;
end