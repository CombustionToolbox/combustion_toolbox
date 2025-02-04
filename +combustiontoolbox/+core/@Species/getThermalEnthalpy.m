function DhT = getThermalEnthalpy(obj, T)
    % Compute thermal enthalpy [J/mol] of the species at the given
    % temperature [K] using piecewise cubic Hermite interpolating
    % polynomials and linear extrapolation
    %
    % Args:
    %     obj (Species): Species object
    %     T (float): Temperature [K]
    %
    % Returns:
    %     DhT (float): Thermal enthalpy in molar basis [J/mol]
    %
    % Example:
    %     DhT = getThermalEnthalpy(obj, 300)

    % Compute enthalpy [J/mol]
    h0 = getEnthalpy(obj, T);

    % Compute thermal enthalpy [J/mol]
    DhT = h0 - obj.hf;
end
