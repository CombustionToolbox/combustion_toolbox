function e0 = getInternalEnergy(obj, T)
    % Compute internal energy [J/mol] of the species at the given
    % temperature [K] using piecewise cubic Hermite interpolating
    % polynomials and linear extrapolation
    %
    % Args:
    %     obj (Species): Species object
    %     T (float): Temperature [K]
    %
    % Returns:
    %     e0 (float): Internal energy in molar basis [J/mol]
    %
    % Example:
    %     e0 = getInternalEnergy(obj, 300)
    
    % Universal gas constant [J/(K mol)]
    R0 = combustiontoolbox.common.Constants.R0;

    % Enthalpy [J/mol]
    h0 = getEnthalpy(obj, T);
    
    % Internal energy [J/mol]
    e0 = h0 - R0 * T;
end