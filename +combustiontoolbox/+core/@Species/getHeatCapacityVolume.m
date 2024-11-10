function cv = getHeatCapacityVolume(obj, T)
    % Compute specific heat at constant volume [J/(mol-K)] of the species
    % at the given temperature [K] using piecewise cubic Hermite
    % interpolating polynomials and linear extrapolation
    %
    % Args:
    %     obj (Species): Species object
    %     T (float): Temperature [K]
    %
    % Returns:
    %     cv (float): Specific heat at constant volume in molar basis [J/(mol-K)]
    %
    % Example:
    %     cv = getHeatCapacityVolume(obj, 300)

    % Universal gas constant [J/(mol-K)]
    R0 = combustiontoolbox.common.Constants.R0;

    % Compute specific heat at constant volume [J/(mol-K)]
    cv = getHeatCapacityPressure(obj, T) - R0;
end
