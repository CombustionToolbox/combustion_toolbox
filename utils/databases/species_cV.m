function cv = species_cV(species, T, DB)
    % Compute specific heat at constant volume [J/(mol-K)] of the species
    % at the given temperature [K] using piecewise cubic Hermite
    % interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     cv (float): Specific heat at constant volume in molar basis [J/(mol-K)]
    %
    % Example:
    %     cv = species_cv('H2O', 300, DB)

    % Universal gas constant [J/(mol-K)]
    R0 = 8.31446261815324;

    % Compute specific heat at constant volume [J/(mol-K)]
    cv = species_cP(species, T, DB) - R0;
end
