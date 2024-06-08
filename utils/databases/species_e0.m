function e0 = species_e0(species, T, DB)
    % Compute internal energy [J/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     e0 (float): Internal energy in molar basis [J/mol]
    %
    % Example:
    %     e0 = species_e0('H2O', 300, DB)
    
    % Universal gas constant [J/(K mol)]
    R0 = 8.31446261815324;

    % Enthalpy [J/mol]
    h0 = species_h0(species, T, DB);
    
    % Internal energy [J/mol]
    e0 = h0 - R0 * T;
end