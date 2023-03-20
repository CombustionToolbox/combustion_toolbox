function e0 = species_e0(species, T, DB)
    % Compute internal energy [kJ/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (str): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     e0 (float): Internal energy [kJ/mol]
    
    % Specific value
    moles = 1;
    % Universal gas constant
    R0 = 8.31446261815324; % [J/(K mol)]. 
    % Enthalpy [J/mol]
    h0 = species_h0(species, T, DB) * 1e3;
    % Internal energy [kJ/mol]
    e0 = (h0 - moles * R0 * T)/ 1000; % [kJ/mol]
end