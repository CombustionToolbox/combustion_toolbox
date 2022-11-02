function DhT = species_DhT_NASA(species, temperature, DB)
    % Compute thermal enthalpy [kJ/mol] of the species at the given
    % temperature [K] using NASA's 9 polynomials
    %
    % Args:
    %     species (str): Chemical species
    %     temperature (float): Range of temperatures to evaluate [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     DhT (float): Thermal enthalpy [kJ/mol]

    [~, DhT] = species_h0_NASA(species, temperature, DB);
end
