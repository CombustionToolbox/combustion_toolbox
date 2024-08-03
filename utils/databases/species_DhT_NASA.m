function DhT = species_DhT_NASA(species, temperature, DB)
    % Compute thermal enthalpy [J/mol] of the species at the given
    % temperature [K] using NASA's 9 polynomials
    %
    % Args:
    %     species (char): Chemical species
    %     temperature (float): Range of temperatures to evaluate [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     DhT (float): Thermal enthalpy in molar basis [J/mol]
    %
    % Example:
    %     DhT = species_DhT_NASA('H2O', 300:100:6000, DB)

    [~, DhT] = species_h0_NASA(species, temperature, DB);
end
