function DeT = species_DeT_NASA(species, temperature, DB)
    % Compute thermal internal energy [J/mol] of the species at the given
    % temperature [K] using NASA's 9 polynomials
    %
    % Args:
    %     species (char): Chemical species
    %     temperature (float): Range of temperatures to evaluate [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     DeT (float): Thermal internal energy in molar basis [J/mol]
    %
    % Example:
    %     DeT = species_DeT_NASA('H2O', 300:100:6000, DB)

    [~, DeT] = species_e0_NASA(species, temperature, DB);
end
