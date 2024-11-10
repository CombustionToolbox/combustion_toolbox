function DhT = set_DhT(listSpecies, T, DB)
    % Function that computes the vector of thermal enthalpy for the given
    % set of species [J/mol]
    %
    % Args:
    %     listSpecies (cell): List of species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     DhT (float): Thermal enthalpy in molar basis [J/mol]
    %
    % Example:
    %     DhT = set_DhT({'H2O', 'CO2'}, 298.15, DB)

    for i = length(listSpecies):-1:1
        species = listSpecies{i};
        DhT(i, 1) = species_DhT(species, T, DB);
    end

end
