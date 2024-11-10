function s0 = set_s0(listSpecies, T, DB)
    % Function that computes the vector of entropy for the given
    % set of species [J/(mol-K)]
    %
    % Args:
    %     listSpecies (cell): List of species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     s0 (float): Entropy in molar basis [J/(mol-K)]
    %
    % Example:
    %     s0 = set_s0({'H2O', 'CO2'}, 298.15, DB)

    for i = length(listSpecies):-1:1
        species = listSpecies{i};
        s0(i, 1) = species_s0(species, T, DB);
    end

end
