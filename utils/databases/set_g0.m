function g0 = set_g0(listSpecies, T, DB)
    % Function that computes the vector of gibbs free energy for the given
    % set of species [J/mol]
    %
    % Args:
    %     listSpecies (cell): List of species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     g0 (float): Gibbs energy in molar basis [J/mol]
    %
    % Example:
    %     g0 = set_g0({'H2O', 'CO2'}, 298.15, DB)

    for i = length(listSpecies):-1:1
        species = listSpecies{i};
        g0(i, 1) = getGibbsEnergy(DB.(species), T); % [J/mol]
    end

end
