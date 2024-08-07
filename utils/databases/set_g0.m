function g0 = set_g0(LS, T, DB)
    % Function that computes the vector of gibbs free energy for the given
    % set of species [J/mol]
    %
    % Args:
    %     LS (cell): List of species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     g0 (float): Gibbs energy in molar basis [J/mol]
    %
    % Example:
    %     g0 = set_g0({'H2O', 'CO2'}, 298.15, DB)

    for i = length(LS):-1:1
        species = LS{i};
        g0(i, 1) = species_g0(species, T, DB);
    end

end
