function h0 = set_h0(LS, T, DB)
    % Function that computes the vector of enthalpies for the given set of
    % species [J/mol]
    %
    % Args:
    %     LS (cell): List of species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     h0 (float): Enthalpy in molar basis [J/mol]
    %
    % Example:
    %     h0 = set_h0({'H2O', 'CO2'}, 298.15, DB)

    for i = length(LS):-1:1
        species = LS{i};
        h0(i, 1) = species_h0(species, T, DB); % [J/mol]
    end

end
