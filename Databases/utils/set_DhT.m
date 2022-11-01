function DhT = set_DhT(LS, T, DB)
    % Function that computes the vector of thermal enthalpy for the given
    % set of species [J/mol]
    %
    % Args:
    %     LS (cell): List of species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     DhT (float): Thermal enthalpy [J/mol]

    for i = length(LS):-1:1
        species = LS{i};
        DhT(i, 1) = species_DhT(species, T, DB) * 1e3;
    end

end
