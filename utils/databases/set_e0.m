function e0 = set_e0(LS, T, DB)
    % Function that computes the vector of internal energy for the given
    % set of species [J/mol]
    %
    % Args:
    %     LS (cell): List of species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     e0 (float): Internal energy [J/mol]

    for i = length(LS):-1:1
        species = LS{i};
        e0(i, 1) = species_DeT(species, T, DB) * 1e3 - DB.(species).ef;
    end

end
