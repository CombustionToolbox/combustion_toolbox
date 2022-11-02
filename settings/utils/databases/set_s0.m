function s0 = set_s0(LS, T, DB)
    % Function that computes the vector of entropy for the given
    % set of species [J/(mol-K)]
    %
    % Args:
    %     LS (cell): List of species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     s0 (float): Entropy [J/(mol-K)]

    for i = length(LS):-1:1
        species = LS{i};
        s0(i, 1) = species_s0(species, T, DB) * 1e3;
    end

end
