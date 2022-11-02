function cP = set_cP(LS, T, DB)
    % Function that computes the vector of specific heats at constant
    % pressure for the given set of species [J/(mol-K)]
    %
    % Args:
    %     LS (cell): List of species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     cP (float): Specific heat at constant pressure [J/(mol-K)]

    for i = length(LS):-1:1
        species = LS{i};
        cP(i, 1) = species_cP(species, T, DB); % [J/(mol-K)]
    end

end
