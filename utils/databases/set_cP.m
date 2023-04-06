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
    %     cP (float): Specific heat at constant pressure in molar basis [J/(mol-K)]
    %
    % Example:
    %     cP = set_cP({'H2O', 'CO2'}, 298.15, DB)

    for i = length(LS):-1:1
        species = LS{i};
        cP(i, 1) = species_cP(species, T, DB); % [J/(mol-K)]
    end

end
