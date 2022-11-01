function cV = species_cV_NASA(species, temperature, DB)
    % Compute specific heat at constant volume [J/(mol-K)] of the species
    % at the given temperature [K] using NASA's 9 polynomials
    %
    % Args:
    %     species (str): Chemical species
    %     temperature (float): Range of temperatures to evaluate [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     cV (float): Specific heat at constant volume [J/(mol-K)]

    [~, cV] = species_cP_NASA(species, temperature, DB); % [J/(mol-K)];
end
