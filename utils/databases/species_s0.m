function s0 = species_s0(species, T, DB)
    % Compute entropy [J/(mol-K)] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     s0 (float): Entropy in molar basis [J/(mol-K)]
    %
    % Example:
    %     s0 = species_s0('H2O', 300, DB)

    persistent cachedSpecies;
    persistent cachedS0curves;
    
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedS0curves = {};
    end
    
    % Check if species data is already cached
    index = find(strcmp(cachedSpecies, species), 1);
    if isempty(index)
        % Load species data and cache it
        s0curve = DB.(species).s0curve;
        cachedSpecies{end+1} = species;
        cachedS0curves{end+1} = s0curve;
    else
        % Retrieve cached data
        s0curve = cachedS0curves{index};
    end
    
    % Compute entropy [J/(mol-K)]
    s0 = s0curve(T);
end
