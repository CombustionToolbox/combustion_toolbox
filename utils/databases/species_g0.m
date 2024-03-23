function g0 = species_g0(species, T, DB)
    % Compute Gibbs energy [kJ/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     g0 (float): Gibbs energy in molar basis [kJ/mol]
    %
    % Example:
    %     g0 = species_g0('H2O', 298.15, DB)
    
    persistent cachedSpecies;
    persistent cachedG0curves;
    
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedG0curves = {};
    end
    
    % Check if species data is already cached
    idx = find(strcmp(cachedSpecies, species), 1);
    if isempty(idx)
        % Load species data and cache it
        g0curve = DB.(species).g0curve;
        cachedSpecies{end+1} = species;
        cachedG0curves{end+1} = g0curve;
    else
        % Retrieve cached data
        g0curve = cachedG0curves{idx};
    end
    
    % Compute Gibbs energy [kJ/mol]
    g0 = g0curve(T) / 1000;
end

