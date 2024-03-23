function h0 = species_h0(species, T, DB)
    % Compute enthalpy [kJ/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     h0 (float): Enthalpy in molar basis [kJ/mol]
    %
    % Example:
    %     h0 = species_h0('H2O', 300, DB)

    persistent cachedSpecies;
    persistent cachedH0curves;
    
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedH0curves = {};
    end
    
    % Check if species data is already cached
    idx = find(strcmp(cachedSpecies, species), 1);
    if isempty(idx)
        % Load species data and cache it
        h0curve = DB.(species).h0curve;
        cachedSpecies{end+1} = species;
        cachedH0curves{end+1} = h0curve;
    else
        % Retrieve cached data
        h0curve = cachedH0curves{idx};
    end
    
    % Compute enthalpy [kJ/mol]
    h0 = h0curve(T) / 1000;
end
