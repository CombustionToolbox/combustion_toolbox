function h0 = species_h0(species, T, DB)
    % Compute enthalpy [J/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     h0 (float): Enthalpy in molar basis [J/mol]
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
    index = find(strcmp(cachedSpecies, species), 1);
    if isempty(index)
        % Load species data and cache it
        h0curve = DB.(species).h0curve;
        cachedSpecies{end+1} = species;
        cachedH0curves{end+1} = h0curve;
    else
        % Retrieve cached data
        h0curve = cachedH0curves{index};
    end
    
    % Compute enthalpy [J/mol]
    h0 = h0curve(T);
end
