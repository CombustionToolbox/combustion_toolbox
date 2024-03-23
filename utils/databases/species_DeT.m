function DeT = species_DeT(species, T, DB)
    % Compute thermal internal energy [kJ/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     DeT (float): Thermal internal energy in molar basis [kJ/mol]
    %
    % Example:
    %     DeT = species_DeT('H2O', 300, DB)

    persistent cachedSpecies;
    persistent cachedDETcurves;
    
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedDETcurves = {};
    end
    
    % Check if species data is already cached
    idx = find(strcmp(cachedSpecies, species), 1);
    if isempty(idx)
        % Load species data and cache it
        DeTcurve = DB.(species).DeTcurve;
        cachedSpecies{end+1} = species;
        cachedDETcurves{end+1} = DeTcurve;
    else
        % Retrieve cached data
        DeTcurve = cachedDETcurves{idx};
    end
    
    % Compute thermal internal energy [kJ/mol]
    DeT = DeTcurve(T) / 1000;
end